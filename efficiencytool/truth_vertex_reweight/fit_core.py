"""
fit_core — shared data-loading and histogram-building for the truth
vertex reweight fit.

Event-level formulation: each MC event contributes directly to the
reco-vertex histogram via a per-event weight that includes the
cross-section scaling, the parametric reweight w(z) evaluated at the
truth vertex (or vertices), and optionally the stage-2 residual
correction c(z).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import uproot


SENTINEL_F = -9999.0
SENTINEL_I = -9999

LOG = logging.getLogger("fit_core")


# ---------------------------------------------------------------------------
# Binning / specs / loaded samples
# ---------------------------------------------------------------------------

@dataclass
class Binning:
    z_min: float
    z_max: float
    n_bins: int
    overflow_clip: bool = False

    @property
    def edges(self) -> np.ndarray:
        return np.linspace(self.z_min, self.z_max, self.n_bins + 1)

    @property
    def centers(self) -> np.ndarray:
        e = self.edges
        return 0.5 * (e[:-1] + e[1:])


@dataclass
class SampleSpec:
    name: str
    path: str
    xs: float
    n_gen: float
    kind: str                # "single_mc", "double_mc", or "data"

    @property
    def event_weight(self) -> float:
        """Per-event cross-section weight = xs / n_gen. For data this is
        unused (data events are not weighted)."""
        if self.n_gen <= 0:
            return 1.0
        return float(self.xs) / float(self.n_gen)


@dataclass
class LoadedSample:
    spec: SampleSpec
    z_h: Optional[np.ndarray]             # vertexz_truth
    z_mb: Optional[np.ndarray]            # vertexz_truth_mb (double MC only)
    z_r: np.ndarray                       # vertexz (reco)
    runnumber: Optional[np.ndarray]       # runnumber (data only, typically)
    trigger_passed: Optional[np.ndarray] = None   # bool[n_events]; data-only.
    # Pre-filtered against the analysis trigger bit (e.g. bit 30 = Photon_4_GeV).
    # Mirrors the pipeline's `scaledtrigger[itrig] != 0` cut at line 1525 of
    # RecoEffCalculator_TTreeReader.C; the fit sees the same data subset the
    # downstream MC is compared against, so reweight closure transfers cleanly.


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

def _resolve_paths(path_spec) -> List[str]:
    """Accept a string, a glob pattern, or a list of strings. Return the
    expanded, sorted concrete file list."""
    import glob as _glob
    if isinstance(path_spec, (list, tuple)):
        out: List[str] = []
        for p in path_spec:
            out.extend(_resolve_paths(p))
        return out
    p = str(path_spec)
    if any(c in p for c in "*?[]"):
        matches = sorted(_glob.glob(p))
        if not matches:
            raise FileNotFoundError(f"no files matched pattern {p!r}")
        return matches
    return [p]


def _max_mtime(paths: Iterable[str]) -> float:
    mt = 0.0
    for p in paths:
        try:
            mt = max(mt, Path(p).stat().st_mtime)
        except FileNotFoundError:
            pass
    return mt


def load_sample(
    spec: SampleSpec,
    branches: Dict[str, str],
    n_events: Optional[int] = None,
    require_mb: bool = False,
    trigger_bit: Optional[int] = None,
) -> LoadedSample:
    """Read one or more slimtree files and return raw vertex arrays.
    ``spec.path`` may be a single path, glob, or list. Missing optional
    branches (vertexz_truth_mb, runnumber) are silently skipped.

    If ``trigger_bit`` is given (and ``scaledtrigger`` is present in the tree),
    a per-event bool array ``trigger_passed = scaledtrigger[:, trigger_bit] != 0``
    is attached to the LoadedSample. Used to mirror the analysis pipeline's
    trigger filter on the data side.
    """
    tree_name = branches["tree_name"]
    bz_h = branches["vertexz_truth"]
    bz_mb = branches["vertexz_truth_mb"]
    bz_r = branches["vertexz_reco"]
    brun = branches.get("runnumber")

    paths = _resolve_paths(spec.path)
    with uproot.open(paths[0]) as f:
        if tree_name not in f:
            raise KeyError(f"tree {tree_name!r} not found in {paths[0]}")
        keys = set(f[tree_name].keys())
    want: List[str] = [bz_r]
    have_zh = bz_h in keys
    have_mb = bz_mb in keys
    have_run = bool(brun and brun in keys)
    have_trig = trigger_bit is not None and "scaledtrigger" in keys
    if have_zh:
        want.append(bz_h)
    if have_mb:
        want.append(bz_mb)
    if have_run:
        want.append(brun)
    if have_trig:
        want.append("scaledtrigger")
    if bz_r not in keys:
        raise KeyError(f"required branch {bz_r!r} missing in {paths[0]}")

    tree_specs = [f"{p}:{tree_name}" for p in paths]
    arrs = uproot.concatenate(
        tree_specs, expressions=want, library="np", entry_stop=n_events,
    )

    z_r = arrs[bz_r].astype(np.float64)
    z_h = arrs[bz_h].astype(np.float64) if have_zh else None
    z_mb = arrs[bz_mb].astype(np.float64) if have_mb else None
    run = arrs[brun].astype(np.int64) if have_run else None
    trig_pass = None
    if have_trig:
        st = arrs["scaledtrigger"]
        trig_pass = (st[:, int(trigger_bit)] != 0).astype(np.bool_)
        LOG.info("trigger filter bit %d: %d / %d events pass (%.1f%%)",
                 int(trigger_bit), int(trig_pass.sum()), len(trig_pass),
                 100.0 * float(trig_pass.sum()) / max(len(trig_pass), 1))

    if require_mb and (z_mb is None or np.all(z_mb <= SENTINEL_F + 1)):
        raise ValueError(
            f"sample {spec.name} at {spec.path} has no usable "
            f"vertexz_truth_mb (require_mb=True)."
        )

    return LoadedSample(
        spec=spec, z_h=z_h, z_mb=z_mb, z_r=z_r, runnumber=run,
        trigger_passed=trig_pass,
    )


def load_sample_cached(
    spec: SampleSpec,
    branches: Dict[str, str],
    cache_dir: Path,
    n_events: Optional[int] = None,
    require_mb: bool = False,
    trigger_bit: Optional[int] = None,
) -> LoadedSample:
    """Same as load_sample, but save/restore event arrays to a per-sample
    .npz file under ``cache_dir``. The cache is considered valid iff its
    mtime is >= the newest source-file mtime AND the cached array
    length matches n_events (for n_events=None, any length is accepted).
    Single source-of-truth for fast iteration.

    If ``trigger_bit`` is given, the cache key includes the bit so different
    trigger filters get separate npz files."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    tag = "all" if n_events is None else f"n{n_events}"
    if trigger_bit is not None:
        tag = f"{tag}_t{int(trigger_bit)}"
    cache_path = cache_dir / f"{spec.name}__{tag}.npz"

    paths = _resolve_paths(spec.path)
    src_mt = _max_mtime(paths)

    if cache_path.exists() and cache_path.stat().st_mtime >= src_mt:
        LOG.info("cache hit: %s", cache_path.name)
        data = np.load(cache_path, allow_pickle=False)
        return LoadedSample(
            spec=spec,
            z_h=data["z_h"] if "z_h" in data.files else None,
            z_mb=data["z_mb"] if "z_mb" in data.files else None,
            z_r=data["z_r"],
            runnumber=data["runnumber"] if "runnumber" in data.files else None,
            trigger_passed=data["trigger_passed"] if "trigger_passed" in data.files else None,
        )

    LOG.info("cache miss: reading %s", spec.name)
    sample = load_sample(spec, branches, n_events=n_events, require_mb=require_mb,
                         trigger_bit=trigger_bit)

    save: Dict[str, np.ndarray] = {"z_r": sample.z_r}
    if sample.z_h is not None:
        save["z_h"] = sample.z_h
    if sample.z_mb is not None:
        save["z_mb"] = sample.z_mb
    if sample.runnumber is not None:
        save["runnumber"] = sample.runnumber
    if sample.trigger_passed is not None:
        save["trigger_passed"] = sample.trigger_passed
    np.savez(cache_path, **save)
    LOG.info("cached %s to %s (%d events)", spec.name, cache_path.name, len(sample.z_r))
    return sample


# ---------------------------------------------------------------------------
# Event-level histogram building (the forward model)
# ---------------------------------------------------------------------------

def _valid_mask(*arrs: Optional[np.ndarray]) -> Optional[np.ndarray]:
    m = None
    for a in arrs:
        if a is None:
            continue
        mm = np.isfinite(a) & (a > SENTINEL_F + 1)
        m = mm if m is None else (m & mm)
    return m


def _period_mask(run: Optional[np.ndarray], run_min: int, run_max: int) -> Optional[np.ndarray]:
    if run is None:
        return None
    return (run >= run_min) & (run < run_max)


def mc_reco_histogram_single(
    sample: LoadedSample,
    binning: Binning,
    w_func: Callable[[np.ndarray], np.ndarray],
    c_func: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> np.ndarray:
    """Build the reco-vertex histogram contribution from one single-MC
    sample, with per-event weight = xs_weight * w(z_h) * [c(z_h)].
    """
    if sample.z_h is None:
        raise ValueError(f"sample {sample.spec.name}: vertexz_truth required for MC histogram")
    mask = _valid_mask(sample.z_h, sample.z_r)
    if mask is None or not mask.any():
        return np.zeros(binning.n_bins)
    z_h = sample.z_h[mask]
    z_r = sample.z_r[mask]
    w_ev = np.full(len(z_h), sample.spec.event_weight) * w_func(z_h)
    if c_func is not None:
        w_ev *= c_func(z_h)
    h, _ = np.histogram(z_r, bins=binning.edges, weights=w_ev)
    return h


def mc_reco_histogram_double(
    sample: LoadedSample,
    binning: Binning,
    w_func: Callable[[np.ndarray], np.ndarray],
    c_func: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> np.ndarray:
    """Build the reco-vertex histogram contribution from one double-MC
    sample, with per-event weight = xs_weight * w(z_h) * w(z_mb)
    * [c(z_h) * c(z_mb)]. This is the "apply w once per vertex" prescription."""
    if sample.z_h is None or sample.z_mb is None:
        raise ValueError(f"sample {sample.spec.name}: double MC requires both vertexz_truth and vertexz_truth_mb")
    mask = _valid_mask(sample.z_h, sample.z_mb, sample.z_r)
    if mask is None or not mask.any():
        return np.zeros(binning.n_bins)
    z_h = sample.z_h[mask]
    z_mb = sample.z_mb[mask]
    z_r = sample.z_r[mask]
    w_ev = np.full(len(z_h), sample.spec.event_weight) * w_func(z_h) * w_func(z_mb)
    if c_func is not None:
        w_ev *= c_func(z_h) * c_func(z_mb)
    h, _ = np.histogram(z_r, bins=binning.edges, weights=w_ev)
    return h


def mixed_mc_histogram(
    single_samples: Sequence[LoadedSample],
    double_samples: Sequence[LoadedSample],
    binning: Binning,
    f_single: float,
    f_double: float,
    w_func: Callable[[np.ndarray], np.ndarray],
    c_func: Optional[Callable[[np.ndarray], np.ndarray]] = None,
    return_components: bool = False,
) -> np.ndarray:
    """Build the mixed-MC reco histogram at physical pile-up fractions.
    If return_components=True, returns (M_single_f_scaled, M_double_f_scaled)
    with the f_s / f_d factors already applied.
    """
    n = binning.n_bins
    M_s = np.zeros(n)
    for s in single_samples:
        M_s += mc_reco_histogram_single(s, binning, w_func, c_func)
    M_d = np.zeros(n)
    for s in double_samples:
        M_d += mc_reco_histogram_double(s, binning, w_func, c_func)
    M_s *= f_single
    M_d *= f_double
    if return_components:
        return M_s, M_d
    return M_s + M_d


def data_reco_histogram(
    data_sample: LoadedSample,
    binning: Binning,
    period: Optional[Tuple[int, int]] = None,
) -> np.ndarray:
    """Build the unweighted data reco-vertex histogram. If period =
    (run_min, run_max) is given, filter events by runnumber.
    If the sample has a ``trigger_passed`` mask attached, also require it
    (mirrors the analysis pipeline's trigger filter)."""
    mask = _valid_mask(data_sample.z_r)
    if period is not None:
        pmask = _period_mask(data_sample.runnumber, period[0], period[1])
        if pmask is None:
            raise ValueError(
                f"period filter requested but data sample has no runnumber"
            )
        mask = (mask & pmask) if mask is not None else pmask
    if data_sample.trigger_passed is not None:
        tmask = data_sample.trigger_passed.astype(bool)
        mask = (mask & tmask) if mask is not None else tmask
    if mask is None or not mask.any():
        return np.zeros(binning.n_bins)
    h, _ = np.histogram(data_sample.z_r[mask], bins=binning.edges)
    return h


# ---------------------------------------------------------------------------
# Parametric weight and σ_gen measurement
# ---------------------------------------------------------------------------

def parametric_weight(
    z: np.ndarray, sigma_gen: float, sigma_tgt: float,
) -> np.ndarray:
    """Gaussian ratio w(z) = exp(0.5 z² (1/σ_gen² − 1/σ_tgt²))."""
    z = np.asarray(z, dtype=np.float64)
    return np.exp(0.5 * z * z * (1.0 / (sigma_gen * sigma_gen) - 1.0 / (sigma_tgt * sigma_tgt)))


def fit_sigma_gen(
    nominal_samples: Iterable[LoadedSample],
) -> Tuple[float, float]:
    """Estimate σ_gen by computing the weighted std of the combined
    single-interaction truth vertex distribution. Returns (sigma, err).
    """
    zs: List[np.ndarray] = []
    ws: List[np.ndarray] = []
    for s in nominal_samples:
        if s.z_h is None:
            continue
        mask = _valid_mask(s.z_h)
        if mask is None:
            continue
        zs.append(s.z_h[mask])
        ws.append(np.full(int(mask.sum()), s.spec.event_weight))
    if not zs:
        raise ValueError("no nominal samples with vertexz_truth")
    z = np.concatenate(zs)
    w = np.concatenate(ws)
    w = w / w.sum()
    mean = float(np.sum(w * z))
    var = float(np.sum(w * (z - mean) ** 2))
    sigma = float(np.sqrt(var))
    n_eff = 1.0 / float(np.sum(w * w))
    err = sigma / np.sqrt(2.0 * max(n_eff, 1.0))
    return sigma, err


def chi2(
    D: np.ndarray,
    M: np.ndarray,
    mask: Optional[np.ndarray] = None,
) -> float:
    """Poisson-error χ² between data D and model M. If mask is given,
    restrict the sum to bins where mask is True."""
    err2 = np.where(D > 0, D, 1.0) + np.where(M > 0, M, 1.0)
    contrib = (D - M) ** 2 / np.where(err2 > 0, err2, 1.0)
    if mask is not None:
        contrib = contrib[mask]
    return float(contrib.sum())


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def c_interpolator(
    centers: np.ndarray, c_values: np.ndarray,
) -> Callable[[np.ndarray], np.ndarray]:
    """Return a callable c(z) that linearly interpolates c_values on
    bin centers. Outside the range, c returns the nearest edge value."""
    def _c(z: np.ndarray) -> np.ndarray:
        return np.interp(z, centers, c_values,
                         left=float(c_values[0]), right=float(c_values[-1]))
    return _c


def summary_stats(arr: np.ndarray, label: str = "") -> str:
    mask = np.isfinite(arr) & (arr > SENTINEL_F + 1)
    if not mask.any():
        return f"{label}: no finite entries"
    a = arr[mask]
    return (
        f"{label}: n={len(a)} min={a.min():+.3f} max={a.max():+.3f} "
        f"mean={a.mean():+.3f} std={a.std(ddof=1):.3f}"
    )
