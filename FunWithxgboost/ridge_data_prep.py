from pathlib import Path

import numpy as np
import pandas as pd


GW_EXCEL = Path("GoyalWelch_Data2024.xlsx")
FF_CSV = Path("F-F_Research_Data_Factors.csv")
OUT_CSV = Path("ridge_panel.csv")


def load_goyal_welch_monthly() -> pd.DataFrame:
    """
    Load the monthly Goyal-Welch predictors from the 2024 Excel file.

    Returns a DataFrame with at least:
      - 'yyyymm' (int)
      - the 14 predictors used in Kelly-Malamud-Zhou:
        dfy, infl, svar, de, lty, tms, tbl, dfr, dp, dy, ltr, ep, b/m, ntis
    """
    if not GW_EXCEL.exists():
        raise FileNotFoundError(f"Cannot find {GW_EXCEL}")

    # Monthly sheet contains the actual time series
    df_gw = pd.read_excel(GW_EXCEL, sheet_name="Monthly")

    # Map paper mnemonics -> column names in this file
    predictor_map = {
        "dfy": "dfy",
        "infl": "infl",
        "svar": "svar",
        "de": "d/e",
        "lty": "lty",
        "tms": "tms",
        "tbl": "tbl",
        "dfr": "dfr",
        "dp": "d/p",
        "dy": "d/y",
        "ltr": "ltr",
        "ep": "e/p",
        "b/m": "b/m",
        "ntis": "ntis",
    }

    missing = [v for v in predictor_map.values() if v not in df_gw.columns]
    if missing:
        raise KeyError(f"Missing expected Goyal-Welch columns: {missing}")

    cols = ["yyyymm"] + list(predictor_map.values())
    out = df_gw[cols].copy()
    # Rename to clean paper-style mnemonics
    out = out.rename(columns={v: k for k, v in predictor_map.items()})

    # Ensure yyyymm is integer
    out["yyyymm"] = out["yyyymm"].astype(int)
    return out


def load_fama_french_monthly() -> pd.DataFrame:
    """
    Load the Fama-French research factors (Mkt-RF, RF) from CSV.

    The file has a short text header; we skip it and parse the numeric panel.
    """
    if not FF_CSV.exists():
        raise FileNotFoundError(f"Cannot find {FF_CSV}")

    # Skip the 3 lines of descriptive header
    df_ff = pd.read_csv(FF_CSV, skiprows=3)

    # First column is the YYYYMM date index with no header name
    if "Unnamed: 0" not in df_ff.columns:
        raise KeyError("Expected first column 'Unnamed: 0' with YYYYMM dates.")

    df_ff = df_ff.rename(columns={"Unnamed: 0": "yyyymm"})

    # Drop non-monthly/footer rows by coercing to numeric
    df_ff["yyyymm"] = pd.to_numeric(df_ff["yyyymm"], errors="coerce")
    df_ff = df_ff.dropna(subset=["yyyymm"])
    df_ff["yyyymm"] = df_ff["yyyymm"].astype(int)

    # Keep only the standard four factors
    keep_cols = ["yyyymm", "Mkt-RF", "SMB", "HML", "RF"]
    missing = [c for c in keep_cols if c not in df_ff.columns]
    if missing:
        raise KeyError(f"Missing expected Fama-French columns: {missing}")
    df_ff = df_ff[keep_cols].copy()

    # Convert returns from strings to floats, if needed
    for c in ["Mkt-RF", "SMB", "HML", "RF"]:
        df_ff[c] = pd.to_numeric(df_ff[c], errors="coerce")

    # Drop rows with missing market excess return
    df_ff = df_ff.dropna(subset=["Mkt-RF"])

    return df_ff


def volatility_standardize_returns(df: pd.DataFrame, ret_col: str, window: int = 12) -> pd.DataFrame:
    """
    Volatility-standardize a return series by its trailing window RMS, as in the paper:
      sigma_t = sqrt( mean_{i=0..window-1} r_{t-i}^2 )
      r_std_t = r_t / sigma_t

    Uses uncentered second moment (RMS) rather than sample variance.
    """
    r = df[ret_col] / 100.0  # convert percent to decimal

    def rms(x: pd.Series) -> float:
        return float(np.sqrt(np.mean(np.square(x))))

    rms_series = r.rolling(window=window, min_periods=window).apply(rms, raw=False)
    df[f"{ret_col}_rms{window}"] = rms_series
    df[f"{ret_col}_std"] = r / rms_series
    return df


def expanding_standardize_predictors(df: pd.DataFrame, cols: list[str], min_periods: int = 36) -> pd.DataFrame:
    """
    Standardize predictors using expanding-window historical standard deviation
    (centered sample std), with an initial burn-in of `min_periods`.
    """
    for c in cols:
        series = df[c].astype(float)
        std_exp = series.expanding(min_periods=min_periods).std()
        df[f"{c}_std"] = series / std_exp
    return df


def build_panel() -> pd.DataFrame:
    """
    Construct the monthly panel used in Kelly-Malamud-Zhou:
      - target: market excess return (Mkt-RF) and RF from Fama-French
      - predictors: 14 Goyal-Welch variables + lagged market return
      - standardizations:
          * returns: trailing 12-month RMS
          * predictors: expanding-window std with 36-month burn-in
    """
    df_gw = load_goyal_welch_monthly()
    df_ff = load_fama_french_monthly()

    # Merge on yyyymm
    df = pd.merge(df_ff, df_gw, on="yyyymm", how="left", validate="one_to_one")

    # Restrict to the paper's typical sample range 1926-2020
    df = df[(df["yyyymm"] >= 192601) & (df["yyyymm"] <= 202012)].copy()
    df = df.sort_values("yyyymm").reset_index(drop=True)

    # Target: market excess return and risk-free
    df = volatility_standardize_returns(df, "Mkt-RF", window=12)

    # Lagged market excess return as a predictor (R_t to forecast R_{t+1})
    # Use the raw excess return (percent) here; it will be standardized as a predictor.
    df["mkt_excess_lag1"] = df["Mkt-RF"].shift(1)

    # Goyal-Welch predictors (paper's 14 variables)
    predictor_cols = [
        "dfy",
        "infl",
        "svar",
        "de",
        "lty",
        "tms",
        "tbl",
        "dfr",
        "dp",
        "dy",
        "ltr",
        "ep",
        "b/m",
        "ntis",
        "mkt_excess_lag1",  # 15th predictor: lagged market return
    ]

    # Standardize predictors with expanding-window std (36-month burn-in)
    df = expanding_standardize_predictors(df, predictor_cols, min_periods=36)

    # Drop rows before the standardization burn-in
    mask_valid = df[[f"{c}_std" for c in predictor_cols]].notna().all(axis=1)
    df = df[mask_valid].copy()

    return df


def main() -> None:
    df = build_panel()

    # Save a compact panel with:
    #   yyyymm, raw and standardized market excess return and RF,
    #   raw and standardized predictors (15 total).
    predictor_cols = [
        "dfy",
        "infl",
        "svar",
        "de",
        "lty",
        "tms",
        "tbl",
        "dfr",
        "dp",
        "dy",
        "ltr",
        "ep",
        "b/m",
        "ntis",
        "mkt_excess_lag1",
    ]

    cols_out: list[str] = [
        "yyyymm",
        "Mkt-RF",
        "Mkt-RF_rms12",
        "Mkt-RF_std",
        "RF",
    ]
    for c in predictor_cols:
        cols_out.append(c)
        cols_out.append(f"{c}_std")

    df_out = df[cols_out].copy()
    df_out.to_csv(OUT_CSV, index=False)

    print(f"Saved panel with {len(df_out)} rows and {len(df_out.columns)} columns to {OUT_CSV}")


if __name__ == "__main__":
    main()



