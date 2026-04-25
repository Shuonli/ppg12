// make_tower_masks.C
//
// Generate 3 binary tower-mask TH2I histograms (one per selection level:
// preselect, common, tight) via the DATA-ONLY phi-symmetry method:
//   - For each 2-eta-row band (48 bands total), fit a Gaussian to the 256
//     per-phi data counts (3-sigma iterative clip).
//   - Flag phi bins with z < -2 within each band as dead towers.
// No MC assumption needed; pp physics is phi-symmetric at fixed eta.
//
//   0 = keep,  1 = veto the cluster if its center tower lands here
//
// Three masks emitted:
//   mask_phisymm_preselect
//   mask_phisymm_common
//   mask_phisymm_tight
//
// Input : /sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root
// Output: /sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root

#include "../plotting/plotcommon.h"

static std::pair<double, double> fit_gauss_3sigma(const std::vector<double> &x)
{
    if (x.empty()) return {0.0, 1.0};
    double mean = 0, sigma = 0;
    for (double v : x) mean += v;
    mean /= x.size();
    for (double v : x) sigma += (v - mean) * (v - mean);
    sigma = std::sqrt(sigma / x.size());
    for (int it = 0; it < 3; ++it) {
        double m2 = 0, s2 = 0; int nk = 0;
        for (double v : x) if (std::fabs(v - mean) < 3.0 * sigma) { m2 += v; nk++; }
        if (nk == 0) break;
        mean = m2 / nk;
        for (double v : x) if (std::fabs(v - mean) < 3.0 * sigma) s2 += (v - mean) * (v - mean);
        sigma = std::sqrt(s2 / nk);
    }
    return {mean, sigma};
}

static TH2I *build_phisymm_mask(TH2F *h_data, const std::string &name,
                                const std::string &title, int eta_group = 2)
{
    const int nx = h_data->GetNbinsX();
    const int ny = h_data->GetNbinsY();
    const int n_bands = nx / eta_group;

    TH2I *h_mask = new TH2I(name.c_str(), title.c_str(),
                            nx, 0, nx, ny, 0, ny);
    h_mask->SetXTitle("cluster i#it{#eta} (tower)");
    h_mask->SetYTitle("cluster i#it{#phi} (tower)");
    h_mask->SetZTitle("mask (1 = veto)");

    int n_masked = 0;
    for (int b = 0; b < n_bands; ++b) {
        // Sum ETA_GROUP consecutive eta rows into one "band"
        std::vector<double> band(ny, 0.0);
        for (int de = 0; de < eta_group; ++de) {
            int ix = b * eta_group + de + 1;  // ROOT bins are 1-indexed
            for (int iy = 1; iy <= ny; ++iy)
                band[iy - 1] += h_data->GetBinContent(ix, iy);
        }
        // Skip bands with too few events
        double sum = 0;
        for (double v : band) sum += v;
        if (sum < 10) continue;
        auto stats = fit_gauss_3sigma(band);
        double mu = stats.first, sig = stats.second;
        if (sig <= 0) continue;
        // Flag phi bins with z < -2
        for (int iy = 0; iy < ny; ++iy) {
            double z = (band[iy] - mu) / sig;
            if (z < -2) {
                // Fill the original eta rows of this band at this iphi
                for (int de = 0; de < eta_group; ++de) {
                    int ix = b * eta_group + de + 1;
                    h_mask->SetBinContent(ix, iy + 1, 1);
                    n_masked++;
                }
            }
        }
    }
    std::cout << "[" << name << "] masked towers = " << n_masked << std::endl;
    return h_mask;
}

void make_tower_masks()
{
    const char *infile =
        "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root";
    const char *outfile =
        "/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root";

    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "cannot open " << infile << std::endl;
        return;
    }

    SetsPhenixStyle();

    TFile *fo = TFile::Open(outfile, "RECREATE");

    std::map<std::string, TH2I *> per_level_masks;
    for (const std::string &lvl : {"preselect", "common", "tight"}) {
        TH2F *h_da = (TH2F *) f->Get(Form("h_etaphi_tower_%s_data", lvl.c_str()));
        if (!h_da) {
            std::cerr << "missing data histo for level " << lvl << std::endl;
            continue;
        }
        std::string name = "mask_phisymm_" + lvl;
        std::string title = "Tower mask (phi-symmetry, z<-2, 2-eta-row bands) " + lvl;
        TH2I *h = build_phisymm_mask(h_da, name, title);
        fo->cd();
        h->Write();
        per_level_masks[lvl] = h;
    }

    // OR mask: union of common + tight phi-symm masks. Broader coverage
    // than either single-level; useful as a conservative systematic test.
    TH2I *hc = per_level_masks["common"];
    TH2I *ht = per_level_masks["tight"];
    if (hc && ht) {
        TH2I *h_or = new TH2I("mask_phisymm_or",
            "Tower mask phi-symmetry OR(common, tight)",
            hc->GetNbinsX(), 0, hc->GetNbinsX(),
            hc->GetNbinsY(), 0, hc->GetNbinsY());
        h_or->SetXTitle("cluster i#it{#eta} (tower)");
        h_or->SetYTitle("cluster i#it{#phi} (tower)");
        h_or->SetZTitle("mask (1 = veto)");
        int n_or = 0;
        for (int ix = 1; ix <= hc->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= hc->GetNbinsY(); ++iy) {
                int v = (hc->GetBinContent(ix, iy) > 0 || ht->GetBinContent(ix, iy) > 0) ? 1 : 0;
                h_or->SetBinContent(ix, iy, v);
                if (v) n_or++;
            }
        std::cout << "[mask_phisymm_or] masked towers = " << n_or << std::endl;
        fo->cd();
        h_or->Write();
    }

    fo->Close();
    f->Close();
    std::cout << "wrote " << outfile << std::endl;
}
