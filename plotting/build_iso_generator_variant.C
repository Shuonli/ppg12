// build_iso_generator_variant.C
//
// Construct a synthetic systematic-variant cross-section file that captures
// ONLY the iso-efficiency difference between the nominal PYTHIA8 Detroit MC
// and the HERWIG7 cross-check (1.5 mrad). The full HERWIG re-run carries
// differences in iso, ID, MBD, and response-matrix; we propagate only the
// iso-eff difference into the systematic quadrature, since ID + unfolding
// are sub-leading cross-checks (Appendix HERWIG xcheck) and MBD-eff is
// covered by the data-driven MBD-coincidence systematic (PPG10).
//
// Sigma scaling: sigma ~ N_data / (eps_iso × ...), so freezing all factors
// except eps_iso gives sigma_iso_only = sigma_nom × eps_iso_P / eps_iso_H.
//
// Output: efficiencytool/results/Photon_final_bdt_iso_generator.root
//         containing a clone of h_unfold_sub_result rescaled bin-by-bin.

#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH1.h>
#include <TEfficiency.h>
#include <iostream>

void build_iso_generator_variant()
{
    const std::string base = "/gpfs/mnt/gpfs02/sphenix/user/shuhangli/ppg12/efficiencytool/results/";
    TFile *f_nom = TFile::Open((base + "Photon_final_bdt_nom.root").c_str(), "READ");
    TFile *f_pyt = TFile::Open((base + "MC_efficiency_bdt_nom_1p5mrad_newrewt.root").c_str(), "READ");
    TFile *f_her = TFile::Open((base + "MC_efficiency_bdt_herwig_1p5mrad_common.root").c_str(), "READ");
    if (!f_nom || f_nom->IsZombie() || !f_pyt || f_pyt->IsZombie() ||
        !f_her || f_her->IsZombie()) {
        std::cerr << "ERROR: cannot open inputs" << std::endl;
        return;
    }

    TH1 *h_xs_nom = (TH1 *)f_nom->Get("h_unfold_sub_result");
    TEfficiency *e_p = (TEfficiency *)f_pyt->Get("eff_iso_eta_0");
    TEfficiency *e_h = (TEfficiency *)f_her->Get("eff_iso_eta_0");
    if (!h_xs_nom || !e_p || !e_h) {
        std::cerr << "ERROR: input histograms missing" << std::endl;
        return;
    }

    TH1 *h_var = (TH1 *)h_xs_nom->Clone("h_unfold_sub_result");
    h_var->SetDirectory(0);

    int nbins = h_var->GetNbinsX();
    TH1 *h_eff_axis = (TH1 *)e_p->GetTotalHistogram()->Clone("h_eff_axis");
    for (int i = 1; i <= nbins; ++i) {
        double pt = h_var->GetBinCenter(i);
        int j = h_eff_axis->FindBin(pt);
        double iso_p = e_p->GetEfficiency(j);
        double iso_h = e_h->GetEfficiency(j);
        if (iso_p <= 0 || iso_h <= 0) continue;
        double scale = iso_p / iso_h;
        h_var->SetBinContent(i, h_xs_nom->GetBinContent(i) * scale);
        h_var->SetBinError(i, h_xs_nom->GetBinError(i) * scale);
        printf("  bin %2d  pT=%5.1f  scale = %.4f / %.4f = %.4f   sigma %.4e -> %.4e\n",
               i, pt, iso_p, iso_h, scale,
               h_xs_nom->GetBinContent(i), h_var->GetBinContent(i));
    }

    const std::string out = base + "Photon_final_bdt_iso_generator.root";
    TFile *fout = new TFile(out.c_str(), "RECREATE");
    h_var->Write();
    fout->Close();
    std::cout << "Wrote " << out << std::endl;
}
