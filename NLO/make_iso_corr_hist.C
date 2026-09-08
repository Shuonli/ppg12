// Build the isolation-correction histograms used to bring the inclusive PHENIX
// points onto the PPG12 isolated fiducial, from the CT18NLO JETPHOX spectra:
//   h_iso_over_incl = sigma_iso(R=0.3, ET<4) / sigma_incl(no iso)  (multiply PHENIX by this)
//   h_incl_over_iso = its inverse
// Binned in the PPG12 truth pT bins. TH1D::Divide propagates the (independent-MC) errors.
void make_iso_corr_hist()
{
    const char *base = "/sphenix/user/shuhangli/ppg12/NLO/rootFiles";
    TFile *fi = TFile::Open(Form("%s/jetPHOX_ct18iso200_10_chunked.root", base)); // inclusive
    TFile *fs = TFile::Open(Form("%s/jetPHOX_ct18_10_chunked.root", base));       // isolated (paper)
    TH1D *hi = (TH1D *) fi->Get("h_truth_pT");
    TH1D *hs = (TH1D *) fs->Get("h_truth_pT");

    TH1D *h_iso_over_incl = (TH1D *) hs->Clone("h_iso_over_incl");
    h_iso_over_incl->SetTitle("sigma_iso / sigma_incl (CT18NLO);p_{T}^{#gamma} [GeV];#sigma_{iso}/#sigma_{incl}");
    h_iso_over_incl->Divide(hi);

    TH1D *h_incl_over_iso = (TH1D *) hi->Clone("h_incl_over_iso");
    h_incl_over_iso->SetTitle("sigma_incl / sigma_iso (CT18NLO);p_{T}^{#gamma} [GeV];#sigma_{incl}/#sigma_{iso}");
    h_incl_over_iso->Divide(hs);

    TFile *fo = new TFile(Form("%s/iso_correction_ct18nlo.root", base), "RECREATE");
    h_iso_over_incl->Write();
    h_incl_over_iso->Write();
    fo->Close();

    printf("wrote %s/iso_correction_ct18nlo.root  (h_iso_over_incl, h_incl_over_iso)\n", base);
    for (int b = 1; b <= h_iso_over_incl->GetNbinsX(); ++b)
        printf("  [%5.1f,%5.1f]  iso/incl = %.4f +- %.4f\n",
               h_iso_over_incl->GetBinLowEdge(b), h_iso_over_incl->GetBinLowEdge(b + 1),
               h_iso_over_incl->GetBinContent(b), h_iso_over_incl->GetBinError(b));
}
