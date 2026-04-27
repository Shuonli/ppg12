#include "plotcommon.h"

#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldBinByBin.h>

RooUnfoldResponse *make_response(TH2D *h2d, std::string name)
{
    TH2D *h2d_clone = (TH2D *)h2d->Clone("h2d_clone");
    h2d_clone->Reset();
    //set the bin based on h2d
    
    for (int i = 1; i <= h2d->GetNbinsX(); i++)
    {
        for (int j = 1; j <= h2d->GetNbinsY(); j++)
        {
            double bincontent = h2d->GetBinContent(i, j);
            double binerror = h2d->GetBinError(i, j);
            
            h2d_clone->SetBinContent(i, j, bincontent);
            h2d_clone->SetBinError(i, j, binerror);
        }
    }
    

    TH1D *h_truth = (TH1D *)h2d_clone->ProjectionY();
    TH1D *h_reco = (TH1D *)h2d_clone->ProjectionX();
    for (int i = 1; i <= h_reco->GetNbinsX(); i++)
    {
       double bincontent = h_reco->GetBinContent(i);
       double res = bincontent  - std::floor(bincontent);
       std::cout << "bincontent: " << bincontent << " res: " << res << std::endl;
    }
    RooUnfoldResponse *response = new RooUnfoldResponse(h_reco, h_truth, h2d_clone, name.c_str(), "", false);
    return response;
}

void plot_unfold_iter()
{
    init_plot();

    string savePath = "figures/";

    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root");
    static const int ntotal_iterations = 10;
    TFile *fresponsein = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_response_bdt_nom.root");

    TFile *fresponsehist = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root");

    // RooUnfoldResponse *response_full = (RooUnfoldResponse *)fresponsein->Get("response_matrix_full_0");
    // RooUnfoldResponse *response_half = (RooUnfoldResponse *)fresponsein->Get("response_matrix_half_0");
    TH1D *h_pT_truth_response = (TH1D *)fresponsein->Get("h_pT_truth_response_0");
    TH1D *h_pT_reco_response = (TH1D *)fresponsein->Get("h_pT_reco_response_0");

    //reset the overflow and underflow bin
    h_pT_truth_response->SetBinContent(0, 0);
    h_pT_reco_response->SetBinContent(0, 0);
    h_pT_truth_response->SetBinContent(h_pT_truth_response->GetNbinsX() + 1, 0);
    h_pT_reco_response->SetBinContent(h_pT_reco_response->GetNbinsX() + 1, 0);

    TH1D *h_pT_truth_secondhalf_response = (TH1D *)fresponsein->Get("h_pT_truth_secondhalf_response_0");
    TH1D *h_pT_reco_secondhalf_response = (TH1D *)fresponsein->Get("h_pT_reco_secondhalf_response_0");

    TH2D *h_response_full = (TH2D *)fresponsehist->Get("h_response_full_0");
    TH2D *h_response_half = (TH2D *)fresponsehist->Get("h_response_half_0");

    TH1D *h_response_full_pfx = h_response_full->ProjectionX("h_response_full_pfx");
    TH1D *h_response_half_pfx = h_response_half->ProjectionX("h_response_half_pfx");

    TH1D *h_response_full_pfy = h_response_full->ProjectionY("h_response_full_pfy");
    TH1D *h_response_half_pfy = h_response_half->ProjectionY("h_response_half_pfy");

    TH1D *h_response_full_pfx_ratio = (TH1D *)h_response_full_pfx->Clone("h_response_full_pfx_ratio");
    h_response_full_pfx_ratio->Divide(h_pT_reco_response);
    TH1D *h_response_half_pfx_ratio = (TH1D *)h_response_half_pfx->Clone("h_response_half_pfx_ratio");
    h_response_half_pfx_ratio->Divide(h_pT_reco_response);

    TH1D *h_response_full_pfy_ratio = (TH1D *)h_response_full_pfy->Clone("h_response_full_pfy_ratio");
    h_response_full_pfy_ratio->Divide(h_pT_truth_response);
    TH1D *h_response_half_pfy_ratio = (TH1D *)h_response_half_pfy->Clone("h_response_half_pfy_ratio");
    h_response_half_pfy_ratio->Divide(h_pT_truth_response);

    std::cout << "h_pT_truth_response: " << h_pT_truth_response->GetEntries() << std::endl;
    std::cout << "h_pT_reco_response: " << h_pT_reco_response->GetEntries() << std::endl;
    RooUnfoldResponse *response_full = make_response(h_response_full, "response_full");
    RooUnfoldResponse *response_half = make_response(h_response_half, "response_half");

    string leg1 = "total relative stat. uncertainty";
    string leg2 = "total relative deviation";

    std::string h_result_name_base = "h_unfold_sub_leak_";

    TH1F *h_stat = new TH1F("h_stat", "Statistical Uncertainty", ntotal_iterations, 0.5, ntotal_iterations + 0.5);

    TH1F *h_iter_delta = new TH1F("h_iter_delta", "Iteration Delta", ntotal_iterations, 0.5, ntotal_iterations + 0.5);
    
        for (int i = 1; i < ntotal_iterations; i++)
        {
            std::string h_prev_name = h_result_name_base + std::to_string(i);
            std::string h_this_name = h_result_name_base + std::to_string(i + 1);

            TH1F *h_prev = (TH1F *)fin->Get(h_prev_name.c_str());
            TH1F *h_this = (TH1F *)fin->Get(h_this_name.c_str());

            TH1F *h_dev_rel = (TH1F *)calcDelta(h_this, h_prev, "h_dev_rel_" + std::to_string(i)).second;

            // loop over bin find the quad sum
            double rel_quad_sum = 0;
            double rel_stat_sum = 0;

            //for (int ibin = 1; ibin <= h_dev_rel->GetNbinsX(); ibin++)
            for (int ibin = 2; ibin <= (h_dev_rel->GetNbinsX() - 2); ibin++)
            {
                double rel = h_dev_rel->GetBinContent(ibin);
                double stat = h_this->GetBinError(ibin) / h_this->GetBinContent(ibin);
                //std::cout << "ibin: " << ibin << " rel: " << rel << " stat: " << stat << std::endl;
                if (rel == rel)
                {
                    rel_quad_sum += rel * rel;
                }
                if (stat == stat)
                {
                    rel_stat_sum += stat * stat;
                }
            }
            rel_quad_sum = sqrt(rel_quad_sum);
            rel_stat_sum = sqrt(rel_stat_sum);

            std::cout << "rel_quad_sum: " << rel_quad_sum << " rel_stat_sum: " << rel_stat_sum << std::endl;

            h_stat->SetBinContent(i + 1, rel_stat_sum);
            h_stat->SetBinError(i + 1, 0);
            h_iter_delta->SetBinContent(i + 1, rel_quad_sum);
            h_iter_delta->SetBinError(i + 1, 0);
        }

        TCanvas *c1 = new TCanvas("can", "", 600, 600);

        frame_iteration->SetYTitle("#sqrt{#delta}");
        frame_iteration->GetYaxis()->SetRangeUser(0, 0.6);
        frame_iteration->GetXaxis()->SetRangeUser(1, 10);
        frame_iteration->Draw("axis");

        h_stat->SetMarkerStyle(20);
        h_stat->SetMarkerColor(kBlack);
        h_stat->SetLineColor(kBlack);
        h_stat->Draw("same p");

        h_iter_delta->SetMarkerStyle(20);
        h_iter_delta->SetMarkerColor(kBlue);
        h_iter_delta->SetLineColor(kBlue);
        h_iter_delta->Draw("same p");

        myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
        myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
        myText(0.5, 0.80, 1, strleg3.c_str(), 0.04);

        myMarkerLineText(0.25, 0.75, 1, kBlack, 20, kBlack, 1, leg1.c_str(), 0.04, true);
        myMarkerLineText(0.25, 0.70, 1, kBlue, 20, kBlue, 1, leg2.c_str(), 0.04, true);

        c1->SaveAs(Form("%s/unfold_iter.pdf", savePath.c_str()));
    

    // closure test
    int niterations = 3;
    // full closure test
    TCanvas *c2 = new TCanvas("can2", "", 800, 889);
    c2->Divide(1, 2);

    TPad *pad_1 = (TPad *)c2->cd(1);
    pad_1->SetPad(0, 0.4, 1, 1);
    pad_1->SetTopMargin(0.05);
    pad_1->SetLeftMargin(0.13);
    pad_1->SetBottomMargin(0.03);
    pad_1->SetRightMargin(0.08);
    pad_1->SetLogy();

    frame_et_rec->SetYTitle("d#sigma/dE_{T} [pb/GeV]");
    frame_et_rec->GetYaxis()->SetRangeUser(1e3, 1e8);
    // Start at 10 GeV — drop the [8,10] truth bin from view: it has no reco
    // counterpart so the unfold there is extrapolation from the prior, and
    // the visible truth/iter offset at [8,10] just clutters the closure
    // verification at the physical truth range [10, 36].
    frame_et_rec->GetXaxis()->SetRangeUser(10, 36);
    frame_et_rec->SetYTitle("counts");

    frame_et_rec->GetXaxis()->SetTitleOffset(0.98);
    frame_et_rec->GetYaxis()->SetTitleOffset(1.15);
    frame_et_rec->GetXaxis()->SetLabelSize(0.045);
    frame_et_rec->GetYaxis()->SetLabelSize(0.045);
    frame_et_rec->GetXaxis()->SetLabelOffset(2);
    // frame_et_rec->GetXaxis()->CenterTitle();
    // frame_et_rec->GetYaxis()->CenterTitle();
    frame_et_rec->GetXaxis()->SetNdivisions(505);

    frame_et_rec->Draw("axis");

    h_pT_truth_response->SetMarkerStyle(20);
    h_pT_truth_response->SetMarkerColor(kBlack);
    h_pT_truth_response->SetLineColor(kBlack);
    h_pT_truth_response->Draw("same");

    h_pT_reco_response->SetLineColor(kRed);
    h_pT_reco_response->SetMarkerColor(kRed);
    h_pT_reco_response->SetMarkerStyle(20);
    h_pT_reco_response->Draw("same");

    // The truth and reco distributions differ at the boundary truth bins
    // [8,10] (where reco has no bin) and [10,12] (where reco accumulates the
    // events whose truth was [8,10]). In a perfect closure, iter k unfolded ==
    // truth for every k -- the iters thus all stack on top of truth in the
    // top panel and on top of each other. This is what the bottom panel
    // shows as ratio = 1. The red reco curve is NOT what the iters chase;
    // it is shown for context only (label "reco (input)" in legend).
    std::cout << "h_pT_reco_response: ";
    for (int ibin = 1; ibin <= h_pT_reco_response->GetNbinsX(); ibin++)
    {
        std::cout << h_pT_reco_response->GetBinContent(ibin) << " ";
    }
    std::cout << std::endl;

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "Full closure test", 0.05);

    myMarkerLineText(0.6, 0.75, 1, kBlack, 20, kBlack, 1, "truth", 0.05, true);
    myMarkerLineText(0.6, 0.70, 1, kRed, 20, kRed, 1, "reco (input)", 0.05, true);

    TPad *pad_2 = (TPad *)c2->cd(2);
    pad_2->SetPad(0, 0, 1, 0.4);
    pad_2->SetTopMargin(0.02);
    pad_2->SetLeftMargin(0.13);
    pad_2->SetBottomMargin(0.25);
    pad_2->SetRightMargin(0.08);

    frame_et_truth->SetYTitle("unfolded / truth");
    frame_et_truth->GetYaxis()->SetNdivisions(506);
    frame_et_truth->GetYaxis()->SetRangeUser(0.9, 1.1);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 36);
    frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
    frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
    frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
    frame_et_truth->GetXaxis()->SetNdivisions(505);
    frame_et_truth->Draw("axis");

    lineone->Draw("L");

    std::vector<TH1D *> h_unfolded;
    std::vector<TH1D *> h_unfolded_ratio;
    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kOrange + 10, kViolet + 3};
    // debug
    // h_pT_reco_response = (TH1D* )response_full->Hresponse()->ProjectionX();

    // Iteration markers: use distinct OPEN styles + size gradient (largest -> smallest)
    // so each iteration is visible even when central values are identical.
    // Without this, iter 3 (filled circle, drawn last) covers iter 1 and iter 2
    // and the user sees only one marker per bin -- no way to verify all
    // iterations agree visually.
    std::vector<int> iter_styles = {24, 25, 26, 32, 27};   // open: circle, square, triangle-up, triangle-dn, diamond
    std::vector<float> iter_sizes = {1.8f, 1.4f, 1.0f, 0.8f, 0.6f};
    for (int i = 0; i < niterations; i++)
    {
        c2->cd();
        RooUnfoldBayes *full_unfold = new RooUnfoldBayes(response_full, h_pT_reco_response, i + 1);
        // Use Bayesian iterative unfolder (matches the cross-section pipeline
        // and properly populates the [32,36) truth bin in the non-square
        // response). Bin-by-bin returns 0 for boundary truth bins under this
        // axis configuration, which made the iter-1/2/3 markers identical and
        // hid the [32,36) ratio off-frame.
        TH1D *h_unfolded_iter = (TH1D *)full_unfold->Hunfold();
        std::cout << "h_unfolded_iter: ";
        for (int ibin = 1; ibin <= h_unfolded_iter->GetNbinsX(); ibin++)
        {
            std::cout << h_unfolded_iter->GetBinContent(ibin) << " ";
        }
        std::cout << std::endl;
        pad_1->cd();
        h_unfolded_iter->SetMarkerStyle(iter_styles[i]);
        h_unfolded_iter->SetMarkerSize(iter_sizes[i]);
        h_unfolded_iter->SetMarkerColor(colors[i]);
        h_unfolded_iter->SetLineColor(colors[i]);
        h_unfolded_iter->Draw("same");
        std::cout << h_unfolded_iter->GetBinError(2) / h_unfolded_iter->GetBinContent(2) << std::endl;
        // h_unfolded.push_back(h_unfolded_iter);
        TH1D *h_unfolded_ratio_iter = (TH1D *)h_unfolded_iter->Clone(Form("h_unfolded_ratio_iter_%d", i));
        h_unfolded_ratio_iter->Divide(h_pT_truth_response);
        // h_unfolded_ratio.push_back(h_unfolded_ratio_iter);
        myMarkerLineText(0.6, 0.65 - 0.05 * i, 1, colors[i], iter_styles[i], colors[i], 1, Form("iteration %d", i + 1), 0.05, true);
        pad_2->cd();
        h_unfolded_ratio_iter->SetMarkerStyle(iter_styles[i]);
        h_unfolded_ratio_iter->SetMarkerSize(iter_sizes[i]);
        h_unfolded_ratio_iter->SetMarkerColor(colors[i]);
        h_unfolded_ratio_iter->SetLineColor(colors[i]);
        h_unfolded_ratio_iter->Draw("same");
    }
    c2->SaveAs(Form("%s/closure_full.pdf", savePath.c_str()));
    TCanvas *c3 = new TCanvas("can3", "", 800, 889);
    c3->Divide(1, 2);

    TPad *pad_3_top = (TPad *)c3->cd(1);
    pad_3_top->SetPad(0, 0.4, 1, 1);
    pad_3_top->SetTopMargin(0.05);
    pad_3_top->SetLeftMargin(0.13);
    pad_3_top->SetBottomMargin(0.03);
    pad_3_top->SetRightMargin(0.08);
    pad_3_top->SetLogy();

    // You can reuse the same frames or define a new pair for clarity.
    // For simplicity, let’s reuse frame_et_rec and frame_et_truth, but you may want
    // to clone or create new TFrames if you want different axis ranges, etc.

    // TOP PAD
    frame_et_rec->Draw("axis");
    frame_et_rec->SetYTitle("counts");
    frame_et_rec->GetYaxis()->SetRangeUser(1e3, 5e7); // covers data over 10-36 GeV
    frame_et_rec->GetXaxis()->SetRangeUser(10, 36);    // physical truth range; [8,10] dropped (no reco counterpart)

    // Draw the “truth” histogram for the half-sample closure
    // (assuming you have h_pT_truth_secondhalf_response as the "truth" distribution)
    h_pT_truth_secondhalf_response->SetMarkerStyle(20);
    h_pT_truth_secondhalf_response->SetMarkerColor(kBlack);
    h_pT_truth_secondhalf_response->SetLineColor(kBlack);
    h_pT_truth_secondhalf_response->Draw("same");

    myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
    myText(0.5, 0.80, 1, "Half closure test", 0.05);
    myMarkerLineText(0.6, 0.75, 1, kBlack, 20, kBlack, 1, "truth", 0.05, true);

    // BOTTOM PAD
    TPad *pad_3_bottom = (TPad *)c3->cd(2);
    pad_3_bottom->SetPad(0, 0, 1, 0.4);
    pad_3_bottom->SetTopMargin(0.02);
    pad_3_bottom->SetLeftMargin(0.13);
    pad_3_bottom->SetBottomMargin(0.25);
    pad_3_bottom->SetRightMargin(0.08);

    // Re-draw ratio axes
    frame_et_truth->SetYTitle("unfolded / truth");
    frame_et_truth->GetYaxis()->SetRangeUser(0.85, 1.15);
    frame_et_truth->GetXaxis()->SetRangeUser(10, 36);
    frame_et_truth->Draw("axis");

    // Horizontal line at y=1
    lineone->Draw("L");

    // Now unfold for the half-sample. Show iterations 1, 2, 4 (skip 3) — three
    // iterations are enough to bracket the unfolding behaviour and the plot
    // stays visually clean. Iter 1 = prior-dominated (loose), iter 2 = nominal
    // analysis count, iter 4 = converged (tight).
    std::vector<int> half_iters_to_show = {1, 2, 4};
    for (size_t k = 0; k < half_iters_to_show.size(); k++)
    {
        int iter_idx = half_iters_to_show[k];                // displayed label (1-based)
        int style_idx = static_cast<int>(k);                 // index into iter_styles/sizes/colors
        // Use response_half, h_pT_reco_secondhalf_response, etc.
        RooUnfoldBayes *half_unfold = new RooUnfoldBayes(response_half, h_pT_reco_secondhalf_response, iter_idx);
        TH1D *h_unfolded_half_iter = (TH1D *)half_unfold->Hunfold();

        // top pad — distinct open marker styles + size gradient (largest -> smallest)
        // so all iterations are visible even when central values agree.
        pad_3_top->cd();
        h_unfolded_half_iter->SetMarkerStyle(iter_styles[style_idx]);
        h_unfolded_half_iter->SetMarkerSize(iter_sizes[style_idx]);
        h_unfolded_half_iter->SetMarkerColor(colors[style_idx]);
        h_unfolded_half_iter->SetLineColor(colors[style_idx]);
        h_unfolded_half_iter->Draw("same");

        myMarkerLineText(0.6, 0.65 - 0.05 * style_idx, 1, colors[style_idx],
                         iter_styles[style_idx], colors[style_idx], 1,
                         Form("iteration %d", iter_idx), 0.05, true);

        // ratio
        TH1D *h_unfolded_half_ratio_iter = (TH1D *)h_unfolded_half_iter->Clone(
            Form("h_unfolded_half_ratio_iter_%d", iter_idx));
        h_unfolded_half_ratio_iter->Divide(h_pT_truth_secondhalf_response);

        pad_3_bottom->cd();
        h_unfolded_half_ratio_iter->SetMarkerStyle(iter_styles[style_idx]);
        h_unfolded_half_ratio_iter->SetMarkerSize(iter_sizes[style_idx]);
        h_unfolded_half_ratio_iter->SetMarkerColor(colors[style_idx]);
        h_unfolded_half_ratio_iter->SetLineColor(colors[style_idx]);
        h_unfolded_half_ratio_iter->Draw("same");
    }

    // Optionally save the half closure plot
    c3->SaveAs(Form("%s/closure_half.pdf", savePath.c_str()));
}
