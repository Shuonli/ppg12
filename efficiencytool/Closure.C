#include <iostream>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <yaml-cpp/yaml.h>
// unfolding
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

#include "draw.C"

/*
 * Closure Test for Unfolding Validation
 *
 * FULL CLOSURE TEST:
 * - Use MC as pseudo-data
 * - Unfold using response matrix from the same MC
 * - Compare unfolded result to MC truth
 * - Should give perfect closure if unfolding is unbiased
 *
 * HALF CLOSURE TEST:
 * - Split MC into two halves
 * - Use first half to build response matrix
 * - Use second half as pseudo-data
 * - Unfold second half with first half's response matrix
 * - Compare to second half truth
 * - Tests for statistical independence and bias
 */

void Closure(const std::string &configname = "config_bdt_none.yaml", int niterations_max = 10, bool save_all_iterations = false)
{
    SetAtlasStyle();
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    std::cout << "\n========================================" << std::endl;
    std::cout << "UNFOLDING CLOSURE TEST" << std::endl;
    std::cout << "========================================\n" << std::endl;

    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    // Get response file from config
    std::string var_type = configYaml["output"]["var_type"].as<std::string>();
    std::string infilename = configYaml["output"]["response_outfile"].as<std::string>() + "_" + var_type + ".root";

    std::cout << "Configuration: " << configname << std::endl;
    std::cout << "Response file: " << infilename << std::endl;
    std::cout << "Max iterations: " << niterations_max << std::endl;

    std::vector<float> eta_bins = configYaml["analysis"]["eta_bins"].as<std::vector<float>>();

    TFile *fresin = new TFile(infilename.c_str(), "READ");
    if (!fresin || fresin->IsZombie()) {
        std::cerr << "ERROR: Cannot open response file: " << infilename << std::endl;
        return;
    }

    // Load response matrices and histograms
    std::vector<RooUnfoldResponse *> responses_full;
    std::vector<RooUnfoldResponse *> responses_half;
    std::vector<TH1D *> h_pT_truth_response;
    std::vector<TH1D *> h_pT_reco_response;
    std::vector<TH1D *> h_pT_truth_half_response;
    std::vector<TH1D *> h_pT_reco_half_response;
    std::vector<TH1D *> h_pT_truth_secondhalf_response;
    std::vector<TH1D *> h_pT_reco_secondhalf_response;

    std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3, kViolet + 3, kOrange + 10};
    std::vector<int> markerstyle = {20, 20, 20, 20, 20};

    // Create output directory if it doesn't exist
    gSystem->Exec("mkdir -p figure");

    // Output ROOT file for detailed results
    std::string output_rootfile = Form("ClosureTest_%s.root", configname.substr(0, configname.find(".yaml")).c_str());
    TFile *fout = new TFile(output_rootfile.c_str(), "RECREATE");

    for (int ieta = 0; ieta < (int)eta_bins.size() - 1; ieta++)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Processing eta bin " << ieta << ": [" << eta_bins[ieta] << ", " << eta_bins[ieta + 1] << "]" << std::endl;
        std::cout << "========================================" << std::endl;

        // Load response matrices
        responses_full.push_back((RooUnfoldResponse *)fresin->Get(Form("response_matrix_full_%d", ieta)));
        responses_half.push_back((RooUnfoldResponse *)fresin->Get(Form("response_matrix_half_%d", ieta)));

        if (!responses_full[ieta] || !responses_half[ieta]) {
            std::cerr << "ERROR: Cannot load response matrices for eta bin " << ieta << std::endl;
            continue;
        }

        // Load histograms
        h_pT_truth_response.push_back((TH1D *)fresin->Get(Form("h_pT_truth_response_%d", ieta)));
        h_pT_reco_response.push_back((TH1D *)fresin->Get(Form("h_pT_reco_response_%d", ieta)));
        h_pT_truth_half_response.push_back((TH1D *)fresin->Get(Form("h_pT_truth_half_response_%d", ieta)));
        h_pT_reco_half_response.push_back((TH1D *)fresin->Get(Form("h_pT_reco_half_response_%d", ieta)));
        h_pT_truth_secondhalf_response.push_back((TH1D *)fresin->Get(Form("h_pT_truth_secondhalf_response_%d", ieta)));
        h_pT_reco_secondhalf_response.push_back((TH1D *)fresin->Get(Form("h_pT_reco_secondhalf_response_%d", ieta)));

        if (!h_pT_truth_response[ieta] || !h_pT_reco_response[ieta]) {
            std::cerr << "ERROR: Cannot load histograms for eta bin " << ieta << std::endl;
            continue;
        }

        // ============================================
        // FULL CLOSURE TEST
        // ============================================
        std::cout << "\n--- FULL CLOSURE TEST ---" << std::endl;

        std::vector<TH1D *> h_unfolded_full_list;
        std::vector<double> chi2_full_list;
        std::vector<int> ndf_full_list;

        for (int iit = 0; iit < niterations_max; iit++)
        {
            RooUnfoldBayes *full_unfold = new RooUnfoldBayes(responses_full[ieta], h_pT_reco_response[ieta], iit + 1);
            TH1D *hRecoPT = (TH1D *)full_unfold->Hunfold();
            hRecoPT->SetName(Form("h_unfolded_full_eta%d_iter%d", ieta, iit + 1));

            // Calculate chi2
            double chi2 = 0;
            int ndf = 0;
            for (int ibin = 1; ibin <= hRecoPT->GetNbinsX(); ibin++) {
                double unf_val = hRecoPT->GetBinContent(ibin);
                double unf_err = hRecoPT->GetBinError(ibin);
                double tru_val = h_pT_truth_response[ieta]->GetBinContent(ibin);

                if (unf_err > 0 && tru_val > 0) {
                    chi2 += TMath::Power((unf_val - tru_val) / unf_err, 2);
                    ndf++;
                }
            }

            chi2_full_list.push_back(chi2);
            ndf_full_list.push_back(ndf);

            std::cout << "  Iteration " << iit + 1 << ": chi2/ndf = " << chi2 << "/" << ndf
                      << " = " << (ndf > 0 ? chi2 / ndf : 0) << std::endl;

            h_unfolded_full_list.push_back(hRecoPT);

            if (save_all_iterations) {
                hRecoPT->Write();
            }

            // Plot for key iterations (1, 2, 4)
            if (iit == 0 || iit == 1 || iit == 3)
            {
                std::vector<TH1F *> h_input;
                h_input.push_back((TH1F *)h_pT_truth_response[ieta]->Clone());
                h_input.push_back((TH1F *)h_pT_reco_response[ieta]->Clone());
                h_input.push_back((TH1F *)hRecoPT->Clone());

                std::vector<std::string> text;
                std::vector<std::string> legend;

                text.push_back("Full closure test");
                text.push_back(Form("Iteration: %d", iit + 1));
                text.push_back(Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf, ndf > 0 ? chi2 / ndf : 0));
                text.push_back(Form("%.1f < #eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]));
                text.push_back("|vtxz|<30cm, tight iso truth #gamma");

                legend.push_back("Truth photon spectrum");
                legend.push_back("Reco photon spectrum");
                legend.push_back("Unfolded photon spectrum");

                draw_1D_multiple_plot_ratio(h_input, colors, markerstyle,
                                            false, 10, true,
                                            false, 0, 40, false,
                                            false, 0, 0.5, true,
                                            true, 0.8, 1.2,
                                            true, "p_{T}^{#gamma} [GeV]", "dN_{#gamma}/dp_{T}^{#gamma}", "Unfolded/Truth",
                                            false, "sPHENIX Simulation",
                                            true, text, 0.25, 0.41, 0.04,
                                            true, legend, 0.58, 0.76, 0.04,
                                            Form("figure/FullClosure_eta%d_iter%d.pdf", ieta, iit + 1));

                for (auto h : h_input) delete h;
            }

            delete full_unfold;
        }

        // Plot chi2 vs iteration
        TCanvas *c_chi2_full = new TCanvas(Form("c_chi2_full_eta%d", ieta), "Chi2 vs Iteration", 800, 600);
        TGraphErrors *g_chi2_full = new TGraphErrors();
        for (int iit = 0; iit < niterations_max; iit++) {
            g_chi2_full->SetPoint(iit, iit + 1, ndf_full_list[iit] > 0 ? chi2_full_list[iit] / ndf_full_list[iit] : 0);
        }
        g_chi2_full->SetTitle(Form("Full Closure: #chi^{2}/ndf vs Iteration (%.1f < #eta < %.1f);Iteration;#chi^{2}/ndf",
                                   eta_bins[ieta], eta_bins[ieta + 1]));
        g_chi2_full->SetMarkerStyle(20);
        g_chi2_full->SetMarkerColor(kBlue);
        g_chi2_full->SetLineColor(kBlue);
        g_chi2_full->Draw("APL");
        TLine *line1 = new TLine(0, 1.0, niterations_max + 1, 1.0);
        line1->SetLineColor(kRed);
        line1->SetLineStyle(2);
        line1->Draw();
        c_chi2_full->SaveAs(Form("figure/FullClosure_chi2_eta%d.pdf", ieta));
        c_chi2_full->Write();
        delete c_chi2_full;

        // ============================================
        // HALF CLOSURE TEST
        // ============================================
        std::cout << "\n--- HALF CLOSURE TEST ---" << std::endl;

        if (!h_pT_truth_secondhalf_response[ieta] || !h_pT_reco_secondhalf_response[ieta]) {
            std::cerr << "ERROR: Second half histograms not available for eta bin " << ieta << std::endl;
            continue;
        }

        std::vector<TH1D *> h_unfolded_half_list;
        std::vector<double> chi2_half_list;
        std::vector<int> ndf_half_list;

        for (int iit = 0; iit < niterations_max; iit++)
        {
            // Unfold second half using first half response matrix
            RooUnfoldBayes *half_unfold = new RooUnfoldBayes(responses_half[ieta], h_pT_reco_secondhalf_response[ieta], iit + 1);
            TH1D *hRecoPT_half = (TH1D *)half_unfold->Hunfold();
            hRecoPT_half->SetName(Form("h_unfolded_half_eta%d_iter%d", ieta, iit + 1));

            // Calculate chi2 against second half truth
            double chi2 = 0;
            int ndf = 0;
            for (int ibin = 1; ibin <= hRecoPT_half->GetNbinsX(); ibin++) {
                double unf_val = hRecoPT_half->GetBinContent(ibin);
                double unf_err = hRecoPT_half->GetBinError(ibin);
                double tru_val = h_pT_truth_secondhalf_response[ieta]->GetBinContent(ibin);

                if (unf_err > 0 && tru_val > 0) {
                    chi2 += TMath::Power((unf_val - tru_val) / unf_err, 2);
                    ndf++;
                }
            }

            chi2_half_list.push_back(chi2);
            ndf_half_list.push_back(ndf);

            std::cout << "  Iteration " << iit + 1 << ": chi2/ndf = " << chi2 << "/" << ndf
                      << " = " << (ndf > 0 ? chi2 / ndf : 0) << std::endl;

            h_unfolded_half_list.push_back(hRecoPT_half);

            if (save_all_iterations) {
                hRecoPT_half->Write();
            }

            // Plot for key iterations (1, 2, 4)
            if (iit == 0 || iit == 1 || iit == 3)
            {
                std::vector<TH1F *> h_input;
                h_input.push_back((TH1F *)h_pT_truth_secondhalf_response[ieta]->Clone());
                h_input.push_back((TH1F *)h_pT_reco_secondhalf_response[ieta]->Clone());
                h_input.push_back((TH1F *)hRecoPT_half->Clone());

                std::vector<std::string> text;
                std::vector<std::string> legend;

                text.push_back("Half closure test");
                text.push_back(Form("Iteration: %d", iit + 1));
                text.push_back(Form("#chi^{2}/ndf = %.2f/%d = %.2f", chi2, ndf, ndf > 0 ? chi2 / ndf : 0));
                text.push_back(Form("%.1f < #eta < %.1f", eta_bins[ieta], eta_bins[ieta + 1]));
                text.push_back("Response from 1st half");

                legend.push_back("2nd half truth");
                legend.push_back("2nd half reco");
                legend.push_back("Unfolded (2nd half)");

                draw_1D_multiple_plot_ratio(h_input, colors, markerstyle,
                                            false, 10, true,
                                            false, 0, 40, false,
                                            false, 0, 0.5, true,
                                            true, 0.8, 1.2,
                                            true, "p_{T}^{#gamma} [GeV]", "dN_{#gamma}/dp_{T}^{#gamma}", "Unfolded/Truth",
                                            false, "sPHENIX Simulation",
                                            true, text, 0.25, 0.41, 0.04,
                                            true, legend, 0.58, 0.76, 0.04,
                                            Form("figure/HalfClosure_eta%d_iter%d.pdf", ieta, iit + 1));

                for (auto h : h_input) delete h;
            }

            delete half_unfold;
        }

        // Plot chi2 vs iteration for half closure
        TCanvas *c_chi2_half = new TCanvas(Form("c_chi2_half_eta%d", ieta), "Chi2 vs Iteration", 800, 600);
        TGraphErrors *g_chi2_half = new TGraphErrors();
        for (int iit = 0; iit < niterations_max; iit++) {
            g_chi2_half->SetPoint(iit, iit + 1, ndf_half_list[iit] > 0 ? chi2_half_list[iit] / ndf_half_list[iit] : 0);
        }
        g_chi2_half->SetTitle(Form("Half Closure: #chi^{2}/ndf vs Iteration (%.1f < #eta < %.1f);Iteration;#chi^{2}/ndf",
                                   eta_bins[ieta], eta_bins[ieta + 1]));
        g_chi2_half->SetMarkerStyle(20);
        g_chi2_half->SetMarkerColor(kGreen + 2);
        g_chi2_half->SetLineColor(kGreen + 2);
        g_chi2_half->Draw("APL");
        TLine *line2 = new TLine(0, 1.0, niterations_max + 1, 1.0);
        line2->SetLineColor(kRed);
        line2->SetLineStyle(2);
        line2->Draw();
        c_chi2_half->SaveAs(Form("figure/HalfClosure_chi2_eta%d.pdf", ieta));
        c_chi2_half->Write();
        delete c_chi2_half;

        // Comparison plot: Full vs Half closure
        TCanvas *c_compare = new TCanvas(Form("c_compare_eta%d", ieta), "Closure Comparison", 800, 600);
        g_chi2_full->SetTitle(Form("Closure Test Comparison (%.1f < #eta < %.1f);Iteration;#chi^{2}/ndf",
                                   eta_bins[ieta], eta_bins[ieta + 1]));
        g_chi2_full->GetYaxis()->SetRangeUser(0, TMath::Max(5.0, 1.2 * TMath::Max(
            *std::max_element(chi2_full_list.begin(), chi2_full_list.end()) / *std::max_element(ndf_full_list.begin(), ndf_full_list.end()),
            *std::max_element(chi2_half_list.begin(), chi2_half_list.end()) / *std::max_element(ndf_half_list.begin(), ndf_half_list.end())
        )));
        g_chi2_full->Draw("APL");
        g_chi2_half->Draw("PL SAME");
        TLine *line3 = new TLine(0, 1.0, niterations_max + 1, 1.0);
        line3->SetLineColor(kRed);
        line3->SetLineStyle(2);
        line3->Draw();
        TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.85);
        leg->AddEntry(g_chi2_full, "Full Closure", "lp");
        leg->AddEntry(g_chi2_half, "Half Closure", "lp");
        leg->AddEntry(line3, "#chi^{2}/ndf = 1", "l");
        leg->Draw();
        c_compare->SaveAs(Form("figure/ClosureComparison_eta%d.pdf", ieta));
        c_compare->Write();
        delete c_compare;
    }

    fout->Close();
    fresin->Close();

    std::cout << "\n========================================" << std::endl;
    std::cout << "Closure test complete!" << std::endl;
    std::cout << "Output ROOT file: " << output_rootfile << std::endl;
    std::cout << "Plots saved to: figure/" << std::endl;
    std::cout << "========================================\n" << std::endl;

    delete fout;
    delete fresin;
}