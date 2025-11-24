#include "plotcommon.h"
#include <algorithm>

void plot_pid()
{
    init_plot();
    int etmin = 10; // GeV
    int etmax = 35; // GeV

    TFile *fin = new TFile("/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_nom.root", "READ");

    TH2F *h_tight_iso_pid_pt_0 = (TH2F *)fin->Get("h_tight_iso_pid_pt_0");

    // Step 1: Limit pT range
    h_tight_iso_pid_pt_0->GetYaxis()->SetRangeUser(etmin, etmax);

    // Step 2: Project onto PID axis (X-axis)
    TH1D *h_pid = h_tight_iso_pid_pt_0->ProjectionX("h_pid_proj");

    // Step 3: Create a new histogram for selected PIDs
    vector<int> pid_list = {-211, 211, -11, 11, 111, 22, 221};
    int n_bins = pid_list.size() + 1; // extra bin for "others"
    TH1D *h_pid_selected = new TH1D("h_pid_selected", ";Particle;Counts", n_bins, 0, n_bins);
    set<int> pid_set(pid_list.begin(), pid_list.end()); // for fast lookup
    double others_content = 0;
    double others_error2 = 0; // sum errors in quadrature

    std::map<int, std::pair<double, double>> pid_counts;

    int nbins_x = h_pid->GetNbinsX();
    for (int i = 1; i <= nbins_x; ++i)
    {
        int pid_bin_center = static_cast<int>(h_pid->GetXaxis()->GetBinCenter(i));
        // plus one for negative bin center
        if (pid_bin_center < 0)
            pid_bin_center -= 1;

        double content = h_pid->GetBinContent(i);
        double error = h_pid->GetBinError(i);

        pid_counts[pid_bin_center].first += content;
        pid_counts[pid_bin_center].second = sqrt(pow(pid_counts[pid_bin_center].second, 2) + error * error);

        auto it = find(pid_list.begin(), pid_list.end(), pid_bin_center);
        if (it != pid_list.end())
        {
            int index = distance(pid_list.begin(), it);
            h_pid_selected->SetBinContent(index + 1, content);
            h_pid_selected->SetBinError(index + 1, error);

            // Set bin label
            string label;
            switch (pid_bin_center)
            {
            case -211:
                label = "#pi^{-}";
                break;
            case 211:
                label = "#pi^{+}";
                break;
            case -11:
                label = "e^{-}";
                break;
            case 11:
                label = "e^{+}";
                break;
            case 111:
                label = "#pi^{0}";
                break;
            case 22:
                label = "#gamma";
                break;
            case 221:
                label = "#eta";
                break;
            }
            h_pid_selected->GetXaxis()->SetBinLabel(index + 1, label.c_str());
        }
        else
        {
            others_content += content;
            others_error2 += error * error;
        }
    }

    // Add the "others" bin
    int others_bin = n_bins;
    h_pid_selected->SetBinContent(others_bin, others_content);
    h_pid_selected->SetBinError(others_bin, sqrt(others_error2));
    h_pid_selected->GetYaxis()->SetRangeUser(0, h_pid_selected->GetMaximum() * 1.2); // Adjust Y-axis range
    h_pid_selected->GetXaxis()->SetBinLabel(others_bin, "other");

    // Optional: Set Y axis title and styling
    h_pid_selected->GetYaxis()->SetTitle("Counts");

    // Step 4: Draw
    TCanvas *c = new TCanvas("c_pid", "Selected PID Projection", 800, 600);
    h_pid_selected->SetMarkerStyle(20);
    h_pid_selected->Draw("E1");

    float xpos(0.2), xpos2(0.835), ypos(0.885), ypos2(0.65), dy(0.054), dy1(0.08), fontsize(0.046), fontsize1(0.048);
    myText(xpos2, ypos - 0 * dy, 1, strIncMC.c_str(), fontsize, 1);
    myText(xpos2, ypos - 2 * dy, 1, "tight iso cluster", fontsize, 1);
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2.c_str(), fontsize, 0);

    double total_count = 0;
    for (const auto &[pid, pair] : pid_counts)
    {
        total_count += pair.first;
    }

    // Print PID fractions
    std::cout << "\n--- PID Fractions ---\n";
    for (const auto &[pid, pair] : pid_counts)
    {
        double fraction = (total_count > 0) ? (pair.first / total_count) : 0;
        if(fraction == 0) continue; // Skip if fraction is zero
        printf("PID %5d : Count = %10.2f ± %.2f (%.3f%%)\n", pid, pair.first, pair.second, 100 * fraction);
    }
}