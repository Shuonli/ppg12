#include <yaml-cpp/yaml.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
struct CountABCD
{
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
};

struct RunStats
{
    double sum = 0.0;
    double sumsq = 0.0;
    int count = 0;
};

double safe_mean(double sum, int count)
{
    return count > 0 ? sum / static_cast<double>(count) : 0.0;
}

double safe_err(double sum, double sumsq, int count)
{
    if (count <= 1)
    {
        return 0.0;
    }
    const double mean = sum / static_cast<double>(count);
    const double var = std::max(0.0, sumsq / static_cast<double>(count) - mean * mean);
    return std::sqrt(var / static_cast<double>(count));
}

double compute_std_time(const TTreeReaderArray<float>& t,
                        const TTreeReaderArray<float>& q,
                        float qmin)
{
    const int n = std::min((int)t.GetSize(), (int)q.GetSize());
    if (n <= 0)
    {
        return 0.0;
    }
    double sum = 0.0;
    double sumsq = 0.0;
    int nsel = 0;
    for (int i = 0; i < n; ++i)
    {
        if (q[i] <= qmin)
            continue;
        const double ti = t[i];
        if (!std::isfinite(ti))
            continue;
        sum += ti;
        sumsq += ti * ti;
        ++nsel;
    }
    if (nsel <= 1)
    {
        return 0.0;
    }
    const double mean = sum / static_cast<double>(nsel);
    const double var = std::max(0.0, (sumsq / static_cast<double>(nsel)) - mean * mean);
    return std::sqrt(var);
}

double compute_mbd_avg_sigma(const TTreeReaderArray<float>& north_t,
                             const TTreeReaderArray<float>& north_q,
                             const TTreeReaderArray<float>& south_t,
                             const TTreeReaderArray<float>& south_q)
{
    const double sigma_north = compute_std_time(north_t, north_q, 0.1f);
    const double sigma_south = compute_std_time(south_t, south_q, 0.1f);
    if (sigma_north > 0.0 && sigma_south > 0.0)
    {
        return 0.5 * (sigma_north + sigma_south);
    }
    if (sigma_north > 0.0)
        return sigma_north;
    if (sigma_south > 0.0)
        return sigma_south;
    return 0.0;
}

double abcd_purity(const CountABCD& c)
{
    if (c.A <= 0.0 || c.D <= 0.0)
    {
        return 0.0;
    }
    const double bkg = c.B * c.C / c.D;
    const double pur = (c.A - bkg) / c.A;
    return pur;
}
} // namespace

void NPB_PurityStudy(const std::string& configname = "config_showershape.yaml",
                     const std::string& filetype = "data")
{
    gSystem->Load("/sphenix/u/shuhang98/install/lib64/libyaml-cpp.so");
    YAML::Node configYaml = YAML::LoadFile(configname);

    const bool issim = (filetype != "data");
    std::string infilename_root_dir = configYaml["input"]["photon_jet_file_root_dir"].as<std::string>();
    std::string infilename_branch_dir = configYaml["input"]["photon_jet_file_branch_dir"].as<std::string>();
    std::string infilename = infilename_root_dir + filetype + infilename_branch_dir;
    if (!issim)
    {
        infilename = configYaml["input"]["data_file"].as<std::string>();
    }

    std::cout << "Input file: " << infilename << std::endl;

    const std::string treename = configYaml["input"]["tree"].as<std::string>();
    TChain chain(treename.c_str());
    chain.Add(infilename.c_str());

    const std::string clusternodename = configYaml["input"]["cluster_node_name"].as<std::string>();
    const std::string bdt_model_name = configYaml["input"]["bdt_model_name"].as<std::string>("base");

    const int conesize = configYaml["analysis"]["cone_size"].as<int>();
    const float reco_min_ET = configYaml["analysis"]["reco_min_ET"].as<float>();
    const float recoiso_min = configYaml["analysis"]["reco_iso_min"].as<float>();
    const float recoiso_max_b = configYaml["analysis"]["reco_iso_max_b"].as<float>();
    const float recoiso_max_s = configYaml["analysis"]["reco_iso_max_s"].as<float>();
    const float recononiso_min_shift = configYaml["analysis"]["reco_noniso_min_shift"].as<float>();
    const float recononiso_max = configYaml["analysis"]["reco_noniso_max"].as<float>();
    const float vertexcut = configYaml["analysis"]["vertex_cut"].as<float>();

    const int iso_threshold = configYaml["analysis"]["iso_threshold"].as<int>(0);
    const int iso_hcalonly = configYaml["analysis"]["iso_hcalonly"].as<int>(0);
    const float iso_emcalinnerr = configYaml["analysis"]["iso_emcalinnerr"].as<float>(0.0);

    const float mc_iso_shift = configYaml["analysis"]["mc_iso_shift"].as<float>(0.0);
    const float mc_iso_scale = configYaml["analysis"]["mc_iso_scale"].as<float>(1.2);

    const int n_nt_fail = configYaml["analysis"]["n_nt_fail"].as<int>(1);
    const int weta_fail = configYaml["analysis"]["weta_fail"].as<int>(0);
    const int wphi_fail = configYaml["analysis"]["wphi_fail"].as<int>(0);
    const int e11_to_e33_fail = configYaml["analysis"]["e11_to_e33_fail"].as<int>(0);
    const int e32_to_e35_fail = configYaml["analysis"]["e32_to_e35_fail"].as<int>(0);
    const int et1_fail = configYaml["analysis"]["et1_fail"].as<int>(0);
    const int bdt_fail = configYaml["analysis"]["bdt_fail"].as<int>(0);
    const int weta_on = configYaml["analysis"]["weta_on"].as<int>(1);
    const int wphi_on = configYaml["analysis"]["wphi_on"].as<int>(1);
    const int e11_to_e33_on = configYaml["analysis"]["e11_to_e33_on"].as<int>(1);
    const int e32_to_e35_on = configYaml["analysis"]["e32_to_e35_on"].as<int>(1);
    const int et1_on = configYaml["analysis"]["et1_on"].as<int>(1);
    const int et2_on = configYaml["analysis"]["et2_on"].as<int>(1);
    const int et3_on = configYaml["analysis"]["et3_on"].as<int>(1);
    const int et4_on = configYaml["analysis"]["et4_on"].as<int>(1);
    const int bdt_on = configYaml["analysis"]["bdt_on"].as<int>(1);

    // Tight cuts
    float tight_weta_cogx_min = configYaml["analysis"]["tight"]["weta_cogx_min"].as<float>();
    float tight_weta_cogx_max = configYaml["analysis"]["tight"]["weta_cogx_max"].as<float>();
    float tight_weta_cogx_max_b = configYaml["analysis"]["tight"]["weta_cogx_max_b"].as<float>();
    float tight_weta_cogx_max_s = configYaml["analysis"]["tight"]["weta_cogx_max_s"].as<float>();
    float tight_wphi_cogx_min = configYaml["analysis"]["tight"]["wphi_cogx_min"].as<float>();
    float tight_wphi_cogx_max = configYaml["analysis"]["tight"]["wphi_cogx_max"].as<float>();
    float tight_wphi_cogx_max_b = configYaml["analysis"]["tight"]["wphi_cogx_max_b"].as<float>();
    float tight_wphi_cogx_max_s = configYaml["analysis"]["tight"]["wphi_cogx_max_s"].as<float>();
    float tight_e11_over_e33_min = configYaml["analysis"]["tight"]["e11_over_e33_min"].as<float>();
    float tight_e11_over_e33_max = configYaml["analysis"]["tight"]["e11_over_e33_max"].as<float>();
    float tight_et1_min = configYaml["analysis"]["tight"]["et1_min"].as<float>();
    float tight_et1_max = configYaml["analysis"]["tight"]["et1_max"].as<float>();
    float tight_et1_min_b = configYaml["analysis"]["tight"]["et1_min_b"].as<float>();
    float tight_et1_min_s = configYaml["analysis"]["tight"]["et1_min_s"].as<float>();
    float tight_et2_min = configYaml["analysis"]["tight"]["et2_min"].as<float>(0.0);
    float tight_et2_max = configYaml["analysis"]["tight"]["et2_max"].as<float>(1.0);
    float tight_et3_min = configYaml["analysis"]["tight"]["et3_min"].as<float>(0.0);
    float tight_et3_max = configYaml["analysis"]["tight"]["et3_max"].as<float>(1.0);
    float tight_e32_over_e35_min = configYaml["analysis"]["tight"]["e32_over_e35_min"].as<float>();
    float tight_e32_over_e35_max = configYaml["analysis"]["tight"]["e32_over_e35_max"].as<float>();
    float tight_prob_min = configYaml["analysis"]["tight"]["prob_min"].as<float>();
    float tight_prob_max = configYaml["analysis"]["tight"]["prob_max"].as<float>();
    float tight_et4_min = configYaml["analysis"]["tight"]["et4_min"].as<float>();
    float tight_et4_max = configYaml["analysis"]["tight"]["et4_max"].as<float>();
    float tight_bdt_min = configYaml["analysis"]["tight"]["bdt_min"].as<float>(0.0);
    float tight_bdt_max = configYaml["analysis"]["tight"]["bdt_max"].as<float>(1.0);

    // Non-tight cuts
    float non_tight_weta_cogx_min = configYaml["analysis"]["non_tight"]["weta_cogx_min"].as<float>();
    float non_tight_weta_cogx_max = configYaml["analysis"]["non_tight"]["weta_cogx_max"].as<float>();
    float non_tight_weta_cogx_max_b = configYaml["analysis"]["non_tight"]["weta_cogx_max_b"].as<float>();
    float non_tight_weta_cogx_max_s = configYaml["analysis"]["non_tight"]["weta_cogx_max_s"].as<float>();
    float non_tight_wphi_cogx_min = configYaml["analysis"]["non_tight"]["wphi_cogx_min"].as<float>();
    float non_tight_wphi_cogx_max = configYaml["analysis"]["non_tight"]["wphi_cogx_max"].as<float>();
    float non_tight_wphi_cogx_max_b = configYaml["analysis"]["non_tight"]["wphi_cogx_max_b"].as<float>();
    float non_tight_wphi_cogx_max_s = configYaml["analysis"]["non_tight"]["wphi_cogx_max_s"].as<float>();
    float non_tight_prob_min = configYaml["analysis"]["non_tight"]["prob_min"].as<float>();
    float non_tight_prob_max = configYaml["analysis"]["non_tight"]["prob_max"].as<float>();
    float non_tight_et1_min = configYaml["analysis"]["non_tight"]["et1_min"].as<float>();
    float non_tight_et1_max = configYaml["analysis"]["non_tight"]["et1_max"].as<float>();
    float non_tight_e11_over_e33_min = configYaml["analysis"]["non_tight"]["e11_over_e33_min"].as<float>();
    float non_tight_e11_over_e33_max = configYaml["analysis"]["non_tight"]["e11_over_e33_max"].as<float>();
    float non_tight_e32_over_e35_min = configYaml["analysis"]["non_tight"]["e32_over_e35_min"].as<float>();
    float non_tight_e32_over_e35_max = configYaml["analysis"]["non_tight"]["e32_over_e35_max"].as<float>();
    float non_tight_et4_min = configYaml["analysis"]["non_tight"]["et4_min"].as<float>();
    float non_tight_et4_max = configYaml["analysis"]["non_tight"]["et4_max"].as<float>();
    float non_tight_bdt_min = configYaml["analysis"]["non_tight"]["bdt_min"].as<float>(0.0);
    float non_tight_bdt_max = configYaml["analysis"]["non_tight"]["bdt_max"].as<float>(1.0);

    // Common cuts
    float common_prob_min = configYaml["analysis"]["common"]["prob_min"].as<float>();
    float common_prob_max = configYaml["analysis"]["common"]["prob_max"].as<float>();
    float common_e11_over_e33_min = configYaml["analysis"]["common"]["e11_over_e33_min"].as<float>();
    float common_e11_over_e33_max = configYaml["analysis"]["common"]["e11_over_e33_max"].as<float>();
    float common_wr_cogx_bound = configYaml["analysis"]["common"]["wr_cogx_bound"].as<float>();
    float common_cluster_weta_cogx_bound = configYaml["analysis"]["common"]["cluster_weta_cogx_bound"].as<float>();

    const float npb_score_cut = configYaml["analysis"]["npb_score_cut"].as<float>(0.5);
    const float npb_score_abcd_cut = 0.5f;
    const int mrad_split_run = configYaml["analysis"]["mrad_split_run"].as<int>(0);

    const std::string toy_double_branch = configYaml["analysis"]["toy_double_branch"].as<std::string>("is_double_interaction");
    const std::string toy_signal_branch = configYaml["analysis"]["toy_signal_branch"].as<std::string>("is_signal");

    // Optional per-cluster NPB score branch:
    std::string npb_score_branch = Form("cluster_npb_score_%s", clusternodename.c_str());
    if (!chain.GetBranch(npb_score_branch.c_str()))
    {
        if (chain.GetBranch("cluster_npb_score"))
        {
            npb_score_branch = "cluster_npb_score";
        }
        else
        {
            npb_score_branch.clear();
            std::cerr << "[NPBScore] WARNING: cannot find branch 'cluster_npb_score_<node>' or 'cluster_npb_score' in input tree." << std::endl;
        }
    }

    TTreeReader reader(&chain);
    TTreeReaderValue<int> ncluster(reader, Form("ncluster_%s", clusternodename.c_str()));
    TTreeReaderValue<float> vertexz(reader, "vertexz");
    TTreeReaderValue<int> runnumber(reader, "runnumber");
    TTreeReaderArray<float> mbdnortht(reader, "mbdnortht");
    TTreeReaderArray<float> mbdsoutht(reader, "mbdsoutht");
    TTreeReaderArray<float> mbdnorthq(reader, "mbdnorthq");
    TTreeReaderArray<float> mbdsouthq(reader, "mbdsouthq");

    TTreeReaderArray<float> cluster_E(reader, Form("cluster_E_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Et(reader, Form("cluster_Et_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Eta(reader, Form("cluster_Eta_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_Phi(reader, Form("cluster_Phi_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_prob(reader, Form("cluster_prob_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et1(reader, Form("cluster_et1_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et2(reader, Form("cluster_et2_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et3(reader, Form("cluster_et3_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_et4(reader, Form("cluster_et4_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_weta_cogx(reader, Form("cluster_weta_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_wphi_cogx(reader, Form("cluster_wphi_cogx_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e11(reader, Form("cluster_e11_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e33(reader, Form("cluster_e33_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e32(reader, Form("cluster_e32_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_e35(reader, Form("cluster_e35_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02(reader, Form("cluster_iso_02_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03(reader, Form("cluster_iso_03_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_04(reader, Form("cluster_iso_04_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_emcal(reader, Form("cluster_iso_03_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalin(reader, Form("cluster_iso_03_70_hcalin_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_03_70_hcalout(reader, Form("cluster_iso_03_70_hcalout_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_02_70_emcal(reader, Form("cluster_iso_02_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_005_70_emcal(reader, Form("cluster_iso_005_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_iso_01_70_emcal(reader, Form("cluster_iso_01_70_emcal_%s", clusternodename.c_str()));
    TTreeReaderArray<float> cluster_bdt(reader, Form("cluster_bdt_%s_%s", clusternodename.c_str(), bdt_model_name.c_str()));

    std::unique_ptr<TTreeReaderArray<float>> cluster_npb_score;
    if (!npb_score_branch.empty())
    {
        cluster_npb_score = std::make_unique<TTreeReaderArray<float>>(reader, npb_score_branch.c_str());
    }

    std::unique_ptr<TTreeReaderValue<int>> toy_double;
    std::unique_ptr<TTreeReaderValue<int>> toy_signal;
    if (chain.GetBranch(toy_double_branch.c_str()))
    {
        toy_double = std::make_unique<TTreeReaderValue<int>>(reader, toy_double_branch.c_str());
    }
    if (chain.GetBranch(toy_signal_branch.c_str()))
    {
        toy_signal = std::make_unique<TTreeReaderValue<int>>(reader, toy_signal_branch.c_str());
    }

    CountABCD total_counts;

    double toy_double_signal = 0.0;
    double toy_double_total = 0.0;
    double toy_single_signal = 0.0;
    double toy_single_total = 0.0;

    std::map<int, RunStats> run_stats;

    const int nentries = chain.GetEntries();
    int ientry = 0;
    while (reader.Next())
    {
        if (ientry % 10000 == 0)
        {
            std::cout << "Processing entry " << ientry << " / " << nentries << std::endl;
        }
        ++ientry;

        if (std::abs(*vertexz) > vertexcut)
        {
            continue;
        }

        const double mbd_avg_sigma = compute_mbd_avg_sigma(mbdnortht, mbdnorthq, mbdsoutht, mbdsouthq);
        auto& stats = run_stats[*runnumber];
        stats.sum += mbd_avg_sigma;
        stats.sumsq += mbd_avg_sigma * mbd_avg_sigma;
        stats.count += 1;

        for (int icluster = 0; icluster < *ncluster; ++icluster)
        {
            const float clusterET = cluster_Et[icluster];
            if (clusterET < reco_min_ET)
                continue;

            if (!std::isfinite(cluster_Eta[icluster]))
                continue;

            const float recoiso_max = recoiso_max_b + recoiso_max_s * clusterET;
            float recoisoET = -999.0f;
            if (conesize == 4)
            {
                recoisoET = cluster_iso_04[icluster];
            }
            else if (conesize == 3)
            {
                recoisoET = cluster_iso_03[icluster];
            }
            else if (conesize == 2)
            {
                recoisoET = cluster_iso_02[icluster];
            }
            else
            {
                std::cout << "Error: conesize not supported" << std::endl;
                continue;
            }

            if (iso_threshold)
            {
                if (iso_hcalonly)
                {
                    recoisoET = cluster_iso_03_70_hcalin[icluster] + cluster_iso_03_70_hcalout[icluster];
                }
                else
                {
                    recoisoET = cluster_iso_03_70_emcal[icluster] +
                                cluster_iso_03_70_hcalin[icluster] +
                                cluster_iso_03_70_hcalout[icluster];
                    if (iso_emcalinnerr > 0.04 && iso_emcalinnerr < 0.06)
                    {
                        recoisoET -= cluster_iso_005_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.09 && iso_emcalinnerr < 0.11)
                    {
                        recoisoET -= cluster_iso_01_70_emcal[icluster];
                    }
                    else if (iso_emcalinnerr > 0.19 && iso_emcalinnerr < 0.21)
                    {
                        recoisoET -= cluster_iso_02_70_emcal[icluster];
                    }
                }
            }

            if (issim)
            {
                recoisoET = recoisoET * mc_iso_scale;
                recoisoET += mc_iso_shift;
            }

            bool iso = (recoisoET > recoiso_min && recoisoET < recoiso_max);
            const float recononiso_min = recoiso_max + recononiso_min_shift;
            bool noniso = (recoisoET > recononiso_min && recoisoET < recononiso_max);

            const float e11_over_e33 = cluster_e33[icluster] > 0 ? cluster_e11[icluster] / cluster_e33[icluster] : 0.0f;
            const float e32_over_e35 = cluster_e35[icluster] > 0 ? cluster_e32[icluster] / cluster_e35[icluster] : 0.0f;
            const float wr_cogx = cluster_weta_cogx[icluster] > 0 ? cluster_wphi_cogx[icluster] / cluster_weta_cogx[icluster] : 0.0f;

            bool common_pass =
                cluster_prob[icluster] > common_prob_min &&
                cluster_prob[icluster] < common_prob_max &&
                e11_over_e33 > common_e11_over_e33_min &&
                e11_over_e33 < common_e11_over_e33_max &&
                wr_cogx < common_wr_cogx_bound &&
                cluster_weta_cogx[icluster] < common_cluster_weta_cogx_bound;

            bool tight = false;
            bool nontight = false;

            if (common_pass)
            {
                tight_weta_cogx_max = tight_weta_cogx_max_b + tight_weta_cogx_max_s * clusterET;
                tight_wphi_cogx_max = tight_wphi_cogx_max_b + tight_wphi_cogx_max_s * clusterET;
                tight_et1_min = tight_et1_min_b + tight_et1_min_s * clusterET;

                const bool is_cluster_weta_cogx_tight =
                    (cluster_weta_cogx[icluster] > tight_weta_cogx_min) &&
                    (cluster_weta_cogx[icluster] < tight_weta_cogx_max);
                const bool is_cluster_wphi_cogx_tight =
                    (cluster_wphi_cogx[icluster] > tight_wphi_cogx_min) &&
                    (cluster_wphi_cogx[icluster] < tight_wphi_cogx_max);
                const bool is_cluster_et1_tight =
                    (cluster_et1[icluster] > tight_et1_min) &&
                    (cluster_et1[icluster] < tight_et1_max);
                const bool is_cluster_et2_tight =
                    (cluster_et2[icluster] > tight_et2_min) &&
                    (cluster_et2[icluster] < tight_et2_max);
                const bool is_cluster_et3_tight =
                    (cluster_et3[icluster] > tight_et3_min) &&
                    (cluster_et3[icluster] < tight_et3_max);
                const bool is_e11_over_e33_tight =
                    (e11_over_e33 > tight_e11_over_e33_min) &&
                    (e11_over_e33 < tight_e11_over_e33_max);
                const bool is_e32_over_e35_tight =
                    (e32_over_e35 > tight_e32_over_e35_min) &&
                    (e32_over_e35 < tight_e32_over_e35_max);
                const bool is_cluster_et4_tight =
                    (cluster_et4[icluster] > tight_et4_min) &&
                    (cluster_et4[icluster] < tight_et4_max);
                const bool is_cluster_prob_tight =
                    (cluster_prob[icluster] > tight_prob_min) &&
                    (cluster_prob[icluster] < tight_prob_max);
                const bool is_bdt_tight =
                    (cluster_bdt[icluster] > tight_bdt_min) &&
                    (cluster_bdt[icluster] < tight_bdt_max);

                if (is_cluster_weta_cogx_tight &&
                    is_cluster_wphi_cogx_tight &&
                    is_cluster_et1_tight &&
                    is_cluster_et2_tight &&
                    is_cluster_et3_tight &&
                    is_e11_over_e33_tight &&
                    is_e32_over_e35_tight &&
                    is_cluster_et4_tight &&
                    is_cluster_prob_tight &&
                    is_bdt_tight)
                {
                    tight = true;
                }

                if (cluster_weta_cogx[icluster] > non_tight_weta_cogx_min &&
                    cluster_weta_cogx[icluster] < non_tight_weta_cogx_max &&
                    cluster_wphi_cogx[icluster] > non_tight_wphi_cogx_min &&
                    cluster_wphi_cogx[icluster] < non_tight_wphi_cogx_max &&
                    cluster_prob[icluster] > non_tight_prob_min &&
                    cluster_prob[icluster] < non_tight_prob_max &&
                    e11_over_e33 > non_tight_e11_over_e33_min &&
                    e11_over_e33 < non_tight_e11_over_e33_max &&
                    e32_over_e35 > non_tight_e32_over_e35_min &&
                    e32_over_e35 < non_tight_e32_over_e35_max &&
                    cluster_et1[icluster] > non_tight_et1_min &&
                    cluster_et1[icluster] < non_tight_et1_max &&
                    cluster_et4[icluster] > non_tight_et4_min &&
                    cluster_et4[icluster] < non_tight_et4_max &&
                    cluster_bdt[icluster] > non_tight_bdt_min &&
                    cluster_bdt[icluster] < non_tight_bdt_max)
                {
                    int nfail = 0;
                    if (!is_cluster_weta_cogx_tight) nfail += weta_on;
                    if (!is_cluster_wphi_cogx_tight) nfail += wphi_on;
                    if (!is_cluster_et1_tight) nfail += et1_on;
                    if (!is_cluster_et2_tight) nfail += et2_on;
                    if (!is_cluster_et3_tight) nfail += et3_on;
                    if (!is_e11_over_e33_tight) nfail += e11_to_e33_on;
                    if (!is_e32_over_e35_tight) nfail += e32_to_e35_on;
                    if (!is_cluster_et4_tight) nfail += et4_on;
                    if (!is_cluster_prob_tight) nfail += 1;
                    if (!is_bdt_tight) nfail += bdt_on;

                    if (nfail > n_nt_fail)
                    {
                        bool all_flags_fail = true;
                        if (weta_fail && is_cluster_weta_cogx_tight) all_flags_fail = false;
                        if (wphi_fail && is_cluster_wphi_cogx_tight) all_flags_fail = false;
                        if (et1_fail && is_cluster_et1_tight) all_flags_fail = false;
                        if (e11_to_e33_fail && is_e11_over_e33_tight) all_flags_fail = false;
                        if (e32_to_e35_fail && is_e32_over_e35_tight) all_flags_fail = false;
                        if (bdt_fail && is_bdt_tight) all_flags_fail = false;
                        if (all_flags_fail)
                        {
                            nontight = true;
                        }
                    }
                }
            }

            const float npb_score_val = (cluster_npb_score ? (*cluster_npb_score)[icluster] : 1.0f);
            if (npb_score_val >= npb_score_abcd_cut)
            {
                continue;
            }

            if (tight && iso) total_counts.A += 1.0;
            if (tight && noniso) total_counts.B += 1.0;
            if (nontight && iso) total_counts.C += 1.0;
            if (nontight && noniso) total_counts.D += 1.0;

            if (toy_double && toy_signal && tight && iso)
            {
                if (**toy_double)
                {
                    toy_double_total += 1.0;
                    if (**toy_signal)
                    {
                        toy_double_signal += 1.0;
                    }
                }
                else
                {
                    toy_single_total += 1.0;
                    if (**toy_signal)
                    {
                        toy_single_signal += 1.0;
                    }
                }
            }
        }
    }

    const std::string outfilename = "results/npb_purity_study_" + filetype + ".root";
    TFile* fout = new TFile(outfilename.c_str(), "RECREATE");

    std::vector<double> runs_0mrad;
    std::vector<double> sigma_0mrad;
    std::vector<double> err_0mrad;
    std::vector<double> runs_15mrad;
    std::vector<double> sigma_15mrad;
    std::vector<double> err_15mrad;

    for (const auto& entry : run_stats)
    {
        const int run = entry.first;
        const RunStats& st = entry.second;
        const double mean = safe_mean(st.sum, st.count);
        const double err = safe_err(st.sum, st.sumsq, st.count);
        if (mrad_split_run > 0 && run < mrad_split_run)
        {
            runs_0mrad.push_back(run);
            sigma_0mrad.push_back(mean);
            err_0mrad.push_back(err);
        }
        else
        {
            runs_15mrad.push_back(run);
            sigma_15mrad.push_back(mean);
            err_15mrad.push_back(err);
        }
    }

    TGraphErrors* gr_0mrad = new TGraphErrors(
        runs_0mrad.size(),
        runs_0mrad.empty() ? nullptr : runs_0mrad.data(),
        sigma_0mrad.empty() ? nullptr : sigma_0mrad.data(),
        nullptr,
        err_0mrad.empty() ? nullptr : err_0mrad.data());
    gr_0mrad->SetName("gr_mbd_avg_sigma_vs_run_0mrad");
    gr_0mrad->SetTitle("MBD Avg #sigma vs Run (0 mrad);Run;MBD Avg #sigma [ns]");
    gr_0mrad->SetMarkerStyle(20);

    TGraphErrors* gr_15mrad = new TGraphErrors(
        runs_15mrad.size(),
        runs_15mrad.empty() ? nullptr : runs_15mrad.data(),
        sigma_15mrad.empty() ? nullptr : sigma_15mrad.data(),
        nullptr,
        err_15mrad.empty() ? nullptr : err_15mrad.data());
    gr_15mrad->SetName("gr_mbd_avg_sigma_vs_run_15mrad");
    gr_15mrad->SetTitle("MBD Avg #sigma vs Run (1.5 mrad);Run;MBD Avg #sigma [ns]");
    gr_15mrad->SetMarkerStyle(21);
    gr_15mrad->SetMarkerColor(kRed + 1);
    gr_15mrad->SetLineColor(kRed + 1);

    TCanvas* c1 = new TCanvas("c_mbd_avg_sigma_vs_run", "MBD Avg Sigma vs Run", 1200, 600);
    c1->SetGrid();
    if (!runs_0mrad.empty())
    {
        gr_0mrad->Draw("AP");
        if (!runs_15mrad.empty())
        {
            gr_15mrad->Draw("P SAME");
        }
    }
    else if (!runs_15mrad.empty())
    {
        gr_15mrad->Draw("AP");
    }

    TLegend* leg = new TLegend(0.12, 0.75, 0.4, 0.88);
    if (!runs_0mrad.empty())
        leg->AddEntry(gr_0mrad, "0 mrad", "p");
    if (!runs_15mrad.empty())
        leg->AddEntry(gr_15mrad, "1.5 mrad", "p");
    leg->Draw();

    gr_0mrad->Write();
    gr_15mrad->Write();
    c1->Write();

    fout->Close();

    std::cout << "\n=== ABCD yields (NPB BDT < 0.5) ===" << std::endl;
    std::cout << "A (tight+iso): " << total_counts.A << std::endl;
    std::cout << "B (tight+noniso): " << total_counts.B << std::endl;
    std::cout << "C (nontight+iso): " << total_counts.C << std::endl;
    std::cout << "D (nontight+noniso): " << total_counts.D << std::endl;

    if (toy_double && toy_signal)
    {
        const double pur_double = (toy_double_total > 0.0) ? (toy_double_signal / toy_double_total) : 0.0;
        const double pur_single = (toy_single_total > 0.0) ? (toy_single_signal / toy_single_total) : 0.0;
        std::cout << "\n=== Toy MC purity (tight+iso, NPB BDT < 0.5) ===" << std::endl;
        std::cout << "Double-interaction purity: " << pur_double << std::endl;
        std::cout << "Single-interaction purity: " << pur_single << std::endl;
    }

    std::cout << "Output written to " << outfilename << std::endl;
}
