#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <limits>

// Lumi file format (space-separated, header lines start with non-digit):
//   RN  Bit10Corr  Bit10UC  Bit18Corr  Bit18UC  Bit22Corr  Bit22UC  Bit30Corr  Bit30UC  (units pb^-1)
// Pileup file format: runnumber(float) pileup_rate
//
// Returns total Bit30Corr luminosity for runs in [run_min, run_max] (inclusive).
// Also prints lumi-weighted pileup rate for each trigger column.

struct LumiEntry {
    float bit10corr, bit10uc;
    float bit18corr, bit18uc;
    float bit22corr, bit22uc;
    float bit30corr, bit30uc;
};

float LumiCalculator(
    //int run_min = 47289, int run_max = 51274,
    int run_min = 51274, int run_max = 54000,
    //int run_min = -1, int run_max = -1,
    const std::string &lumifilename = "/sphenix/user/shuhangli/ppg12/lumi/60cmLumi_fromJoey.list",
    const std::string &pileupfilename = "/sphenix/user/shuhangli/ppg12/efficiencytool/pileup.list")
{
    std::set<int> skiprunnumbers = {0}; //51274

    // --- Load lumi file ---
    std::map<int, LumiEntry> lumimap;
    {
        std::ifstream lumifile(lumifilename);
        if (!lumifile.is_open())
        {
            std::cerr << "ERROR: cannot open lumi file: " << lumifilename << std::endl;
            return -1.0f;
        }
        std::string line;
        while (std::getline(lumifile, line))
        {
            if (line.empty()) continue;
            if (!std::isdigit(static_cast<unsigned char>(line[0]))) continue;
            std::istringstream iss(line);
            int runnumber;
            LumiEntry e;
            if (!(iss >> runnumber
                      >> e.bit10corr >> e.bit10uc
                      >> e.bit18corr >> e.bit18uc
                      >> e.bit22corr >> e.bit22uc
                      >> e.bit30corr >> e.bit30uc))
                continue;
            lumimap[runnumber] = e;
        }
    }

    // --- Load pileup file ---
    std::map<int, float> pileupmap;
    {
        std::ifstream pileupfile(pileupfilename);
        if (!pileupfile.is_open())
        {
            std::cerr << "WARNING: cannot open pileup file: " << pileupfilename << std::endl;
        }
        else
        {
            std::string line;
            while (std::getline(pileupfile, line))
            {
                if (line.empty()) continue;
                std::istringstream iss(line);
                float run_f, pileup;
                if (!(iss >> run_f >> pileup)) continue;
                pileupmap[static_cast<int>(run_f)] = pileup;
            }
        }
    }

    int effective_min = (run_min > 0) ? run_min : 0;
    int effective_max = (run_max > 0) ? run_max : std::numeric_limits<int>::max();

    // Lumi accumulators
    float tot10c = 0, tot10uc = 0;
    float tot18c = 0, tot18uc = 0;
    float tot22c = 0, tot22uc = 0;
    float tot30c = 0, tot30uc = 0;

    // Lumi-weighted pileup accumulators: sum(lumi_i * pileup_i)
    float wp10c = 0, wp10uc = 0;
    float wp18c = 0, wp18uc = 0;
    float wp22c = 0, wp22uc = 0;
    float wp30c = 0, wp30uc = 0;

    int nruns = 0;

    std::cout << "=== Lumi Calculator ===" << std::endl;
    std::cout << "Lumi file:   " << lumifilename << std::endl;
    std::cout << "Pileup file: " << pileupfilename << std::endl;
    std::cout << "Run range: [" << effective_min << ", " << effective_max << "]" << std::endl;
    std::cout << std::endl;

    // Per-run table
    std::string sep(127, '-');
    printf("%-8s  %10s %10s  %10s %10s  %10s %10s  %10s %10s  %8s\n",
           "Run",
           "Bit10Corr", "Bit10UC",
           "Bit18Corr", "Bit18UC",
           "Bit22Corr", "Bit22UC",
           "Bit30Corr", "Bit30UC",
           "Pileup");
    std::cout << sep << std::endl;

    for (auto const &kv : lumimap)
    {
        int run = kv.first;
        const LumiEntry &e = kv.second;

        if (!(run >= effective_min && run <= effective_max)) continue;
        if (skiprunnumbers.count(run)) continue;
        if (e.bit30corr <= 0 && e.bit18corr <= 0) continue; // skip dead runs

        float pileup = 0.0f;
        bool has_pileup = pileupmap.count(run);
        if (has_pileup) pileup = pileupmap[run];

        printf("%-8d  %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f  %8.5f%s\n",
               run,
               e.bit10corr, e.bit10uc,
               e.bit18corr, e.bit18uc,
               e.bit22corr, e.bit22uc,
               e.bit30corr, e.bit30uc,
               pileup, has_pileup ? "" : " (no pileup)");

        tot10c  += e.bit10corr;  tot10uc  += e.bit10uc;
        tot18c  += e.bit18corr;  tot18uc  += e.bit18uc;
        tot22c  += e.bit22corr;  tot22uc  += e.bit22uc;
        tot30c  += e.bit30corr;  tot30uc  += e.bit30uc;

        if (has_pileup)
        {
            wp10c  += e.bit10corr * pileup;  wp10uc  += e.bit10uc * pileup;
            wp18c  += e.bit18corr * pileup;  wp18uc  += e.bit18uc * pileup;
            wp22c  += e.bit22corr * pileup;  wp22uc  += e.bit22uc * pileup;
            wp30c  += e.bit30corr * pileup;  wp30uc  += e.bit30uc * pileup;
        }

        nruns++;
    }

    std::cout << sep << std::endl;
    printf("%-8s  %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f  %10.4f %10.4f\n",
           "TOTAL",
           tot10c, tot10uc,
           tot18c, tot18uc,
           tot22c, tot22uc,
           tot30c, tot30uc);

    std::cout << std::endl;
    std::cout << "=== Lumi-weighted pileup rate ===" << std::endl;
    printf("%-12s  %10s %10s  %10s %10s  %10s %10s  %10s %10s\n",
           "",
           "Bit10Corr", "Bit10UC",
           "Bit18Corr", "Bit18UC",
           "Bit22Corr", "Bit22UC",
           "Bit30Corr", "Bit30UC");
    printf("%-12s  %10.5f %10.5f  %10.5f %10.5f  %10.5f %10.5f  %10.5f %10.5f\n",
           "<pileup>",
           (tot10c  > 0 ? wp10c  / tot10c  : 0),
           (tot10uc > 0 ? wp10uc / tot10uc : 0),
           (tot18c  > 0 ? wp18c  / tot18c  : 0),
           (tot18uc > 0 ? wp18uc / tot18uc : 0),
           (tot22c  > 0 ? wp22c  / tot22c  : 0),
           (tot22uc > 0 ? wp22uc / tot22uc : 0),
           (tot30c  > 0 ? wp30c  / tot30c  : 0),
           (tot30uc > 0 ? wp30uc / tot30uc : 0));

    std::cout << "Runs counted: " << nruns << std::endl;
    std::cout << "=================================" << std::endl;

    return tot30c; // return Bit30Corr as default
}
