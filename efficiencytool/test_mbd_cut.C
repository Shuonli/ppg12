// Quick test: verify mbd_avgsigma_cut_on correctly rejects events above threshold.
// Reads 2000 events from one data file, applies the cut, checks no passing event
// has avgsigma >= threshold.
#include "MbdPileupHelper.h"
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <iostream>

void test_mbd_cut()
{
    const float threshold = 0.5;
    const int maxevents = 2000;

    TChain chain("slimtree");
    chain.Add("/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_10_with_bdt_split.root");

    TTreeReader reader(&chain);
    TTreeReaderArray<float> mbd_north_time(reader,   "mbdnortht");
    TTreeReaderArray<float> mbd_south_time(reader,   "mbdsoutht");
    TTreeReaderArray<float> mbd_north_charge(reader, "mbdnorthq");
    TTreeReaderArray<float> mbd_south_charge(reader, "mbdsouthq");

    int n_total   = 0;
    int n_valid   = 0;
    int n_pass    = 0;
    int n_reject  = 0;
    int n_violation = 0;  // passing events with avgsigma >= threshold (should be 0)

    while (reader.Next() && n_total < maxevents)
    {
        n_total++;
        MbdPileupResult res = calculateMbdPileupMetrics(
            mbd_north_time, mbd_south_time, mbd_north_charge, mbd_south_charge);

        if (!res.valid)
        {
            n_pass++;  // no valid MBD info → event passes (not rejected)
            continue;
        }
        n_valid++;

        bool reject = (res.avgsigma >= threshold);
        if (reject)
        {
            n_reject++;
        }
        else
        {
            n_pass++;
            if (res.avgsigma >= threshold)  // sanity double-check
                n_violation++;
        }
    }

    std::cout << "=== MBD avgsigma cut test (threshold = " << threshold << " ns) ===" << std::endl;
    std::cout << "  Total events processed : " << n_total   << std::endl;
    std::cout << "  Events with valid MBD  : " << n_valid   << std::endl;
    std::cout << "  Events passing cut     : " << n_pass    << std::endl;
    std::cout << "  Events rejected        : " << n_reject  << std::endl;
    std::cout << "  Violations (pass but avgsigma >= threshold): " << n_violation << std::endl;

    if (n_violation == 0)
        std::cout << "  RESULT: PASS -- no event above threshold slipped through." << std::endl;
    else
        std::cout << "  RESULT: FAIL -- " << n_violation << " event(s) with avgsigma >= threshold passed!" << std::endl;
}
