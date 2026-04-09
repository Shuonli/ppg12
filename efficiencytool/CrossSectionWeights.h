#ifndef CROSSSECTIONWEIGHTS_H
#define CROSSSECTIONWEIGHTS_H

// Single source of truth for MC sample cross-section weights (pb) and
// per-sample kinematic windows.
// Include this header instead of copy-pasting constants per file.

#include <string>

namespace PPG12 {

// Photon-jet samples
constexpr float photon5cross  = 146359.3f;   // 0-14 GeV truth photon pT window
constexpr float photon10cross = 6944.675f;    // 14-30 GeV
constexpr float photon20cross = 130.4461f;    // 30+ GeV

// QCD jet samples
constexpr float jet5cross  = 1.3878e+08f;
constexpr float jet8cross  = 1.15e+07f;
constexpr float jet10cross = 3.997e+06f;
constexpr float jet12cross = 1.4903e+06f;
constexpr float jet15cross = 4.073e+05f;
constexpr float jet20cross = 6.2623e+04f;
constexpr float jet30cross = 2.5298e+03f;
constexpr float jet40cross = 1.3553e+02f;
constexpr float jet50cross = 7.3113f;

// Number of generated events per sample (default)
constexpr float default_nsimevents = 1e7f;

// ---------------------------------------------------------------------------
// SampleConfig: per-filetype kinematic windows and cross-section weight.
// Replaces the ~120-line if/else chains duplicated across analysis macros.
// ---------------------------------------------------------------------------
struct SampleConfig {
    float photon_pt_lower = 0;     // truth photon pT window lower bound
    float photon_pt_upper = 0;     // truth photon pT window upper bound
    float jet_pt_lower    = 0;     // truth jet pT window lower bound
    float jet_pt_upper    = 100;   // truth jet pT window upper bound
    float cluster_ET_upper = 100;  // reco cluster ET upper bound
    float weight          = 1.0f;  // cross-section weight (relative to reference sample)
    bool  isbackground    = false; // true for QCD jet samples
    bool  valid           = false; // false if filetype was not recognized
};

// Return kinematic windows and weight for a given filetype string.
// Photon weights are normalized to photon20cross; jet weights to jet50cross.
inline SampleConfig GetSampleConfig(const std::string &filetype)
{
    SampleConfig c;
    c.valid = true;

    // --- Photon samples ---
    if (filetype == "photon5") {
        c.photon_pt_lower = 0;   c.photon_pt_upper = 14;
        c.weight = photon5cross / photon20cross;
    }
    else if (filetype == "photon10") {
        c.photon_pt_lower = 14;  c.photon_pt_upper = 30;
        c.weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon10_double" || filetype == "photon10_nom") {
        c.photon_pt_lower = 10;  c.photon_pt_upper = 100;
        c.weight = photon10cross / photon20cross;
    }
    else if (filetype == "photon20") {
        c.photon_pt_lower = 30;  c.photon_pt_upper = 200;
        c.weight = 1.0f;
    }
    // --- Jet samples ---
    else if (filetype == "jet5") {
        c.jet_pt_lower = 7;   c.jet_pt_upper = 9;   c.cluster_ET_upper = 10;
        c.weight = jet5cross / jet50cross;   c.isbackground = true;
    }
    else if (filetype == "jet8") {
        c.jet_pt_lower = 9;   c.jet_pt_upper = 14;  c.cluster_ET_upper = 15;
        c.weight = jet8cross / jet50cross;   c.isbackground = true;
    }
    else if (filetype == "jet10") {
        c.jet_pt_lower = 10;  c.jet_pt_upper = 15;  c.cluster_ET_upper = 18;
        c.weight = jet10cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet12") {
        c.jet_pt_lower = 14;  c.jet_pt_upper = 21;  c.cluster_ET_upper = 23;
        c.weight = jet12cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet12_double" || filetype == "jet12_nom") {
        c.jet_pt_lower = 10;  c.jet_pt_upper = 100; c.cluster_ET_upper = 100;
        c.weight = jet12cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet15") {
        c.jet_pt_lower = 15;  c.jet_pt_upper = 21;  c.cluster_ET_upper = 23;
        c.weight = jet15cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet20") {
        c.jet_pt_lower = 21;  c.jet_pt_upper = 32;  c.cluster_ET_upper = 35;
        c.weight = jet20cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet30") {
        c.jet_pt_lower = 32;  c.jet_pt_upper = 42;  c.cluster_ET_upper = 45;
        c.weight = jet30cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet40") {
        c.jet_pt_lower = 42;  c.jet_pt_upper = 100; c.cluster_ET_upper = 100;
        c.weight = jet40cross / jet50cross;  c.isbackground = true;
    }
    else if (filetype == "jet50") {
        c.jet_pt_lower = 52;  c.jet_pt_upper = 100;
        c.weight = jet50cross / jet50cross;  c.isbackground = true;
    }
    else {
        c.valid = false;
    }
    return c;
}

}  // namespace PPG12

#endif  // CROSSSECTIONWEIGHTS_H
