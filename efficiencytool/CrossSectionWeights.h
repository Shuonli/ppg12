#ifndef CROSSSECTIONWEIGHTS_H
#define CROSSSECTIONWEIGHTS_H

// Single source of truth for MC sample cross-section weights (pb).
// Authoritative reference: RecoEffCalculator_TTreeReader.C values.
// Include this header instead of copy-pasting constants per file.

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

}  // namespace PPG12

#endif  // CROSSSECTIONWEIGHTS_H
