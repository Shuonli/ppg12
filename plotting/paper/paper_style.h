// paper_style.h -- shared header for the PPG12 final-paper figure macros.
//
// Re-exports plotcommon.h (frames, ptRanges, legend strings, init_plot()) and
// adds two paper-specific pieces:
//
//   1. paperFiguresDir / paper_savepath()   -- canonical output location so
//      every plot_paper_*.C macro writes directly into PPG12-Paper/figures/.
//   2. paper_init(bool internal = true)     -- wraps init_plot() and
//      conditionally strips the "Internal" qualifier from strleg1 once the
//      paper goes from draft to publication.
//
// The journal lumi/legend strings already exist in plotcommon.h as
//   strleg2_1  = "#it{p}+#it{p} ... 64.37 pb^{-1}"
//   strleg5    = "         = 64.37 pb^{-1}"
// so paper macros just pick between strleg2 (no lumi) and strleg2_1 (with
// lumi) depending on the figure layout. No further override is needed here.

#ifndef PAPER_STYLE_H
#define PAPER_STYLE_H

#include "../plotcommon.h"

#include <string>

// Output directory for every paper PDF. Relative to plotting/paper/ at exec
// time; the make_paper_figures.sh driver always cd's there before running.
inline const std::string &paper_savepath()
{
    static const std::string p = "../../PPG12-Paper/figures";
    return p;
}

// Default analysis tune used everywhere in the paper figures. Centralised
// here so a single constant flip rebuilds every figure when the nominal
// changes.
inline const std::string &paper_tune()
{
    static const std::string t = "bdt_nom";
    return t;
}

// init_plot() then optionally drop the "Internal" tag from strleg1.
// keep_internal = false flips strleg1 to the publication form.
inline void paper_init(bool keep_internal = true)
{
    init_plot();
    if (!keep_internal) strleg1 = "#bf{#it{sPHENIX}}";
}

#endif // PAPER_STYLE_H
