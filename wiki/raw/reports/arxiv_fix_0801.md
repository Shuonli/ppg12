# arXiv ID Fix: 0801.4020 -> 2005.14270

## Problem

The PPG12 reading list (`wiki/proposed_reading_list.md`) entry:

| Priority | arXiv | Description | Collaboration | System | Relevance |
|----------|-------|-------------|---------------|--------|-----------|
| 1 | 0801.4020 | Azimuthal gamma-hadron dAu+AuAu 200 GeV | PHENIX | dAu+AuAu 200 | Early photon-hadron azimuthal study |

**arXiv:0801.4020 is the WRONG paper.** It is actually:

- **Title:** "Suppression pattern of neutral pions at high transverse momentum in Au+Au collisions at sqrt(s_NN) = 200 GeV and constraints on medium transport coefficients"
- **Journal:** Phys. Rev. Lett. 101, 232301 (2008)
- **Content:** Single-hadron pi0 R_AA measurement, NOT gamma-hadron correlations

## Correct Paper

The intended paper is **arXiv:2005.14270**:

- **Title:** "Measurement of jet-medium interactions via direct photon-hadron correlations in Au+Au and d+Au collisions at sqrt(s_NN) = 200 GeV"
- **Collaboration:** PHENIX
- **arXiv:** 2005.14270 [nucl-ex], 28 May 2020
- **Journal:** Phys. Rev. C 102, 054910 (2020)
- **DOI:** 10.1103/PhysRevC.102.054910
- **INSPIRE texkey:** PHENIX:2020alr
- **Abstract:** Presents direct photon-hadron correlations in 200 GeV/A Au+Au, d+Au, and p+p collisions, for direct photon pT from 5-12 GeV/c, collected 2006-2011. Observes no significant modification in d+Au (cold nuclear matter effects small). Suppression of high-zT hadrons in Au+Au; excess at low-zT and large angles attributed to medium response.

## Other Candidate Papers Considered

During the search, the following related PHENIX direct photon-hadron papers were identified:

1. **arXiv:0903.3399** -- "Photon-Hadron Jet Correlations in p+p and Au+Au Collisions at sqrt(s) = 200 GeV"
   - Phys. Rev. C 80, 024908 (2009)
   - First observation of away-side suppression in direct photon+jet channel in Au+Au vs p+p
   - Does NOT include d+Au data

2. **arXiv:1212.3323** -- "Medium modification of jet fragmentation in Au+Au collisions at sqrt(s_NN) = 200 GeV measured in direct photon-hadron correlations"
   - Phys. Rev. Lett. 111, 032301 (2013)
   - Jet fragmentation function via zT in Au+Au vs p+p
   - Does NOT include d+Au data

3. **arXiv:1006.1347** -- "High pT direct photon and pi0 triggered azimuthal jet correlations and measurement of kT for isolated direct photons in p+p collisions at sqrt(s) = 200 GeV"
   - Phys. Rev. D 82, 072001 (2010)
   - p+p only, kT measurement

Only **arXiv:2005.14270** includes d+Au data alongside Au+Au, matching the reading list description.

## Note on Publication Date

The reading list description ("early photon-hadron azimuthal study") suggests an older paper was expected. However, 2005.14270 (published 2020) is the ONLY PHENIX journal article presenting direct photon-hadron correlations in both d+Au and Au+Au at 200 GeV. The earlier PHENIX papers (0903.3399 and 1212.3323) cover only p+p and Au+Au. The d+Au measurement was first published in 2020 using data collected 2006-2011.

## Action Required

Replace `0801.4020` with `2005.14270` in the reading list. Update the description to:

| Priority | arXiv | Description | Collaboration | System | Relevance |
|----------|-------|-------------|---------------|--------|-----------|
| 1 | 2005.14270 | Direct photon-hadron correlations in AuAu+dAu+pp 200 GeV (PRC 102, 2020) | PHENIX | pp+dAu+AuAu 200 | Photon-hadron jet-medium interactions; d+Au baseline for cold nuclear matter |

## Search Method

Searched INSPIRE-HEP API with queries:
- `find cn phenix and t photon hadron and d > 2008 and d < 2015`
- `collaboration phenix and title "photon-hadron"`
- `find eprint 0801.4020` (confirmed wrong paper)
- `find eprint 2005.14270` (confirmed correct paper)
