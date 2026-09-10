#!/usr/bin/env bash
# Fetch the pre-arXiv direct/isolated-photon datasets that complete the Bock HP2018
# Fig.1 compilation. These live ONLY on HEPData, whose download API is Cloudflare-
# blocked from the analysis cluster. RUN THIS ON A MACHINE WHERE HEPData IS REACHABLE
# (your laptop browser network), then copy plotting/xtscaling/data/hepdata/ back here.
#
# If curl is also Cloudflare-challenged on your laptop, just open each "record" URL
# below in a browser and click "Download All -> YAML", saving into data/hepdata/.
#
# After the files are in place:  ask me to ingest them and I will convert + replot.

set -u
OUT="$(cd "$(dirname "$0")" && pwd)/data/hepdata"
mkdir -p "$OUT"

# name                inspire-id   sqrt(s)[GeV]  eta-range          iso  (from PHENIX PRD86 072008 appendix)
records=(
  "CDF_1800           ins375582    1800   |eta|<0.9     iso"
  "D0_1800            ins417044    1800   |eta|<0.9     iso"
  "UA1_630_546        ins261356     630   |eta|<0.8     iso"
  "UA2_630            ins336186     630   |eta|<0.76    iso"
  "R110_CMOR_63       ins279988      63   |eta|<0.8     iso"
  "R806_63            ins176956      63   |eta|<0.2     iso"
  "R108_CCOR_62       ins153483    62.4   |eta|<0.45    iso"
  "E704_19p4          ins399938    19.4   |xF|<0.15     iso"
  "NA24_23p8          ins236248    23.8   -0.65<eta<0.52 noiso"
  "WA70_23            ins252633    23.0   |xF|<0.05     noiso"
  "UA6_24p3           ins476703    24.3   -0.1<eta<0.9  noiso"
)
# NOTE: R807 (AFS, Akesson Sov.J.Nucl.Phys.51) has no HEPData record -> not fetchable.

for r in "${records[@]}"; do
  name=$(echo "$r" | awk '{print $1}')
  ins=$(echo  "$r" | awk '{print $2}')
  url="https://www.hepdata.net/download/submission/${ins}/1/yaml"
  echo "[$name] $ins"
  echo "   record : https://www.hepdata.net/record/${ins}"
  echo "   yaml   : $url"
  curl -fsSL "$url" -o "$OUT/${name}_${ins}.tar.gz" \
       -A "Mozilla/5.0 (X11; Linux x86_64) Gecko/20100101 Firefox/128.0" \
    && echo "   -> saved $OUT/${name}_${ins}.tar.gz" \
    || echo "   -> FAILED (Cloudflare?); download the 'record' URL in a browser instead"
done
echo
echo "Done. If files landed in $OUT, copy that folder to the cluster and ask me to ingest."
