#!/usr/bin/env bash

source /sphenix/u/shuhang98/setup.sh

# Check if exactly one argument is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <argument>"
  exit 1
fi

# Store the argument
arg="$1"

# Execute the ROOT command with the argument
root -b -q "plot_sideband_selection.C(\"${arg}\")" &
root -b -q "plot_sideband_sim_selection.C(\"${arg}\")" &
root -b -q "plot_purity_selection.C(\"${arg}\")" &
root -b -q "plot_purity_sim_selection.C(\"${arg}\")" &
root -b -q "plot_isoET.C(\"${arg}\")" &
root -b -q "plot_efficiency.C(\"${arg}\")" &
root -b -q "plot_final_selection.C(\"${arg}\")" &
root -b -q "CONF_plots.C(\"${arg}\", 0, 2)" &
root -b -q "CONF_plots.C(\"${arg}\", 3, 5)" &
root -b -q "CONF_plots.C(\"${arg}\", 6, 8)" &
root -b -q "CONF_plots.C(\"${arg}\", 9, 9)" &

wait


#gs "figures/purity_sim_${arg}.pdf" 
#gs "figures/purity_${arg}.pdf" 
#gs "figures/et_sbs_${arg}.pdf" 
#gs "figures/et_sbs_ratio_${arg}.pdf"  

