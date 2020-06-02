# plot_spectra
The script enables to plot IR spectra with arbitrary x-axis interruptions and more.

## Prerequisites
- Numpy
- Scipy
- matplotlib

## Installation

Copy both files anywere you like or make an alias to it.

## Usage

plot.py --help

Most options are self explanatory, but some require knowledge about the argument parsing. Go to the example section for more details

## Examples

For the main figure of the paper *insert here when published*, I used the following options:

```
  ./plot.py spectra_after_tunnelling.csv computed_spectra/cc_dihydroxycarbene_AE-CCSDT_aug_cc-pCVTZ_vpt2.dpt \
  computed_spectra/ct_dihydroxycarbene_AE-CCSDT_aug_cc-pCVTZ_vpt2.dpt Figure1.pdf --plotLimitsX '3800-3200 1600-400' \
  --colors black blue darkgoldenrod \
  --labels 'Diff. spectrum after Tunneling' '\textbf{1cc} VPT2-AE-CCSD(T)/cc-pCVTZ' \
  '\textbf{1ct} VPT2-AE-CCSD(T)/cc-pCVTZ' \
  --invert no no yes --yshift -0.0 -10000.0 4000.0 --yaxisIndex 0 1 1 \
  --plotLimitsY '-0.2-0.2 -10500-4500.0' \
  --assignments 'None None None 1260.0 None None 730.1 None None 3592.4 None 1386.17 1295.6 1148.4 None None 670.0 None' \
  '3286.4 3238.3 1325.6 1275.0 1104.7 1096.7 705.7 658.8 532.6' \
  '3638.608 3362.839 1420.397 1293.912 1130.269 1078.932 740.880 625.794 622.062'
```

## Known issues
- The assignment does not work reliable
- Interactive option not fully developted and is just a preview
