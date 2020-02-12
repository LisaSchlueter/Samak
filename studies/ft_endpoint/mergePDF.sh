#!/bin/bash
# this is a executable to merge pdfs (created in the FitRunList function in MultiRunAnalysis)
cd plots/tmp
gs -sDEVICE=pdfwrite -sOutputFile="RunSummary.pdf" -dNOPAUSE -dEPSCrop -c "<</Orientation 0>> setpagedevice" -f KATRIN_FT_AllPixels_Samak_*.pdf -c quit
cd ../..
