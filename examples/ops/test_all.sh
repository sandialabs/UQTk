#!/bin/sh 

# Test elementary operations with PCEs
# See also UQTk manual, Chapter on examples for more details.

# Run elementary operations program
./Ops.x

# Plot pdfs of a and log(a)
./plot_pdf.py samples.a.dat
./plot_pdf.py samples.loga.dat

# Compare the Taylor and integration approach for log(a)
./LogComp.x
./plot_logs.py
