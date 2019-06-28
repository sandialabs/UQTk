#!/bin/bash

# python
/bin/rm -rf *.pyc

# pdf files
/bin/rm -rf forUQ*.pdf

# log files
/bin/rm -rf log*.dat

# intermediate files
/bin/rm -rf jointsens.dat mainsens.dat totsens.dat varfrac.dat
/bin/rm -rf PCEspls_100k.dat output_val.dat output_val_pc.dat
/bin/rm -rf solution.dat sp_mindex.1.dat PCcoeff_quad.dat mipc.dat
/bin/rm -rf surf_rxn.in.parsed.xml surf_rxn.in.xml
/bin/rm -rf input.dat output.dat input_fcn.dat output_fcn.dat input_val.dat
/bin/rm -rf xy*.dat

