#!/bin/bash



# Test function with 3 inputs and 4 outputs
#awk '{print $3*log($1**2)+$1+($2**2)*(1-exp(-$3**2)), $1+$3**2, $2+$3, $1*$3}' $1 > $2

# Another test function
#./model.py $1 $2

# Sensitivity
rm -rf ms ts
./transpose_file.x $1 > pccfs.dat

# HARDWIRED!!
#gen_mi -x TO -p 5 -q 5; cp mindex.dat mi_prop.txt

NSAM=`echo | awk 'END{print NR}' $1`

for (( i=1; i<=$NSAM ; i++ )); do

	awk '{print $i}' i=$i pccfs.dat > pcf_this.dat
	${UQTK_INS}/bin/pce_sens -x LU -m mi.dat -f pcf_this.dat > pcsens.log
	./transpose_file.x mainsens.dat >> ms
	./transpose_file.x totsens.dat >> ts
done

paste ms ts > $2
