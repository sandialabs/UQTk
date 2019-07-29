#!/bin/bash
# =====================================================================================
#                      The UQ Toolkit (UQTk) version @UQTKVERSION@
#                     Copyright (@UQTKYEAR@) Sandia Corporation
#                      http://www.sandia.gov/UQToolkit/
#
#     Copyright (@UQTKYEAR@) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
#     with Sandia Corporation, the U.S. Government retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public License
#     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
#
#     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
# =====================================================================================
if [ "$#" -ne 4 ]; then
  echo "usage: $0 nspl npar modeval modpar                       "
  echo "        nspl    : no. of samples                         "
  echo "        npar    : no. of parameters                      "
  echo "        modeval : filename for model evaluations         "
  echo "        modtrim : filename for trimmed model evaluations "
  echo "  Note: the id's of unwanted model evaluations should be "
  echo "        in file \"spllst.dat\""
  exit
fi

nspl=$1
npar=$2
modeval=$3
modtrim=$4

# compute no. of model outputs
neval=`echo "$nspl $npar" | awk '{print $1*($2+2)}'`

# get id's for bad values
spllst=( `cat "spllst.dat"`)
extlst=("${spllst[@]}")
#echo ${spllst[@]}

# create extended list
for i in ${spllst[@]}; do
  j=$[$i-$nspl]
  while [ $j -gt 0 ]; do
    extlst+=("$j")
    j=$[$j-$nspl]
    #echo ${extlst[@]}
  done
  j=$[$i+$nspl]
  while [ $j -le $neval ]; do
    extlst+=("$j")
    j=$[$j+$nspl]
    #echo ${extlst[@]}
  done
done


# sort list of id's to be eliminated and eliminate duplicates id's
extlst=(`sort -n <(printf "%s\n" "${extlst[@]}")`)
extlst=(`echo ${extlst[@]} | awk 'BEGIN{RS=ORS=" "}{ if (a[$0] == 0){ a[$0] += 1; print $0}}'`)

# Eliminate unwanted evaluations
/bin/cp $modeval $modtrim
for ((i=${#extlst[@]} - 1;i >= 0;i--)); do
    #echo "\"${extlst[i]}\""
    awk -v j=${extlst[i]} 'NR!=j {print $0}' $modtrim > tmp.tmp
    /bin/mv tmp.tmp $modtrim
done

# Compute number of effective samples left after trimming
nspl=`wc -l $modtrim | awk -v np=$npar '{print $1/(np+2)}'`
echo "Number of samples left: $nspl"
