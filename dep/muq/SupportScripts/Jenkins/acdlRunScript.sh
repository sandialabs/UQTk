#!/bin/bash
shopt -s nocasematch

#######################################
##### Was MUQ built with Python?
#######################################

dir=`pwd`
bin_dir=$(echo "$dir" | tr '[:lower:]' '[:upper:]')

if echo $bin_dir| grep -q "PYTHON"; then
  with_python=1
else
  with_python=0
fi

if echo $bin_dir| grep -q "NLOPT"; then
  with_nlopt=1
else
  with_nlopt=0
fi

#######################################
##### Update the library path if on OSX
#######################################
if [[ `hostname` == "macys.mit.edu" ]]; then
  export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/jenkins/util/gtest/lib/:$WORKSPACE/install/lib:$WORKSPACE/install/muq_external/lib
  PYTEST=py.test
fi

#######################################
##### Update the library path if on Linux
#######################################
if [[ `hostname`=="reynolds.mit.edu" ]]; then
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WORKSPACE/install/lib:$WORKSPACE/install/muq_external/lib
  PYTEST=/usr/local/bin/py.test
fi

#######################################
##### Run the tests
#######################################
build/RunAllTests --gtest_output=xml:results/tests/TestResults.xml
if [ $with_python -eq 1 ]; then
  export PYTHONPATH=$PYTHONPATH:$WORKSPACE/install/lib

  if [ $with_nlopt -eq 1 ]; then
    $PYTEST -v --junitxml results/tests/PythonTestResults.xml modules/RunPythonTests.py
  else
    $PYTEST -v -k "not Nlopt and not PolynomialApproximator" --junitxml results/tests/PythonTestResults.xml modules/RunPythonTests.py
  fi
fi

#######################################
##### Run the examples
#######################################
for f in $(find $WORKSPACE/examples/ -name build -prune -o -type d -name example-* -print); do
  cd $f
  # run any executables that were created in the build directory
  for e in $(find build/ -maxdepth 1 -type f -perm +0111); do
    $e
  done
done

cd $dir

exit 0
