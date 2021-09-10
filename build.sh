mkdir -p build
pushd build >/dev/null 2>&1
  ../config/config-gcc-Python.sh
  make -j 8
  ctest
  make install
popd >/dev/null 2>&1
