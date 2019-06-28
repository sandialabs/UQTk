#!/bin/bash

if [ ! -f "delta.dat" ]
then
  touch delta.dat
  echo 2 > delta.dat
fi

./modelLL.x
