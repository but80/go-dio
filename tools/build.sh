#!/bin/bash

cd "$(dirname "$0")"

# Download original World source
if ! [ -d world-cpp ]; then
  git clone https://github.com/mmorise/World.git world-cpp || exit $?
fi
cd world-cpp || exit $?
git reset --hard ad770cc00c45ef1f00a8c1ed2c1d2dce8d26bada || exit $?
cd .. || exit $?

# Compile
gcc -I. -I./world-cpp/src -o make-testdata *.cpp ./world-cpp/src/*.cpp -lm -lstdc++ -lsndfile || exit $?
