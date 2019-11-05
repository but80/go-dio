#!/bin/bash

cd "$(dirname "$0")"

# Download original World source
if ! [ -d .original ]; then
  git clone https://github.com/mmorise/World.git .original || exit $?
fi
cd .original || exit $?
git reset --hard ad770cc00c45ef1f00a8c1ed2c1d2dce8d26bada || exit $?
cd .. || exit $?

# Compile
gcc -I. -I./.original/src -o make-testdata *.cpp ./.original/src/*.cpp -lm -lstdc++ -lsndfile || exit $?
