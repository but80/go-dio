#!/bin/bash

cd "$(dirname "$0")"

# Download original World source
if ! [ -d .original ]; then
  git clone https://github.com/mmorise/World.git .original || exit $?
fi
cd .original || exit $?
git reset --hard ad770cc00c45ef1f00a8c1ed2c1d2dce8d26bada || exit $?
cd .. || exit $?

# Convert wav file (for test) to C++/Go source code
sox testdata.wav testdata.dat || exit $?
echo 'double test_data[] = {' > testdata.h || exit $?
cat testdata.dat | tr -d '\r' | grep -v '^;' | sed 's/^ *[^ ]* */\t/' | sed 's/ *$/,/' >> testdata.h || exit $?
echo '};' >> testdata.h || exit $?

# Compile
gcc -I. -I./.original/src -o make-testdata *.cpp ./.original/src/*.cpp -lm -lstdc++ || exit $?
