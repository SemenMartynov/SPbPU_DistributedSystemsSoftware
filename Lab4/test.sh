#!/bin/bash

# Check for arguments
if [ $# -lt 1 ]; then
  echo "please use: $0 <file_name.cpp>"
  exit 1
fi
# Get the file name from the first argument
filename="$1"

# Create an array of parameters
params=("-O0" "-Os" "-O1" "-O2" "-O3" "-O2 -march=native" "-O3 -march=native" "-O2 -march=native -funroll-loops" "-O3 -march=native -funroll-loops")

# Iterate through the array and execute the process for each parameter
for param in "${params[@]}"; do
  echo "executing 'g++ $param $filename -o test'"
  g++ $param $filename -o test
  time ./test
  du -h ./test
  rm ./test
done
