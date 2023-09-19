#!/bin/bash

# Check for arguments
if [ $# -lt 1 ]; then
  echo "please use: $0 <file_name.cpp>"
  exit 1
fi
# Get the file name from the first argument
filename="$1"

# Create an array of parameters
params=("-O3 -march=native -funroll-loops -fipa-pta" "-O3 -march=native -funroll-loops -flto=auto" "-O3 -march=native -funroll-loops -fipa-pta -flto=auto")

# Iterate through the array and execute the process for each parameter
for param in "${params[@]}"; do
  echo "executing 'g++ $param $filename -o test'"
  g++ $param $filename -o test
  time ./test
  du -h ./test
  rm ./test
done
