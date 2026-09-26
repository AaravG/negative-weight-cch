#!/bin/sh
# Final benchmark runs for the paper (medians over 5 repetitions).
set -e
cd "$(dirname "$0")"
./cch.exe data/NY  --queries 1000 --check 1000 --shifted-cch --repeat 5 --parallel 16 > ../results/cpp_NY.txt 2>&1
./cch.exe data/BAY --queries 1000 --check 1000 --shifted-cch --repeat 5 --parallel 16 > ../results/cpp_BAY.txt 2>&1
./cch.exe data/USA --order data/USA/order.i32 --queries 1000 --check 100 --shifted-cch --repeat 5 --parallel 16 > ../results/cpp_USA.txt 2>&1
echo done
