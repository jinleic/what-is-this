#!/bin/bash
# benchmark.sh — rebuild + equivalence + timing for the beam8 C++ port.
# Run from THIS directory (src/beam_cpp).
set -e
PY=${PY:-/Users/jinleic/jinleic-workspace/physics/.venv/bin/python}
SHOTS=${SHOTS:-1000}
SEED=${SEED:-20260830}

echo "== build =="
clang++ -O3 -std=c++17 -Wall -Wextra beam8.cpp -o beam8_cpp

echo
echo "== equivalence + timing: p=3e-3 (companion reference: frozen beam8 2362 ms/shot over 2000 shots) =="
PYTHONPATH=../../src "$PY" equiv.py --p 3e-3 --shots "$SHOTS" --seed "$SEED" --cpp-repeats 3

echo
echo "== equivalence + timing: p=1e-3 (pinned Gate B point) =="
PYTHONPATH=../../src "$PY" equiv.py --p 1e-3 --shots "$SHOTS" --seed "$SEED" --cpp-repeats 3

echo
echo "== reports =="
ls -l scratch/equiv_p3e-3_n${SHOTS}_s${SEED}.json scratch/equiv_p1e-3_n${SHOTS}_s${SEED}.json
