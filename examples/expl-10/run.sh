#!/bin/bash
# Constrained iMTD-GC conformer search of 1-propanol:
# the C-C-C-O backbone (atoms 1-4) is frozen; only the OH dihedral is sampled.
# Expected output: 2 conformers (different OH rotamers) in crest_conformers.xyz.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
# Step 1: generate a constraint template (.xcontrol.sample):
crest struc.xyz -constrain 1-4
# Step 2: run the constrained conformer search:
crest struc.xyz -gfnff -cinp .xcontrol.sample -ewin 2.0

# --- TOML run (constraints defined directly in the input file) ---
# crest input.toml
