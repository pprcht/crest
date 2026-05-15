#!/bin/bash
# Standalone CREGEN sorting of an ensemble/trajectory file.
# Removes duplicate structures (by RMSD and energy) and sorts by energy.
# Output: crest_conformers.xyz (unique), crest_rotamers.xyz (all), xtb.trj.sorted
# Note: -cregen standalone mode has no TOML runtype; use the CLI form below.

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run ---
crest -cregen xtb.trj -ewin 100.0

# (no TOML runtype equivalent for standalone CREGEN; see input.toml for
#  the [cregen] settings used during a full iMTD-GC run)
