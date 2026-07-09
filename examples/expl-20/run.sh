#!/bin/bash
# Geometry optimization of n-pentane using an external ORCA subprocess as the
# energy+gradient backend (r2scan-3c as defined in the ORCA.in template).
#
# CREST runs its ANCOPT optimizer and calls ORCA for every engrad evaluation.
# The ORCA calculator is configured through the TOML input only (there is no
# dedicated CLI flag for it).
#
# NOTE: the ORCA executable path in input.toml is a placeholder
# (/PATH/TO/ORCA/orca) and must be set to your local ORCA installation for
# this example to actually run.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- TOML run ---
crest input.toml
