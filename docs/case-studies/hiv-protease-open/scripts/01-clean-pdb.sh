#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export PDB_ID="2PC0"
SAVE_DIR="../files/structures"

# Cleanup files from previous run
rm -rf $SAVE_DIR
mkdir -p $SAVE_DIR

# Get PDB file

)
