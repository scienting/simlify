#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
PDB_ID="2PC0"
SAVE_DIR="../files/structures"

cp $SAVE_DIR/$PDB_ID.pdb $SAVE_DIR/cleaned.pdb

simlify-pdb-filter $SAVE_DIR/cleaned.pdb --output $SAVE_DIR/cleaned.pdb
simlify-pdb-reorder $SAVE_DIR/cleaned.pdb --output $SAVE_DIR/cleaned.pdb

)
