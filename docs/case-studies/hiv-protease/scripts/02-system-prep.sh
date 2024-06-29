#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export SIMLIFY_LOG=True
export SIMLIFY_STDOUT=False
export SIMLIFY_LOG_LEVEL=10
export SIMLIFY_LOG_FILE_PATH="02-amber-prep.log"

rm $SIMLIFY_LOG_FILE_PATH

PDB_PATH="../files/structures/2PC0-cleaned.pdb"

CONTEXTS_DIR="../files/contexts"
YAML_PATH="$CONTEXTS_DIR/base.yml"

SAVE_DIR="../simulations/02-prep"
TOPO_PATH="$SAVE_DIR/mol.prmtop"
COORD_PATH="$SAVE_DIR/mol.inpcrd"
OUTPUT_PDB_PATH="$SAVE_DIR/mol.pdb"



# Cleanup files from previous run
rm -rf "$SAVE_DIR"
mkdir -p "$SAVE_DIR"
rm -f $SIMLIFY_LOG_FILE_PATH

simlify-tleap $PDB_PATH $TOPO_PATH $COORD_PATH --yaml $YAML_PATH --work_dir "$(dirname "$0")"
simlify-pdb $OUTPUT_PDB_PATH --files $TOPO_PATH $COORD_PATH

export SIMLIFY_LOG=False
)
