#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
PDB_ID="2PC0"
SAVE_DIR="../files/structures"

cp $SAVE_DIR/$PDB_ID.pdb $SAVE_DIR/cleaned.pdb

simlify-pdb-filter $SAVE_DIR/cleaned.pdb --output $SAVE_DIR/cleaned.pdb

# Remove water molecules.
sed -i "/HOH/d" "$SAVE_DIR/cleaned.pdb"
sed -i "/ MG /d" "$SAVE_DIR/cleaned.pdb"
sed -i "/ PGR /d" "$SAVE_DIR/cleaned.pdb"

# Remove duplicate atom lines.
sed -i "/BARG/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/AARG/ ARG/g" "$SAVE_DIR/cleaned.pdb"

sed -i "/BLYS/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/ALYS/ LYS/g" "$SAVE_DIR/cleaned.pdb"

sed -i "/BGLU/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/AGLU/ GLU/g" "$SAVE_DIR/cleaned.pdb"

sed -i "/BVAL/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/AVAL/ VAL/g" "$SAVE_DIR/cleaned.pdb"

sed -i "/BMET/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/AMET/ MET/g" "$SAVE_DIR/cleaned.pdb"

sed -i "/BILE/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/AILE/ ILE/g" "$SAVE_DIR/cleaned.pdb"

sed -i "/BTHR/d" "$SAVE_DIR/cleaned.pdb"
sed -i "/CTHR/d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/ATHR/ THR/g" "$SAVE_DIR/cleaned.pdb"

# Turn models into chains
sed -i "/MODEL   /d" "$SAVE_DIR/cleaned.pdb"
sed -i "s/ENDMDL/TER/g" "$SAVE_DIR/cleaned.pdb"
simlify-pdb-reorder $SAVE_DIR/cleaned.pdb --output $SAVE_DIR/cleaned.pdb


)
