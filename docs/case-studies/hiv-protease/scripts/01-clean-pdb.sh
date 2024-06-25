#!/usr/bin/env bash

(
# Ensure we are in the directory of this script.
cd "$(dirname "$0")"

# Set environmental variables
export SIMLIFY_LOG=True
export SIMLIFY_STDOUT=False
export SIMLIFY_LOG_LEVEL=0
export SIMLIFY_LOG_FILE_PATH="01-clean-pdb.log"

rm 01-clean-pdb.log

# Set environmental variables
PDB_ID="2PC0"
SAVE_DIR="../files/structures"

INPUT_PATH=$SAVE_DIR/$PDB_ID.pdb
OUTPUT_PATH=$SAVE_DIR/$PDB_ID-cleaned.pdb

cp $INPUT_PATH $OUTPUT_PATH

simlify-pdb-filter $OUTPUT_PATH --output $OUTPUT_PATH

# Remove water molecules and other ligands.
sed -i "/HOH/d" "$OUTPUT_PATH"
sed -i "/ MG /d" "$OUTPUT_PATH"
sed -i "/ PGR /d" "$OUTPUT_PATH"

# Remove duplicate atom lines.
sed -i "/BARG/d" "$OUTPUT_PATH"
sed -i "s/AARG/ ARG/g" "$OUTPUT_PATH"

sed -i "/BLYS/d" "$OUTPUT_PATH"
sed -i "s/ALYS/ LYS/g" "$OUTPUT_PATH"

sed -i "/BGLU/d" "$OUTPUT_PATH"
sed -i "s/AGLU/ GLU/g" "$OUTPUT_PATH"

sed -i "/BVAL/d" "$OUTPUT_PATH"
sed -i "s/AVAL/ VAL/g" "$OUTPUT_PATH"

sed -i "/BMET/d" "$OUTPUT_PATH"
sed -i "s/AMET/ MET/g" "$OUTPUT_PATH"

sed -i "/BILE/d" "$OUTPUT_PATH"
sed -i "s/AILE/ ILE/g" "$OUTPUT_PATH"

sed -i "/BTHR/d" "$OUTPUT_PATH"
sed -i "/CTHR/d" "$OUTPUT_PATH"
sed -i "s/ATHR/ THR/g" "$OUTPUT_PATH"

# Turn models into chains
sed -i "/MODEL   /d" "$OUTPUT_PATH"
sed -i "s/ENDMDL/TER/g" "$OUTPUT_PATH"

simlify-pdb-renumber $OUTPUT_PATH --output $OUTPUT_PATH

)
