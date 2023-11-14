#!/bin/bash

source activate sister

pge_dir=$(cd "$(dirname "$0")" ; pwd -P)

echo "Creating runconfig"
python ${pge_dir}/generate_runconfig.py "${@:1}"

echo "Running L2B snow grain size PGE"
python ${pge_dir}/grainsize.py runconfig.json
