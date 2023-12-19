#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""

import argparse
import sys
import json

def main():
    '''
        This function takes as input the path to an inputs.json file and exports a run config json
        containing the arguments needed to run the L2A snow grain size PGE.

    '''

    parser = argparse.ArgumentParser(description='Parse inputs to create runconfig.json')
    parser.add_argument('--reflectance_dataset', help='Path to reflectance dataset')
    parser.add_argument('--frcov_dataset', help='Path to uncertainty dataset')
    parser.add_argument('--snow_cover', help='Minimum snow cover threshold')
    parser.add_argument('--crid', help='CRID value')
    parser.add_argument('--experimental', help='If true then designates data as experiemntal')
    args = parser.parse_args()

    run_config = {
        "inputs": {
            "reflectance_dataset": args.reflectance_dataset,
            "frcov_dataset": args.frcov_dataset,
            "snow_cover": float(args.snow_cover),
            "crid": args.crid
        }
    }
    run_config["inputs"]["experimental"] = True if args.experimental.lower() == "true" else False

    # Add metadata to runconfig
    corfl_basename = os.path.basename(run_config["inputs"]["reflectance_dataset"])

    met_json_path = os.path.join("input", corfl_basename, f"{corfl_basename}.met.json")
    with open(met_json_path, "r") as f:
        metadata = json.load(f)
    run_config["metadata"] = metadata

    config_file = 'runconfig.json'

    with open(config_file, 'w') as outfile:
        json.dump(run_config,outfile,indent=3)


if __name__=='__main__':
    main()
