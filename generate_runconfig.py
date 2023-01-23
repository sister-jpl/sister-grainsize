#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""

import sys
import json

def main():
    '''
        This function takes as input the path to an inputs.json file and exports a run config json
        containing the arguments needed to run the L2A snow grain size PGE.

    '''

    inputs_json  = sys.argv[1]

    with open(inputs_json, 'r') as in_file:
        inputs =json.load(in_file)

    run_config = {"inputs":{}}

    for file_dict in inputs["file"]:
        for key,value in file_dict.items():
            run_config["inputs"][key] = value

    run_config["inputs"].update(inputs["config"])

    config_file = 'runconfig.json'

    with open(config_file, 'w') as outfile:
        json.dump(run_config,outfile,indent=3)


if __name__=='__main__':
    main()
