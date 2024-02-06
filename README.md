# SISTER Snow grain size PGE Documentation

## Description
The L2B snow grain size PGE takes as input surface reflectance and calculates snow grain size using the method of Nolin and Dozier (2000). Snow grain size is modeled as a function of scaled band area centered at the 1030 nm ice absorption feature.

![](./examples/prisma_snow_spectrum.png)

![](./examples/prisma_grainsize_example1.png)

### References 
- Nolin, A. W., & Dozier, J. (2000).
A hyperspectral method for remotely sensing the grain size of snow.
Remote sensing of Environment, 74(2), 207-216.
[doi.org/10.1016/S0034-4257(00)00111-5](https://doi.org/10.1016/S0034-4257(00)00111-5)

## PGE Arguments

The L2B snow grain size PGE takes the following argument(s):


| Argument            | Description                          | Default |
|---------------------|--------------------------------------|---------|
| reflectance_dataset | L2A reflectance dataset              | -       |
| frcover_dataset     | L2B fractional cover dataset         | -       |
| snow_cover          | Snow fractional cover threshold      | 0.9     |
| experimental        | Designates outputs as "experimental" | 'True'  |

## Outputs

The outputs of the L2B snow grainsize PGE use the following naming convention:

    (EXPERIMENTAL-)SISTER_<SENSOR>_L2B_GRAINSIZE_<YYYYMMDDTHHMMSS>_<CRID>

Note that the "EXPERIMENTAL-" prefix is optional and is only added when the "experimental" flag is set to True.

The following data products are produced:

| Product description                      | Units   | Example filename                                                  |
|------------------------------------------|---------|-------------------------------------------------------------------|
| Snow grain size COGeotiff                | -       | SISTER\_AVCL\_L2B\_GRAINSIZE\_20110513T175417\_001.tif            |
| 1. Snow grain size                       | microns | -                                                                 |
| 2. Quality assurance mask                | -       | -                                                                 |
| Snow grainsize metadata (STAC formatted) | -       | SISTER\_AVCL\_L2B\_GRAINSIZE\_20110513T175417\_001.json           |
| Snow grainsize quicklook                 | -       | SISTER\_AVCL\_L2B\_GRAINSIZE\_20110513T175417\_001.png            |
| PGE runconfig                            | -       | SISTER\_AVCL\_L2B\_GRAINSIZE\_20110513T175417\_001.runconfig.json |
| PGE log                                  | -       | SISTER\_AVCL\_L2B\_GRAINSIZE\_20110513T175417\_001.log            |

Metadata files are [STAC formatted](https://stacspec.org/en) and compatible with tools in the [STAC ecosystem](https://stacindex.org/ecosystem).

## Executing the Algorithm

This algorithm requires [Anaconda Python](https://www.anaconda.com/download)

To install and run the code, first clone the repository and execute the install script:

    git clone https://github.com/sister-jpl/sister-grainsize.git
    cd sister-grainsize
    ./install.sh
    cd ..

Then, create a working directory and enter it:

    mkdir WORK_DIR
    cd WORK_DIR

Copy input files to the work directory. For each "dataset" input, create a folder with the dataset name, then download 
the data file(s) and STAC JSON file into the folder.  For example, the reflectance dataset input would look like this:

    WORK_DIR/SISTER_AVCL_L2A_CORFL_20110513T175417_000/SISTER_AVCL_L2A_CORFL_20110513T175417_000.bin
    WORK_DIR/SISTER_AVCL_L2A_CORFL_20110513T175417_000/SISTER_AVCL_L2A_CORFL_20110513T175417_000.hdr
    WORK_DIR/SISTER_AVCL_L2A_CORFL_20110513T175417_000/SISTER_AVCL_L2A_CORFL_20110513T175417_000.json

Finally, run the code 

    ../sister-grainsize/pge_run.sh --reflectance_dataset SISTER_AVCL_L2A_CORFL_20110513T175417_000 --frcov_dataset SISTER_AVCL_L2B_FRCOV_20110513T175417_000