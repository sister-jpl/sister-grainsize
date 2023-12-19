# -*- coding: utf-8 -*-
"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus
"""

import datetime as dt
import glob
import json
import os
import sys
import shutil
import hytools as ht
import numpy as np
import matplotlib.pyplot as plt
from osgeo import gdal
from scipy.interpolate import interp1d
from PIL import Image
import pystac

def main():

    ''' Estimate snow grain size from rfl data using method of Nolin and Dozier (2000).

    Nolin, A. W., & Dozier, J. (2000).
    A hyperspectral method for remotely sensing the grain size of snow.
    Remote sensing of Environment, 74(2), 207-216.
    https://doi.org/10.1016/S0034-4257(00)00111-5

    '''

    pge_path = os.path.dirname(os.path.realpath(__file__))

    run_config_json = sys.argv[1]

    with open(run_config_json, 'r') as in_file:
        run_config =json.load(in_file)

    experimental = run_config['inputs']['experimental']
    if experimental:
        disclaimer = "(DISCLAIMER: THIS DATA IS EXPERIMENTAL AND NOT INTENDED FOR SCIENTIFIC USE) "
    else:
        disclaimer = ""

    os.mkdir('output')
    os.mkdir('temp')

    crid = run_config["inputs"]["crid"]

    rfl_base_name = os.path.basename(run_config['inputs']['reflectance_dataset'])
    sister,sensor,level,product,datetime,in_crid = rfl_base_name.split('_')
    rfl_file = f'input/{rfl_base_name}/{rfl_base_name}.bin'
    rfl_met = rfl_file.replace('.bin','.met.json')

    fc_base_name = os.path.basename(run_config['inputs']['frcov_dataset'])
    fc_file = f'input/{fc_base_name}/{fc_base_name}.tif'

    grain_met = f'output/SISTER_{sensor}_L2B_SNOWGRAIN_{datetime}_{crid}.met.json'

    data = f'{pge_path}/data/grainsize_v_bandarea_nolin_dozier_2000_interp.csv'
    interp_data = np.loadtxt(data,
                             delimiter = ',').T
    interpolator = interp1d(interp_data[0],interp_data[1],
                            kind = 'cubic',fill_value ='extrapolate')

    # Set fractional cover mask
    fc_obj = gdal.Open(fc_file)
    snow_mask = fc_obj.GetRasterBand(4).ReadAsArray() >= run_config['inputs']['snow_cover']

    rfl = ht.HyTools()
    rfl.read_file(rfl_file,'envi')

    # Get continuum start and end bands
    wave1, wave2  = 950,1080
    wave1_ind = rfl.wave_to_band(wave1)
    wave2_ind = rfl.wave_to_band(wave2)

    grain_size = np.full((rfl.lines,rfl.columns),-1)
    iterator =rfl.iterate(by = 'chunk',chunk_size = (200,200))
    i = 0

    while not iterator.complete:
        chunk = iterator.read_next()[:,:,wave1_ind:wave2_ind+1]

        # Calculate continuum removed spectrum
        slope = (chunk[:,:,-1]-chunk[:,:,0])/(wave2-wave1)
        intercept = chunk[:,:,-1]- wave2*slope
        continuum = np.einsum('i,mj->mji', rfl.wavelengths[wave1_ind:wave2_ind+1],slope) + np.expand_dims(intercept,axis=2)

        depth  = (continuum-chunk)/continuum
        band_area =np.trapz(depth,axis=2)

        grain_chunk = interpolator(band_area)

        grain_size[iterator.current_line:iterator.current_line+chunk.shape[0],
              iterator.current_column:iterator.current_column+chunk.shape[1]] = grain_chunk
        i+=grain_chunk.shape[0]*grain_chunk.shape[1]

    #Mask pixels outside of bounds
    grain_size[~rfl.mask['no_data']] = -9999
    grain_size[~snow_mask] = -9999
    qa_mask = (grain_size > interp_data[1].min()) & (grain_size  < interp_data[1].max())

    temp_file =  f'temp/SISTER_{sensor}_L2B_SNOWGRAIN_{datetime}_{crid}.tif'
    grain_file =  temp_file.replace('temp','output')

    band_names = ["snowgrain_size",
                  "snowgrain_qa_mask"]

    units= ["MICRONS",
            "NA"]

    descriptions= ["SNOWGRAIN SIZE MICRONS",
                  "QUALITY ASSURANCE MASK"]

    in_file = gdal.Open(rfl.file_name)

    # Set the output raster transform and projection properties
    driver = gdal.GetDriverByName("GTIFF")
    tiff = driver.Create(temp_file,
                         rfl.columns,
                         rfl.lines,
                         2,
                         gdal.GDT_Float32)

    tiff.SetGeoTransform(in_file.GetGeoTransform())
    tiff.SetProjection(in_file.GetProjection())
    tiff.SetMetadataItem("DESCRIPTION",f"{disclaimer}SNOWGRAIN SIZE")

    # Write bands to file
    for i,band_name in enumerate(band_names,start=1):
        band = tiff.GetRasterBand(i)
        if i == 1:
            band.WriteArray(grain_size)
        else:
            band.WriteArray(qa_mask)
        band.SetDescription(band_name)
        band.SetNoDataValue(rfl.no_data)
        band.SetMetadataItem("UNITS",units[i-1])
        band.SetMetadataItem("DESCRIPTION",descriptions[i-1])
    del tiff, driver

    os.system(f"gdaladdo -minsize 900 {temp_file}")
    os.system(f"gdal_translate {temp_file} {grain_file} -co COMPRESS=LZW -co TILED=YES -co COPY_SRC_OVERVIEWS=YES")

    qlook = np.copy(grain_size)
    qlook[(qlook < interp_data[1].min()) | (qlook  > interp_data[1].max())] = 0
    qlook = (qlook-interp_data[1].min())/(interp_data[1].max()-interp_data[1].min())

    cmap = plt.get_cmap('cool')
    qlook = cmap(qlook)[:,:,:3]
    qlook = (255 * qlook).astype(np.uint8)
    qlook[grain_size == -9999] = 0

    im = Image.fromarray(qlook, 'RGB')
    im.save(grain_file.replace('tif','png'))

    shutil.copyfile(run_config_json,
                    grain_file.replace('.tif','.runconfig.json'))

    if os.path.exists("run.log"):
        shutil.copyfile('run.log',
                        grain_file.replace('.tif','.log'))

    # If experimental, prefix filenames with "EXPERIMENTAL-"
    if experimental:
        for file in glob.glob(f"output/SISTER*"):
            shutil.move(file, f"output/EXPERIMENTAL-{os.path.basename(file)}")

    # Update the path variables if now experimental
    grain_file = glob.glob("output/*%s.tif" % run_config['inputs']['crid'])[0]
    out_runconfig = glob.glob("output/*%s.runconfig.json" % run_config['inputs']['crid'])[0]
    log_path = out_runconfig.replace(".runconfig.json", ".log")
    grain_basename = os.path.basename(grain_file)[:-4]

    # Generate STAC
    catalog = pystac.Catalog(id=grain_basename,
                             description=f'{disclaimer}This catalog contains the output data products of the SISTER '
                                         f'snow grain size PGE, including a snow grain size cloud-optimized GeoTIFF. '
                                         f'Execution artifacts including the runconfig file and execution '
                                         f'log file are included with the snow grain size data.')

    # Add items for data products
    tif_files = glob.glob("output/*SISTER*.tif")
    tif_files.sort()
    description = f'{disclaimer}Snow grain size, microns'
    for tif_file in tif_files:
        metadata = generate_stac_metadata(grain_basename, description, run_config["metadata"])
        assets = {
            "cog": f"./{os.path.basename(tif_file)}",
        }
        # If it's the snow grain size product, then add png, runconfig, and log
        if os.path.basename(tif_file) == f"{grain_basename}.tif":
            png_file = tif_file.replace(".tif", ".png")
            assets["browse"] = f"./{os.path.basename(png_file)}"
            assets["runconfig"] = f"./{os.path.basename(out_runconfig)}"
            if os.path.exists(log_path):
                assets["log"] = f"./{os.path.basename(log_path)}"
        item = create_item(metadata, assets)
        catalog.add_item(item)

    # set catalog hrefs
    catalog.normalize_hrefs(f"./output/{grain_basename}")

    # save the catalog
    catalog.describe()
    catalog.save(catalog_type=pystac.CatalogType.SELF_CONTAINED)
    print("Catalog HREF: ", catalog.get_self_href())

    # Move the assets from the output directory to the stac item directories
    for item in catalog.get_items():
        for asset in item.assets.values():
            fname = os.path.basename(asset.href)
            shutil.move(f"output/{fname}", f"output/{grain_basename}/{item.id}/{fname}")


def generate_stac_metadata(basename, trait, description, in_meta):

    out_meta = {}
    out_meta['id'] = basename
    out_meta['start_datetime'] = dt.datetime.strptime(in_meta['start_time'], "%Y-%m-%dT%H:%M:%SZ")
    out_meta['end_datetime'] = dt.datetime.strptime(in_meta['end_time'], "%Y-%m-%dT%H:%M:%SZ")
    # Split corner coordinates string into list
    geometry = in_meta['bounding_box'].copy()
    # Add first coord to the end of the list to close the polygon
    geometry.append(geometry[0])
    out_meta['geometry'] = geometry
    product = basename.split('_')[3]
    if trait is not None:
        product += f"_{trait}"
    out_meta['properties'] = {
        'sensor': in_meta['sensor'],
        'description': description,
        'product': product,
        'processing_level': basename.split('_')[2]
    }
    return out_meta


def create_item(metadata, assets):
    item = pystac.Item(
        id=metadata['id'],
        datetime=metadata['start_datetime'],
        start_datetime=metadata['start_datetime'],
        end_datetime=metadata['end_datetime'],
        geometry=metadata['geometry'],
        bbox=None,
        properties=metadata['properties']
    )
    # Add assets
    for key, href in assets.items():
        item.add_asset(key=key, asset=pystac.Asset(href=href))
    return item


if __name__ == "__main__":
    main()
