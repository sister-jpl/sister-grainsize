"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus

"""

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

    os.mkdir('output')
    os.mkdir('temp')

    CRID = run_config["inputs"]["CRID"]

    rfl_base_name = os.path.basename(run_config['inputs']['l2a_rfl'])
    sister,sensor,level,product,datetime,in_CRID = rfl_base_name.split('_')
    rfl_file = f'input/{rfl_base_name}/{rfl_base_name}.bin'
    rfl_met = rfl_file.replace('.bin','.met.json')

    fc_base_name = os.path.basename(run_config['inputs']['l2b_frcov'])
    fc_file = f'input/{fc_base_name}/{fc_base_name}.tif'

    grain_met = f'output/SISTER_{sensor}_L2B_SNOWGRAIN_{datetime}_{CRID}.met.json'

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

    temp_file =  f'temp/SISTER_{sensor}_L2B_SNOWGRAIN_{datetime}_{CRID}.tif'
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
    tiff.SetMetadataItem("DESCRIPTION","SNOWGRAIN SIZE")

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

    generate_metadata(rfl_met,
                      grain_met,
                      {'product': 'SNOWGRAIN',
                      'processing_level': 'L2B',
                      'description' : 'Snow grain size, microns'})

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

    shutil.copyfile('run.log',
                    grain_file.replace('.tif','.log'))



def generate_metadata(in_file,out_file,metadata):

    with open(in_file, 'r') as in_obj:
        in_met =json.load(in_obj)

    for key,value in metadata.items():
        in_met[key] = value

    with open(out_file, 'w') as out_obj:
        json.dump(in_met,out_obj,indent=3)



if __name__ == "__main__":
    main()
