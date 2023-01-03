"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus

"""

import argparse
import os
import hytools_lite as htl
from hytools_lite.io.envi import WriteENVI
import numpy as np
from scipy.interpolate import interp1d

def main():

    ''' Estimate snow grain size from reflectance data using method of Nolin and Dozier (2000).

    Nolin, A. W., & Dozier, J. (2000).
    A hyperspectral method for remotely sensing the grain size of snow.
    Remote sensing of Environment, 74(2), 207-216.
    https://doi.org/10.1016/S0034-4257(00)00111-5

    '''

    desc = "Estimate snow grain size"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('rfl_file', type=str,
                        help='Input reflectance image')
    parser.add_argument('frac_cover', type=str, default=None,
                        help='Fractional cover map')
    parser.add_argument('out_dir', type=str,
                        help='Output directory')
    parser.add_argument('--frac_threshold', type=float, default=0.6,
                        help='Snow fractional cover threshold')
    args = parser.parse_args()


    reflectance = htl.HyTools()
    reflectance.read_file(args.rfl_file,'envi')

    #Check if input wavelengths cover absorption feat
    if (reflectance.wavelengths.min() <= 940) & (reflectance.wavelengths.max() >= 1090):

        # Load band depth and snow size data
        root_dir =os.path.realpath(os.path.split(__file__)[0])
        interp_data = np.loadtxt("%s/data/grainsize_v_bandarea_nolin_dozier_2000_interp.csv" % root_dir,
                                 delimiter = ',').T
        interpolator = interp1d(interp_data[0],interp_data[1],
                                kind = 'cubic',fill_value ='extrapolate')

        # Use external mask
        if args.frac_cover:
            msk = htl.HyTools()
            msk.read_file(args.frac_cover,'envi')
            mask = msk.get_band(3) < args.frac_threshold

        # Use NDSI threshold for mask
        else:
            mask = reflectance.ndi(670,1500) < .9

        mask[~reflectance.mask['no_data']] =False

        # Get continuum start and end bands
        wave1, wave2  = 950,1080
        wave1_ind = reflectance.wave_to_band(wave1)
        wave2_ind = reflectance.wave_to_band(wave2)

        grain_size = np.full((reflectance.lines,reflectance.columns),-1)
        iterator =reflectance.iterate(by = 'chunk',chunk_size = (200,200))
        i = 0

        while not iterator.complete:
            chunk = iterator.read_next()[:,:,wave1_ind:wave2_ind+1]

            # Calculate continuum removed spectrum
            slope = (chunk[:,:,-1]-chunk[:,:,0])/(wave2-wave1)
            intercept = chunk[:,:,-1]- wave2*slope
            continuum = np.einsum('i,mj->mji', reflectance.wavelengths[wave1_ind:wave2_ind+1],slope) + np.expand_dims(intercept,axis=2)

            depth  = (continuum-chunk)/continuum
            band_area =np.trapz(depth,axis=2)

            grain_chunk = interpolator(band_area)

            grain_size[iterator.current_line:iterator.current_line+chunk.shape[0],
                  iterator.current_column:iterator.current_column+chunk.shape[1]] = grain_chunk
            i+=grain_chunk.shape[0]*grain_chunk.shape[1]

        #Mask pixels outside of bounds
        grain_size[(grain_size < 60) | (grain_size > 1000)] = 0
        grain_size[~reflectance.mask['no_data']] = -9999
        grain_size[mask] = -9999

        # Export grain size map
        grain_header = reflectance.get_header()
        grain_header['description']= 'Snow grain size map'
        grain_header['bands']= 1
        grain_header['band names']= ['grain_size_um']
        grain_header['wavelength']= None
        grain_header['fwhm']= None
        grain_header['wavelength units']= None

        out_file = args.out_dir + reflectance.base_name + '_grainsize'
        writer = WriteENVI(out_file,grain_header)
        writer.write_band(grain_size,0)

    else:
        print("Wavelength range does not cover ice absorption feature.")

if __name__ == "__main__":
    main()
