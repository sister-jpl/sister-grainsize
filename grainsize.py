"""
SISTER
Space-based Imaging Spectroscopy and Thermal PathfindER
Author: Adam Chlus

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import os
import hytools_lite as htl
from hytools_lite.io.envi import WriteENVI
import numpy as np
from scipy.interpolate import interp1d

def progbar(curr, total, full_progbar = 100):
    '''Display progress bar.
    Gist from:
    https://gist.github.com/marzukr/3ca9e0a1b5881597ce0bcb7fb0adc549

    Args:
        curr (int, float): Current task level.
        total (int, float): Task level at completion.
        full_progbar (TYPE): Defaults to 100.
    Returns:
        None.
    '''
    frac = curr/total
    filled_progbar = round(frac*full_progbar)
    print('\r', '#'*filled_progbar + '-'*(full_progbar-filled_progbar), '[{:>7.2%}]'.format(frac), end='')

def main():
    ''' Estimate snow grain size from reflectance data using method of Nolin and Dozier (2000).

    Nolin, A. W., & Dozier, J. (2000).
    A hyperspectral method for remotely sensing the grain size of snow.
    Remote sensing of Environment, 74(2), 207-216.
    https://doi.org/10.1016/S0034-4257(00)00111-5

    '''

    desc = "Estimate snow grain size"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('input', type=str,
                        help='Input reflectance image')
    parser.add_argument('out_dir', type=str,
                        help='Output directory')
    parser.add_argument('--root_dir', action='store',
                        help='Root of source directory',
                        default=os.path.realpath(os.path.split(__file__)[0]))
    parser.add_argument('--verbose',      action='store_true')
    parser.add_argument('--mask', type=str, default=None,
                        help='External snow mask')

    args = parser.parse_args()

    #Testing
    # parser = argparse.ArgumentParser()
    # args = parser.parse_args([])
    # args.input = '/data2/prisma/l2/PRS_20210302180507_20210302180511_0001_srtm/PRS_20210302180507_20210302180511_0001_rfl_prj_rfl_10nm'
    # args.out_dir ='/data2/temp/snow_testing/'
    # args.root_dir = '/home/adam/Dropbox/rs/sister/repos/sister-grainsize/'
    # args.mask = '/data2/prisma/l2/PRS_20210302180507_20210302180511_0001_srtm/PRS_20210302180507_20210302180511_0001_rdn_prj_cls'
    # args.verbose = True

    # Load band depth and snow size data
    interp_data = np.loadtxt("%s/data/grainsize_v_bandarea_nolin_dozier_2000_interp.csv" % args.root_dir,
                             delimiter = ',').T
    interpolator = interp1d(interp_data[0],interp_data[1],
                            kind = 'cubic',fill_value ='extrapolate')

    img = htl.HyTools()
    img.read_file(args.input,'envi')

    # Use external mask
    if args.mask:
        msk = htl.HyTools()
        msk.read_file(args.mask,'envi')
        mask = msk.get_band(0)==2
    # Use NDSI threshold for mask
    else:
        mask = img.ndi(670,1500) < .9
        mask[~img.mask['no_data']] =False

    # Get continuum start and end bands
    wave1, wave2  = 950,1080
    wave1_ind = img.wave_to_band(wave1)
    wave2_ind = img.wave_to_band(wave2)

    grain_size = np.full((img.lines,img.columns),-1)
    iterator =img.iterate(by = 'chunk',chunk_size = (200,200))
    i = 0

    while not iterator.complete:

        chunk = iterator.read_next()[:,:,wave1_ind:wave2_ind+1]

        # Calculate continuum removed spectrum
        slope = (chunk[:,:,-1]-chunk[:,:,0])/(wave2-wave1)
        intercept = chunk[:,:,-1]- wave2*slope
        continuum = np.einsum('i,mj->mji', img.wavelengths[wave1_ind:wave2_ind+1],slope) + np.expand_dims(intercept,axis=2)

        depth  = (continuum-chunk)/continuum
        band_area =np.trapz(depth,axis=2)

        grain_chunk = interpolator(band_area)

        grain_size[iterator.current_line:iterator.current_line+chunk.shape[0],
              iterator.current_column:iterator.current_column+chunk.shape[1]] = grain_chunk
        i+=grain_chunk.shape[0]*grain_chunk.shape[1]
        if args.verbose:
            progbar(i,img.lines*img.columns, full_progbar = 100)

    print('\n')

    #Mask pixels outside of bounds
    grain_size[(grain_size < 60) | (grain_size >= 1000)] = 0
    grain_size[~img.mask['no_data']] = -9999
    grain_size[mask] = -9999

    # Export grain size map
    grain_header = img.get_header()
    grain_header['bands']= 1
    grain_header['band names']= ['grain_size_um']
    grain_header['wavelength']= []
    grain_header['fwhm']= []

    out_file = args.out_dir + img.base_name + '_grainsize'
    writer = WriteENVI(out_file,grain_header)
    writer.write_band(grain_size,0)

if __name__ == "__main__":
    main()
