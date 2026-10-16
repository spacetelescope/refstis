"""Functions to create weekly superdarks for the STIS instrument.
"""

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
from scipy.signal import medfilt
import shutil
import argparse

from . import functions
from .functions import SIDE2_MJD_START


def create_superdark(crj_filename, basedark):
    """Create a superdark from the crj and basedark

    Create the science portion of the forthcoming reference dark by adding the
    'only baseline dark current' image to the 'only hot pixels' image.

    .. note:: input file will be updated in-place.

    Parameters
    ----------
    crj_filename : str
        filename of the cosmic-ray rejected file
    basedark : str
        basedark name
    """
    with fits.open(crj_filename, mode='update') as crj_hdu:

        # Perform iterative statistics on this normalized superdark
        data_mean, data_median, data_std = sigma_clipped_stats(crj_hdu[('sci', 1)].data,
                                                               sigma=5,
                                                               maxiters=40)

        p_five_sigma = data_median + (5 * data_std)
        print('hot pixels are defined as above: ', p_five_sigma)
        basedark_hdu = fits.open(basedark)

        base_mean, base_median, base_std = sigma_clipped_stats(basedark_hdu[('sci', 1)].data,
                                                               sigma=5,
                                                               maxiters=40)

        fivesig = base_median + 5.0 * base_std
        zerodark = crj_hdu[('sci', 1)].data - base_median
        only_hotpix = np.where(crj_hdu[('sci', 1)].data >= p_five_sigma,
                               zerodark,
                               0.0)

        basedark_med = medfilt(np.array(basedark_hdu[('sci', 1)].data, dtype=np.float32), (5, 5))
        only_dark = np.where(basedark_hdu[('sci', 1)].data >= p_five_sigma,
                             basedark_med,
                             basedark_hdu[('sci', 1)].data)

        crj_hdu[('sci', 1)].data = only_dark + only_hotpix

        # update DQ extension
        crj_hdu[('dq', 1)].data = np.where(only_hotpix >= p_five_sigma,
                                           16,
                                           crj_hdu[('dq', 1)].data)

        # Update Error
        crj_hdu[('err', 1)].data = np.where(only_hotpix == 0,
                                            basedark_hdu[('err', 1)].data,
                                            crj_hdu[('err', 1)].data)


def make_weekdark(input_list, refdark_name, thebasedark, thebiasfile=None,
                  dark_v_temp=0.07, s2ref_temp=18.):
    """Create a weekly dark reference file

    1. If not already done, run basic2d with blevcorr, biascorr, and dqicorr
       set to perform
    2. Apply temperature correction to the data
    3. split all raw images into their imsets
    4. join imsets together into a single file
    5. combine and cr-reject
    6. normalize to e/s by dividing by (exptime/gain)
    7. do hot pixel things

    .. NOTE::
       Update ERR extension of new superdark by assigning the ERR values of the
       basedark except for the new hot pixels that are updated from the weekly
       superdark, for which the error extension of the weekly superdark is taken.

    Parameters
    ----------
    input_list : list
        list of input STIS dark files
    refdark_name : str
        output name of the reference dark file
    thebasedark : str
        Monthly basedark
    thebiasfile : str, bool, optional
        biasfile to use for calibration
    dark_v_temp : float, default=0.07
        slope of temperature correction for side-2 (1/°C)
    s2ref_temp : float, default=18.
        reference temperature for side-2 (°C)
    """
    print('#-------------------------------#')
    print('#        Running weekdark       #')
    print('#-------------------------------#')
    if not thebiasfile:
        thebiasfile = fits.getval(input_list[0], 'biasfile', 0)

    print(f'Making weekdark {refdark_name}')
    print(f'With : {thebiasfile}')
    print(f'     : {thebasedark}')

    flt_list = [functions.bias_subtract_data(item, thebiasfile) for item in input_list]

    for filename in flt_list:
        texpstrt = fits.getval(filename, 'texpstrt', 0)
        # Side 1 operations ended on May 16, 2001.
        # Side 2 operations started on July 10, 2001,
        # 52091.0 corresponds to July 1, 2001
        if texpstrt > SIDE2_MJD_START:
            functions.apply_dark_correction(filename, texpstrt,
                                            dark_v_temp=dark_v_temp, s2ref_temp=s2ref_temp)

    joined_out = refdark_name.replace('.fits', '_joined.fits')
    print(f'Joining images to {joined_out}')
    functions.msjoin(flt_list, joined_out)

    crdone = functions.bd_crreject(joined_out)
    print("## crdone is ", crdone)
    if not crdone:
        functions.bd_calstis(joined_out, thebiasfile)

    crj_filename = joined_out.replace('.fits', '_crj.fits')
    shutil.copy(crj_filename, refdark_name)
    functions.normalize_crj(refdark_name)

    create_superdark(refdark_name, thebasedark)

    functions.update_header_from_input(refdark_name, input_list,
                                       dark_v_temp=dark_v_temp, s2ref_temp=s2ref_temp)
    fits.setval(refdark_name, 'TASKNAME', ext=0, value='WEEKDARK')

    print('Cleaning up...')
    functions.RemoveIfThere(crj_filename)
    functions.RemoveIfThere(joined_out)
    # map(functions.RemoveIfThere, flt_list)

    print('Weekdark done for {}'.format(refdark_name))


def call_make_weekdark():
    """Parse command line arguments and call ``weekdark``.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('input_list',
                        metavar='FILE',
                        nargs='*',
                        help='input files to turn into reference file')

    parser.add_argument('-o',
                        dest='refdark_name',
                        metavar='OUTNAME',
                        type=str,
                        default='basebias.fits',
                        help='output name for the reference file; default="basebias.fits"')

    parser.add_argument('-d',
                        dest='thebasedark',
                        metavar='DARKNAME',
                        type=str,
                        required=True,
                        help='basedark file for calibration')

    parser.add_argument('-b',
                        dest='thebiasfile',
                        metavar='BIASNAME',
                        type=str,
                        default=None,
                        help='bias file for calibration; default=None')

    parser.add_argument('-s',
                        dest='dark_v_temp',
                        type=float,
                        default=0.07,
                        help='slope of temperature correction for side-2; default=0.07 / °C')

    parser.add_argument('-r',
                        dest='s2ref_temp',
                        type=float,
                        default=18.,
                        help='reference temperature for side-2; default=18 °C')

    kwargs = vars(parser.parse_args())
    make_weekdark(**kwargs)
