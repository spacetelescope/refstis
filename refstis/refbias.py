"""Procedure to create a suberbias using the 'refbias' method.

Basic premise for making a refbias.

1. join imsets from each datset together into one large file
2. combine, cosmic ray screen, and normalize to number of imsets
3. set header keywords
4. clean up intermediate files

"""

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
from scipy.signal import medfilt
import shutil

from . import functions

#-------------------------------------------------------------------------------

def flag_hot_pixels(refbias_name):
    """Flag hotpixels in the DQ array

    Pixels more than (mean + 3*sigma) away from a median smoothed image
    are flagged as DQ=16 in the DQ array.

    .. note:: The input file is updated in-place.

    .. note::
      The IRAF version of this pipeline specified a 2x15 pixel median filter
      to calculate the smoothed imaged, but the IRAF task documentation says that
      even sizes are increased by 1 for the computation.  A 3x5 pixel filter is
      used directly here.

    Parameters
    ----------
    refbias_name : str
        name of the reference file to flag

    """

    with fits.open(refbias_name, mode='update') as refbias_hdu:
        smooth_bias = medfilt(np.array(refbias_hdu[('sci', 1)].data, dtype=np.float32), (3, 15))

        smooth_bias_mean, smooth_bias_med, smooth_bias_std = sigma_clipped_stats(smooth_bias, sigma=3, maxiters=30)
        bias_mean, bias_median, bias_std = sigma_clipped_stats(refbias_hdu[('sci', 1)].data, sigma=3, maxiters=30)


        smooth_bias += (bias_mean - smooth_bias_mean)

        bias_residual = refbias_hdu[('sci', 1)].data - smooth_bias

        resid_mean, resid_median, resid_std = sigma_clipped_stats(bias_residual,
                                                               sigma=3,
                                                               maxiters=30)
        r_five_sigma = resid_mean + 5.0 * resid_std

        print('Updating DQ values of hot pixels above a level of ', r_five_sigma)
        refbias_hdu[('dq', 1)].data = np.where(bias_residual > r_five_sigma,
                                               16,
                                               refbias_hdu[('dq', 1)].data)

#-------------------------------------------------------------------------------

def make_refbias(input_list, refbias_name='refbias.fits'):
    """Create a refbias FITS file

    Parameters
    ----------
    input_list : list
        list of input bias files
    refbias_name : str
        name of the output bias reference file

    """

    print('#-------------------------------#')
    print('#        Running refbias        #')
    print('#-------------------------------#')
    print('Making refbias %s' % (refbias_name))
    joined_out = refbias_name.replace('.fits', '_joined.fits')

    print('Joining images to %s' % joined_out)
    functions.msjoin(input_list, joined_out)

    print('Checking for cosmic ray rejection')
    crj_filename = functions.crreject(joined_out)

    shutil.copy(crj_filename, refbias_name)
    flag_hot_pixels(refbias_name)
    functions.update_header_from_input(refbias_name, input_list)
    fits.setval(refbias_name, 'TASKNAME', ext=0, value='refbias')

    print('Cleaning up...')
    functions.RemoveIfThere(crj_filename)
    functions.RemoveIfThere(joined_out)

    print('refbias done for {}'.format(refbias_name))

#-------------------------------------------------------------------------------
