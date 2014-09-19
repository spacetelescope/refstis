"""

Script used to create monthly reference biases for the STIS Darks and Biases 
reference file pipeline

"""

from astropy.io import fits
import numpy as np
from scipy.signal import medfilt
import shutil

import support
import functions

def flag_hot_pixels(refbias_name):
    """Flag hotpixels in the DQ array

    Parameters:
    -----------
    refbias_name : str
        name of the reference file to flag

    """
    
    with fits.open(refbias_name, mode = 'update') as refbias_hdu:
        #--iraf used  a 2 pixel binning, but python insists that:
        #--  ValueError: Each element of kernel_size should be odd.
        smooth_bias = medfilt(refbias_hdu[('sci', 1)].data, 
                                     (3, 15))
        smooth_bias_med, smooth_bias_mean, smooth_bias_std = support.sigma_clip(smooth_bias, sigma = 3, iterations = 30)

        bias_median, bias_mean, bias_std = support.sigma_clip(refbias_hdu[('sci', 1)].data , sigma=3, iterations=30)

        diff_med = bias_mean - smooth_bias_mean
        smooth_bias = smooth_bias + diff_med

        bias_residual = refbias_hdu[('sci', 1)].data - smooth_bias
        resid_median, resid_mean, resid_std = support.sigma_clip(bias_residual,
                                                                 sigma = 3, 
                                                                 iterations = 30)
        r_five_sigma = resid_mean + 5.0 * resid_std

        print 'Updating DQ values of hot pixels above a level of ', r_five_sigma
        refbias_hdu[('dq', 1)].data = np.where(bias_residual > r_five_sigma,
                                               16, 
                                               refbias_hdu[('dq', 1)].data)

#-------------------------------------------------------------------------------

def make_refbias(input_list, refbias_name='refbias.fits'):
    """ Create a refbias

    Basic premise for making a refbias
    1- join imsets from each datset together into one large file
    2- combine and cosmic ray screen joined imset
    3- divide by number of imsets
    4- set header keywords

    Parameters:
    -----------
    input_list : list
        list of input bias files
    refbias_name : str
        name of the output bias reference file

    """

    print '#-------------------------------#'
    print '#        Running refbias        #'
    print '#-------------------------------#'
    print 'Making refbias %s' % (refbias_name)
    joined_out = refbias_name.replace('.fits', '_joined.fits')

    print 'Joining images to %s' % joined_out
    functions.msjoin(input_list, joined_out)

    print 'Checking for cosmic ray rejection'
    crj_filename = functions.crreject(joined_out)

    shutil.copy(crj_filename, refbias_name)
    flag_hot_pixels(refbias_name)
    functions.update_header_from_input(refbias_name, input_list)
    
    print 'Cleaning up...'
    functions.RemoveIfThere(crj_filename)
    functions.RemoveIfThere(joined_out)
    
    print 'refbias done for {}'.format(refbias_name)

#-------------------------------------------------------------------------------
