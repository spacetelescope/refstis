"""

Script used to create monthly reference biases for the STIS Darks and Biases 
reference file pipeline

"""

import REFSTIS_functions
import shutil
try:
    import pyfits
except:
    from astropy.io import fits as pyfits
from scipy.signal import medfilt
import numpy as np
import support
import pdb


def flag_hot_pixels(refbias_name):
    refbias_hdu = pyfits.open(refbias_name, mode = 'update')
    bias_median_smooth = medfilt( refbias_hdu[('sci', 1)].data, (3, 15) )  #iraf used 2, but python insists that:ValueError: Each element of kernel_size should be odd.
    bias_median_smooth_med, bias_median_smooth_mean, bias_median_smooth_std = support.sigma_clip(bias_median_smooth, sigma = 3, iterations = 30)
    bias_median, bias_mean, bias_std = support.sigma_clip( refbias_hdu[ ('sci', 1) ].data , sigma=3, iterations=30 )
    diff_med = bias_mean - bias_median_smooth_mean
    bias_median_smooth = bias_median_smooth + diff_med
    bias_residual = refbias_hdu[('sci', 1)].data - bias_median_smooth
    resid_median, resid_mean, resid_std = support.sigma_clip(bias_residual, sigma = 3, iterations = 30)
    r_five_sigma = resid_mean + 5.0*resid_std
    print 'Updating DQ values of hot pixels above a level of ', r_five_sigma
    refbias_hdu[('dq', 1)].data = np.where(bias_residual > r_five_sigma, 16, refbias_hdu[('dq', 1)].data)
    pdb.set_trace()
    refbias_hdu.flush()
    refbias_hdu.close()

#-------------------------------------------------------------------------------

def make_refbias( input_list, refbias_name='refbias.fits' ):
    """ Basic premise for making a refbias:
    1- join imsets from each datset together into one large file
    2- combine and cosmic ray screen joined imset
    3- divide by number of imsets
    4- set header keywords

    """

    print '#-------------------------------#'
    print '#        Running refbias        #'
    print '#-------------------------------#'
    print 'Making refbias %s' % (refbias_name)
    joined_out = refbias_name.replace('.fits', '_joined.fits' )

    print 'Joining images to %s' % joined_out
    REFSTIS_functions.msjoin( input_list, joined_out)

    print 'Checking for cosmic ray rejection'
    crj_filename = REFSTIS_functions.crreject( joined_out )

    shutil.copy( crj_filename, refbias_name)
    flag_hot_pixels(refbias_name)
    REFSTIS_functions.update_header_from_input( refbias_name, input_list )
    
    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )
    
    print 'refbias done'

#-------------------------------------------------------------------------------
