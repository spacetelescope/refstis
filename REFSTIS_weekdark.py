"""
Functions to create weekly superdarks for the STIS darks and bias reference
file pipeline.

"""

import shutil
import pyfits
from scipy.signal import medfilt
import numpy as np

import support
import REFSTIS_functions

#-------------------------------------------------------------------------------

def create_superdark( crj_filename, basedark ):
    """

    files will be updated in place

    # Add 'only baseline dark current' image to 'only hot pixels' image. "
    # This creates the science portion of the forthcoming reference dark. "
    """
    crj_hdu = pyfits.open( crj_filename, mode='update' )
    data_median, data_mean, data_std = support.sigma_clip( crj_hdu[ ('sci', 1) ].data , sigma=3, iterations=40 )
    ### Not used:  five_sigma = data_median + 5 * data_sigma

    basedark_hdu = pyfits.open( basedark )
    basedark_med = medfilt( basedark_hdu[('sci', 1)].data, (5, 5) )
    base_median, base_mean, base_std = support.sigma_clip( basedark_hdu[ ('sci', 1) ].data, sigma=3, iterations=40 )
    five_sigma = base_median + 5 * base_std

    zerodark = crj_hdu[ ('sci', 1) ].data - base_median
    only_hotpix = np.where( basedark_hdu[ ('sci', 1) ].data >= five_sigma, zerodark, 0 )

    only_dark = np.where( basedark_hdu[ ('sci', 1) ].data >= five_sigma, basedark_med, basedark_hdu[ ('sci', 1) ].data )

    superdark_im = only_dark + only_hotpix

    crj_hdu[ ('sci', 1) ].data = superdark_im

    #- update DQ extension 
    index = np.where( only_hotpix >= five_sigma )[0]
    crj_hdu[ ('dq', 1) ].data[index] = 16

    #- Update Error
    index_to_replace = np.where( only_hotpix == 0 )
    crj_hdu[ ('err', 1) ].data[ index_to_replace ] = basedark_hdu[ ('err', 1) ].data[ index_to_replace ]

#-------------------------------------------------------------------------------

def make_weekdark( input_list, refdark_name, thebiasfile, thebasedark ):
    """
    1- split all raw images into their imsets
    2- join imsets together into a single file
    3- combine and cr-reject
    4- normalize to e/s by dividing by (exptime/gain)
    5- do hot pixel things
    6- 
    7- 

    # Update ERR extension of new superdark by assigning the ERR values of the
    # basedark except for the new hot pixels that are updated from the weekly
    # superdark, for which the error extension of the weekly superdark is taken.

    """

    print '#-------------------------------#'
    print '#        Running weekdark       #'
    print '#-------------------------------#'
    print 'Making weekdark %s' % (refdark_name)
    print 'With : %s' % (thebiasfile)
    print '     : %s' % (thebasedark)
    joined_out = refdark_name.replace('.fits', '_joined.fits' )

    print 'Joining images to %s' % joined_out
    REFSTIS_functions.msjoin( input_list, joined_out)

    crdone = REFSTIS_functions.bd_crreject( joined_out )
    print "## crdone is ", crdone
    if (not crdone):
        REFSTIS_functions.bd_calstis(joined_out, thebiasfile)

    crj_filename = joined_out.replace('.fits', '_crj.fits')
    REFSTIS_functions.normalize_crj( crj_filename )

    create_superdark( crj_filename, thebasedark )
    
    shutil.copy( crj_filename, refdark_name )

    REFSTIS_functions.update_header_from_input( refdark_name, input_list )

    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )

    print 'Weekdark done'

#-------------------------------------------------------------------------------
