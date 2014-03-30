"""
Functions to create a weekly bias for the STIS Darks and Biases reference file
pipeline

"""

import numpy as np
import shutil

from astropy.io import fits as pyfits

import support
import functions

#-------------------------------------------------------------------------------

def make_weekbias( input_list, refbias_name, basebias ):
    """ Make 'weekly' bias from list of input bias files

    1- join imsets from each datset together into one large file
    2- combine and cosmic ray screen joined imset
    3- find hot colums
    4- add hot colums in to basebias sci data as output science data
    5- update error, dq, and headers

    # update SCI
    #***************************************************************************
    #  the baseline bias rate is taken from the aptly named basebias
    #  anything that is "hot" is assumed to be potentially transient and 
    #  is thus taken from the weekly biases.

    # update err
    #***************************************************************************
    # Update ERR extension of new superbias by assigning the ERR values of the
    # baseline superbias except for the new hot pixels that are updated from 
    # the weekly superbias, for which the error extension of the weekly
    # superbias is taken. Put the result in temporary ERR image.
    """

    print '#-------------------------------#'
    print '#        Running weekbias       #'
    print '#-------------------------------#'
    print 'Output to %s' % (refbias_name)

    joined_out = refbias_name.replace('.fits', '_joined.fits')
    functions.msjoin( input_list, joined_out)

    crj_filename = functions.crreject( joined_out )
    residual_image, median_image = functions.make_residual( crj_filename )

    residual_columns_2d = functions.make_resicols_image( residual_image )
    resi_cols_median, resi_cols_mean, resi_cols_std = support.sigma_clip( residual_columns_2d[0], sigma=3, iterations=20 )
    replval = resi_cols_mean + 3.0 * resi_cols_std
    only_hotcols = np.where( residual_columns_2d > replval, residual_image, 0 )

    #-- update sci
    hdu = pyfits.open( crj_filename, mode='update')
    hdu[ ('sci', 1) ].data += only_hotcols 

    #-- update DQ
    hot_index = np.where( only_hotcols > 0 )[0]
    hdu[ ('dq', 1) ].data[ hot_index ] = 16

    #- update ERR
    baseline_err = pyfits.getdata( basebias, ext=('err', 1) )
    no_hot_index = np.where( only_hotcols == 0 )[0]
    hdu[ ('err', 1) ].data[ no_hot_index ] = baseline_err[ no_hot_index ]

    hdu.flush()
    hdu.close()
    del hdu

    shutil.copy( crj_filename, refbias_name )
    functions.update_header_from_input( refbias_name, input_list )

    print 'Cleaning up...'
    functions.RemoveIfThere( crj_filename )
    functions.RemoveIfThere( joined_out )

    print 'weekbias done'

#-------------------------------------------------------------------------------
