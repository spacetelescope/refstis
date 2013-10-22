"""

Functions to produce a monthly basedark for the STIS Darks and Biases
reference file pipeline

"""

import numpy as np
try:
    import pyfits
except:
    from astropy.io import fits as pyfits

import shutil

import REFSTIS_functions
import support

#--------------------------------------------------------------------------

def find_hotpix( filename ):
    """ find pixels hotter that median + 5*sigma and update DQ array with 
    DQ == 16

    Updates file in place.
    
    """

    hdu = pyfits.open( filename, mode='update' )

    data_median, data_mean, data_std = support.sigma_clip( hdu[ ('sci', 1) ].data, 
                                                           sigma=3, 
                                                           iterations=40 )

    five_sigma = data_median + 5 * data_std
    index = np.where( hdu[ ('SCI', 1) ].data >= five_sigma + .1 )
    hdu[ ('DQ', 1) ].data[index] = 16

    hdu.flush()
    hdu.close()

#--------------------------------------------------------------------------

def make_basedark( input_list, refdark_name='basedark.fits', bias_file=None ):
    """
    Make a monthly baseline dark from the input list of raw dark files and
    the given bias file.

    1- Join all imsets from input list into single file
    2- combine and cr-reject
    3- normalize to e/s by dividing by (exptime/gain)
    4- update DQ array with hot pixel information

    """

    print '#-------------------------------#'
    print '#        Running basedark       #'
    print '#-------------------------------#'
    print 'output to: %s' % refdark_name
    print 'with biasfile %s' % bias_file

    joined_filename = refdark_name.replace('.fits', '_joined.fits') 
    crj_filename = joined_filename.replace('.fits', '_crj.fits')

    if not bias_file: raise IOError( 'No biasfile specified, this task needs one to run' )

    print 'Joining images'
    REFSTIS_functions.msjoin( input_list, joined_filename  )

    print 'CRREJECT'
    crdone = REFSTIS_functions.bd_crreject( joined_filename )
    if (not crdone):
        REFSTIS_functions.bd_calstis( joined_filename, bias_file )

    REFSTIS_functions.normalize_crj( crj_filename)

    shutil.copy( crj_filename, refdark_name )

    find_hotpix( refdark_name )

    REFSTIS_functions.update_header_from_input( refdark_name, input_list )

    print 'Cleaning...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_filename )

    print 'basedark done'
    
#--------------------------------------------------------------------------

    
    
