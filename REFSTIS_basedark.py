"""

Functions to produce a monthly basedark for the STIS Darks and Biases
reference file pipeline

"""

import numpy as np
import pyfits
import shutil

import REFSTIS_functions

#--------------------------------------------------------------------------

def find_hotpix( filename ):
    """ find pixels hotter that median + 5*sigma and update DQ array with 
    DQ == 16

    Updates file in place.
    
    """

    iter_count, median, sigma, npx, med, mod, data_min, data_max = REFSTIS_functions.iterate( filename )
    five_sigma = median + 5*sigma

    hdu = pyfits.open( filename, mode='update' )
    index = np.where( hdu[ ('SCI', 1) ].data >= five_sigma + .1)
    hdu[ ('DQ', 1) ].data[index] = 16

    hdu.flush()
    hdu.close()

#--------------------------------------------------------------------------

def update_header( filename, xbin, ybin ):
    """ Update header information for the reference bias"""
   
    ### NOT USED

    pyfits.setval( filename, 'FILENAME', value=filename )
    pyfits.setval( filename, 'FILETYPE', value='DARK IMAGE' )
    pyfits.setval( filename, 'DETECTOR', value='CCD' )
    pyfits.setval( filename, 'CCDAMP', value='ANY' )
    pyfits.setval( filename, 'CCDGAIN', value='-1' )
    pyfits.setval( filename, 'BINAXIS1', value=xbin )
    pyfits.setval( filename, 'BINAXIS2', value=ybin )
    pyfits.setval( filename, 'USEAFTER', value='' )
    pyfits.setval( filename, 'PEDIGREE', value='INFLIGHT' )
    pyfits.setval( filename, 'DESCRIP', value='Monthly superdark created by J. Ely' )
    pyfits.setval( filename, 'NEXTEND', value=3 )
    pyfits.setval( filename, 'COMMENT', value='Reference file created by the STIS DARK reference file pipeline')

#--------------------------------------------------------------------------


def normalize_crj( filename ):
    """ Normalize the input filename by exptim/gain and flush hdu """

    hdu = pyfits.open( filename, mode='update' )

    exptime = hdu[0].header[ 'TEXPTIME' ]
    gain = hdu[0].header[ 'ATODGAIN' ]

    hdu[ ('sci', 1) ].data /= (float(exptime) / gain)

    hdu[0].header['TEXPTIME'] = 1

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

    normalize_crj( crj_filename)

    shutil.copy( crj_filename, refdark_name )

    find_hotpix( refdark_name )

    print 'Cleaning...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_filename )

    ### Do i need any of this?
    #hot_data = pyfits.getdata( norm_filename,ext=1 )
    #np.where( hot_data > 5*median_level, hot_data - median_level, 0 )

    #median_image = norm_filename + '_med.fits'
    #iraf.median( norm_filename, median_image, xwindow=2, ywindow=2,verbose=no)
    
    #only_dark = norm_filename+'_onlydark.fits'
    #med_hdu = pyfits.getdata( median_image,ext=1 )

    print '#-------------------------------#'
    print '#        Finished basedark      #'
    print '#-------------------------------#'
    
#--------------------------------------------------------------------------

    
    
