"""
Functions to create weekly superdarks for the STIS darks and bias reference
file pipeline.

"""

import os
import shutil
import pyfits
from scipy.signal import medfilt
import numpy as np

import REFSTIS_functions

#-------------------------------------------------------------------------------

def make_weekdark( input_list, refdark_name, thebiasfile, thebasedark ):
    """
    1- split all raw images into their imsets
    2- join imsets together into a single file
    3- combine and cr-reject
    4- normalize to e/s by dividing by (exptime/gain)
    5- do hot pixel things

    """

    print '#-------------------------------#'
    print '#        Running weekdark       #'
    print '#-------------------------------#'
    print 'Making weekdark %s' % (refdark_name)
    print 'With : %s' % (thebiasfile)
    print '     : %s' % (thebasedark)
    joined_out = refdark_name.replace('.fits', '_joined.fits' )

    refdark_path = os.path.split( refdark_name )[0] or './'

    print 'Joining images to %s' % joined_out
    REFSTIS_functions.msjoin( input_list, joined_out)

    crdone = REFSTIS_functions.bd_crreject( joined_out )
    print "## crdone is ", crdone
    if (not crdone):
        REFSTIS_functions.bd_calstis(joined_out, thebiasfile)

    crj_filename = joined_out.replace('.fits', '_crj.fits')
    REFSTIS_functions.normalize_crj( crj_filename )


    #
    #***************************************************************************
    # Perform iterative statistics on the normalized superdark
    #   (i.e., neglecting hot pixels in the process)
    #
    iter_count, median, sigma, npx, med, mod, data_min, data_max = REFSTIS_functions.iterate( crj_filename )
    five_sigma = median + 5 * sigma

    # save hot pixel level and the name of the baseline dark
    # for use in updating history
    #weekoutfile = os.path.join(refdark_path,'five_sig.txt')
    #out_fd = open(weekoutfile, 'w')
    #out_fd.write(str(p_fivesig)+' '+basedark)
    #out_fd.close()

    #
    #***************************************************************************
    # Perform iterative statistics on the baseline superdark
    #   (i.e., neglecting hot pixels in the process)
    #
    # 1- norm_file - median = zerodark
    # 2- only_hotpix = 
    print "## Perform iterative statistics on the baseline superdark (thebasedark)"
    iter_count, base_median, base_sigma, npx, basemed, mod, d_min, d_max = REFSTIS_functions.iterate( crj_filename )
    five_sigma = base_median + 5 * base_sigma

    print "## Create median-filtered version of super-de-buper dark "
    basedark_hdu = pyfits.open( thebasedark )
    basedark_med = medfilt( basedark_hdu[('sci', 1)].data, (5, 5) )

    crj_hdu = pyfits.open( crj_filename )

    zerodark = crj_hdu[ ('sci', 1) ].data - basemed
    only_hotpix = np.where( basedark_hdu[ ('sci', 1) ].data >= five_sigma, zerodark, 0 )

    print "## Create 'only baseline dark current' image from these two "
    only_dark = np.where( basedark_hdu[ ('sci', 1) ].data >= five_sigma, basedark_med, basedark_hdu[ ('sci', 1) ].data )

    print "## Add 'only baseline dark current' image to 'only hot pixels' image. "
    print "## This creates the science portion of the forthcoming reference dark. "
    superdark = only_dark + only_hotpix
    crj_hdu[ ('sci', 1) ].data = superdark

    #
    #***************************************************************************
    # Update DQ extension of normalized dark by assigning the value 16 to the
    # hot pixels, and put the result in temporary image.
    #
    index = np.where( only_hotpix >= five_sigma )[0]
    crj_hdu[ ('dq', 1) ].data[index] = 16

    #***************************************************************************
    # Update ERR extension of new superdark by assigning the ERR values of the
    # basedark except for the new hot pixels that are updated from the weekly
    # superdark, for which the error extension of the weekly superdark is taken.
    # Put the result in temporary ERR image.
    #
    print "## Use imcalc to update ERR extension of new superdark "
    ### This is obviously wrong
    index = np.where( only_hotpix == 0 )
    crj_hdu[ ('err', 1) ].data[index] = basedark_hdu[ ('err', 1) ].data[index]
    index = np.where( only_hotpix != 0 )
    crj_hdu[ ('err', 1) ].data[index] = basedark_hdu[ ('err', 1) ].data[index]


    pyfits.setval(crj_filename, "FILENAME", value=crj_filename)
    pyfits.setval(crj_filename, "FILETYPE", value="DARK IMAGE")
    pyfits.setval(crj_filename, "CCDAMP", value="ANY")
    pyfits.setval(crj_filename, "CCDGAIN", value="-1")
    #pyfits.setval(crj_filename, "BINAXIS1", value=xbin)
    #pyfits.setval(crj_filename, "BINAXIS2", value=ybin)
    pyfits.setval(crj_filename, "USEAFTER", value=" ")
    pyfits.setval(crj_filename, "PEDIGREE", value="INFLIGHT")
    pyfits.setval(crj_filename, "DESCRIP",
                  value="Weekly superdark created by J. Ely")
    pyfits.setval(crj_filename, "NEXTEND", value=3)
    pyfits.setval(crj_filename, "COMMENT", 
                  value="created by the STIS weekdark task in the reference file pipeline")

    shutil.copy( crj_filename, refdark_name )

    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )

    print '#-------------------------------#'
    print '#        Finished weekdark      #'
    print '#-------------------------------#'

#-------------------------------------------------------------------------------
