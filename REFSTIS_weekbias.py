"""
Functions to create a weekly bias for the STIS Darks and Biases reference file
pipeline

"""

from pyraf import iraf
from iraf import stsdas,hst_calib,nicmos,stis,imgtools,mstools,ttools
from pyraf.irafglobals import *
import numpy as np
import os
import sys
import glob
import math
import shutil
import tempfile
import pyfits
import traceback

import REFSTIS_functions

#--------------------------------------------------------------------------

def make_weekbias( input_list, refbias_name, basebias ):
    """ Make 'weekly' bias from list of input bias files

    """

    from pyraf import iraf
    from iraf import stsdas,toolbox,imgtools,mstools
    import os

    print '#-------------------------------#'
    print '#        Running weekbias       #'
    print '#-------------------------------#'
    print 'Making weekbias %s' % (refbias_name)
    joined_out = refbias_name.replace('.fits', '_joined.fits')
    REFSTIS_functions.RemoveIfThere( refbias_name )

    bias_path = os.path.split( input_list[0] )[0] or './'

    print 'Joining images to %s' % joined_out
    REFSTIS_functions.msjoin( input_list, joined_out)

    crj_filename = REFSTIS_functions.crreject( joined_out )
    bias_median = os.path.join( bias_path, 'median.fits')
    REFSTIS_functions.RemoveIfThere( bias_median )

    iraf.median( crj_filename + '[1]', bias_median, xwindow = 15, ywindow = 2, verb=yes)

    iraf.iterstat( crj_filename + '[1]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
    iraf.iterstat( bias_median + '[0]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
    diffmean = float(iraf.iterstat.mean) - float(iraf.iterstat.mean)

    bias_residual = os.path.join( bias_path, 'residual.fits' )
    REFSTIS_functions.RemoveIfThere( bias_residual )
    iraf.imarith( crj_filename + '[1]', '-', bias_median + '[0]', bias_residual, verb=no)

    #
    # STEP TO IDENTIFY HOT COLUMNS:
    # Average only the lower 20% of all rows together and stretch resulting
    # image to match input image.
    #
    resi_cols = os.path.join( bias_path, "resi_cols.fits")
    REFSTIS_functions.RemoveIfThere(resi_cols)
    resi_cols2d = os.path.join( bias_path, "resi_cols2d.fits" )
    REFSTIS_functions.RemoveIfThere(resi_cols2d)
    ysize = pyfits.getval( crj_filename, 'NAXIS2',ext=1)

    #
    # Only use lower 25% of the Y range to check for hot columns
    #
    hotrange = int(math.floor(float(ysize) * 0.25 + 0.5))
    iraf.blkavg(bias_residual+'[*,1:' + str(hotrange) + ']',
                resi_cols, 1, hotrange)
    iraf.blkrep(resi_cols, resi_cols2d, 1, ysize)
    #
    #***************************************************************************
    # Determine x coordinates of "hot" columns and put them in data file.
    # For now, the threshold above which a column is called "hot" is defined as
    # 3*sigma above the mean (using iterative statistics) in "resi_cols.fits".
    #
    '''
    hotcolfile = workdir + os.sep + 'hotcols.dat'
    REFSTIS_functions.RemoveIfThere(hotcolfile)
    iraf.iterstat(resi_cols, nsigrej = 3., maxiter = 20, PYprint=no, verbose=no)
    replval = float(iraf.iterstat.mean) + 5 * (float(iraf.iterstat.sigma))

    hotcols_fd = open(hotcolfile, 'a')
    savestreams = sys.stdin, sys.stdout
    sys.stdout = hotcols_fd
    iraf.pixlocate (resi_cols+"[*,1]", lower=replval, upper=INDEF, maxvals=1024, border=0, outside=no)
    hotcols_fd.close()
    sys.stdin, sys.stdout = savestreams
    '''

    #
    #***************************************************************************
    # Create image containing only hot columns (defined as > 5 sigma above the
    # local median value) from the "residual bias" image that was created before
    #
    iraf.iterstat(resi_cols, nsigrej = 3., maxiter = 20, PYprint=no, verbose=no)
    replval = float(iraf.iterstat.mean) + 5 * (float(iraf.iterstat.sigma))

    only_hotcols = os.path.join( bias_path, 'only_hotcols.fits' )
    REFSTIS_functions.RemoveIfThere(only_hotcols)
    iraf.imcalc (resi_cols2d+"[0],"+bias_residual+"[0]", only_hotcols, 
             "if im1 .ge. "+str(replval)+" then im2 else 0.0", verb=no)
    #
    #***************************************************************************
    # Add "baseline superbias" image to "only hot columns" image.
    # This creates the science portion of the forthcoming weekly reference bias.
    #

    # update sci
    hdu = pyfits.open( crj_filename, mode='update')

    ccdgain = hdu[0].header['ATODGAIN']
    xbin = hdu[0].header['BINAXIS1']
    ybin = hdu[0].header['BINAXIS2']

    hdu[ ('sci',1) ].data += pyfits.getdata( only_hotcols, ext=0)

    # update dq
    hot_data = pyfits.getdata( only_hotcols, ext=0 )
    hot_index = np.where( hot_data > 0 )[0]
    hdu[ ('dq',1) ].data[ hot_index ] = 16

    # update err
    #***************************************************************************
    # Update ERR extension of new superbias by assigning the ERR values of the
    # baseline superbias except for the new hot pixels that are updated from 
    # the weekly superbias, for which the error extension of the weekly
    # superbias is taken. Put the result in temporary ERR image.
    baseline_err = pyfits.getdata( basebias, ext=('err',1) )
    no_hot_index = np.where( hot_data == 0 )[0]
    hdu[ ('err',1) ].data[no_hot_index] = baseline_err[no_hot_index]

    hdu.flush()
    hdu.close()
    del hdu

    pyfits.setval( crj_filename, 'FILENAME', value=crj_filename)
    pyfits.setval( crj_filename, 'FILETYPE', value='CCD BIAS IMAGE')
    pyfits.setval( crj_filename, 'CCDGAIN', value=ccdgain)
    pyfits.setval( crj_filename, 'BINAXIS1', value=xbin)
    pyfits.setval( crj_filename, 'BINAXIS2', value=ybin)
    pyfits.setval( crj_filename, 'USEAFTER', value=' ')
    pyfits.setval( crj_filename, 'PEDIGREE', value='INFLIGHT')
    pyfits.setval( crj_filename, 'DESCRIP', value='Weekly Bias')
    pyfits.setval( crj_filename, 'NEXTEND', value='3')
    pyfits.setval( crj_filename, 'COMMENT', value='Reference file created by J. Ely')

    shutil.copy( crj_filename, refbias_name )

    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )

    print '#-------------------------------#'
    print '#        Finished weekbias      #'
    print '#-------------------------------#'

#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_weekbias( glob.glob(sys.argv[1]), sys.argv[2], sys.argv[3] )
