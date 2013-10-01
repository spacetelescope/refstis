
#!/usr/bin/env python
#--------------------------------------------------------------------------
#
# Name: weekbias
#
# Description:
#       This script produces a bias reference file and other products.
#
# Return:
#
# Usage:
#       % weekbias
#
# History:
# Date     OPR      Who         Reason
# -------- -------- ---------- ---------------------------------------------
# 05/22/01 45222    MSwam      from stlocal.teststis.weekbias.cl
#--------------------------------------------------------------------------
#
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

def make_weekbias( bias_list, refbias_name, basebias ):

  from pyraf import iraf
  from iraf import stsdas,toolbox,imgtools,mstools
  import os
  maxiter = 30 
  version = "08Nov99"
  verbose = 1
 
  for item in bias_list:
    for imset in glob.glob( item.replace('.fits','??.fits') ):
      REFSTIS_functions.RemoveIfThere( imset )
 
  bias_path = os.path.split( bias_list[0] )[0]

  refbias_name.replace('.fits','') # ensure its just a rootname
  out_ref = os.path.join( bias_path, refbias_name+'.fits' )
  REFSTIS_functions.RemoveIfThere( out_ref )


  print 'Splitting images'
  nimsets = REFSTIS_functions.split_images( bias_list )
  
  print 'Joining images'
  msjoin_list = ','.join( [ item for item in glob.glob('*raw??.fits') ] )# if item[:9] in bias_list] )
  print msjoin_list
  joined_out = refbias_name+ '_joined' +'.fits' 
  REFSTIS_functions.RemoveIfThere( joined_out )
  print joined_out

  msjoin_file = open('msjoin.txt','w')
  msjoin_file.write( '\n'.join(msjoin_list.split(',')) )
  msjoin_file.close()
    
  iraf.msjoin( inimg='@msjoin.txt', outimg=joined_out, Stderr='dev$null')
    
  crj_filename = REFSTIS_functions.crreject( joined_out )
  mean_bias = crj_filename
  bias_median = os.path.join( bias_path, 'median.fits')
  REFSTIS_functions.RemoveIfThere( bias_median )

  iraf.median( mean_bias + '[1]', bias_median, xwindow = 15, ywindow = 2, verb=yes)

  iraf.iterstat( mean_bias + '[1]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
  iraf.iterstat( bias_median + '[0]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
  diffmean = float(iraf.iterstat.mean) - float(iraf.iterstat.mean)

  bias_residual = os.path.join( bias_path, 'residual.fits' )
  REFSTIS_functions.RemoveIfThere( bias_residual )
  iraf.imarith( mean_bias + '[1]', '-', bias_median + '[0]', bias_residual, verb=no)

#
# STEP TO IDENTIFY HOT COLUMNS:
# Average only the lower 20% of all rows together and stretch resulting
# image to match input image.
#
  resi_cols = os.path.join( bias_path, "resi_cols.fits")
  REFSTIS_functions.RemoveIfThere(resi_cols)
  resi_cols2d = os.path.join( bias_path, "resi_cols2d.fits" )
  REFSTIS_functions.RemoveIfThere(resi_cols2d)
  ysize = pyfits.getval( mean_bias, 'NAXIS2',ext=1)

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
  hdu = pyfits.open( mean_bias, mode='update')
  
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




#***************************************************************************
# Build reference superbias file
#
# Copy "NULL reference bias" to current directory, taking into account
#  the x- and y-binning (xsize and ysize) (the "NULL" ref. bias is one
#  without a HISTORY keyword)
#
  #ref_template = os.path.expandvars('$oref/ref_null_bia' + str(xbin) + 'x' + str(ybin) + '.fits')
  #shutil.copyfile(ref_template, thereffile + '.fits')

  

  # fill header keywords
  pyfits.setval(mean_bias, 'FILENAME', value=mean_bias)
  pyfits.setval(mean_bias, 'FILETYPE', value='CCD BIAS IMAGE')
  pyfits.setval(mean_bias, 'CCDGAIN', value=ccdgain)
  pyfits.setval(mean_bias, 'BINAXIS1', value=xbin)
  pyfits.setval(mean_bias, 'BINAXIS2', value=ybin)
  pyfits.setval(mean_bias, 'USEAFTER', value=' ')
  pyfits.setval(mean_bias, 'PEDIGREE', value='INFLIGHT')
  pyfits.setval(mean_bias, 'DESCRIP', value='Superbias created by R. de los Santos from proposals 7949/8409/8439')
  pyfits.setval(mean_bias, 'NEXTEND', value='3')
  pyfits.setval(mean_bias, 'COMMENT', value='Reference file created by J. Ely')


#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_weekbias( glob.glob(sys.argv[1]), sys.argv[2], sys.argv[3] )
