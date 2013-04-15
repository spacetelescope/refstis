#!/usr/bin/env python
#--------------------------------------------------------------------------
#
# Name: basejoint
#
# Description:
#       Python translation of an original IRAF cl script to create superbias 
# reference file from a list of bias frames.
#
# METHOD:
#  The input image (with multiple extensions) is overscan-subtracted and
#  cosmic-ray-rejected using the CALSTIS algorithms within STSDAS. The
#  cosmic-ray-rejected image is divided by the number of imsets present
#  in the input image (since ocrreject adds up the individual imsets).
#  After that, the superbias is median filtered using a window of 15 x 1
#  pixels. The median-filtered is subtracted from the superbias to produce a
#  "residual" image containing hot columns and such. This "residual" image
#  is averaged along rows and replicated back to the original image size, so
#  that hot columns clearly show up. After that, the image values in
#  hot columns and pixels of the original superbias image (defined as those
#  pixels having values greater than (mean + 5 sigma of Poisson noise) are
#  replaced by those in the median-filtered bias image.
#  Plots are made of the row- and column-averaged superbias, with plotting
#  scales appropriate to the gain and binning settings of the superbias.
#
# Return:
#
# Usage:
#       % basejoint input_biases
#
# History:
# Date     OPR      Who         Reason
# -------- -------- ---------- ---------------------------------------------
# 09/25/97 ?????    P.Goudfrooij original version
# 10/28/97 ?????    P.Goudfrooij Included automated plotting of averaged columns
#                                and rows to PS files.
# 05/22/01 45222    MSwam        ported to Python/Pyraf
# 02/02/10 62549    Sherbert     remove extraneous strings, fix python
# 02/18/10 64454    Sherbert     Inc factor to 3 for min_imsets (TIR 2000-05)
# 07/01/11 68438    Sherbert     long path fix; more debugs; go absent
# 07/13/11 68438    Sherbert     python update again
#--------------------------------------------------------------------------
#
from pyraf import iraf
from iraf import stsdas,hst_calib,nicmos,stis,imgtools,ttools
from pyraf.irafglobals import *
import os
import sys
import math
import glob
import shutil
import tempfile
import pyfits
import stiref
import opusutil

import REFSTI_functions

#---------------------------------------------------------------------------

def average_biases( bias_list ):
  '''
  Create a weighted sum of the individual input files.
  First make sure all individual input files have been ocrrejected.

  returns the filename of the averaged file
  '''

  totalweight = 0
  file_path,file_name = os.path.split( bias_list[0] )
  sum_file_tmp = os.path.join( file_path, 'sum_tmp.fits' )
  sum_file = os.path.join( file_path, 'sum.fits' )
  mean_file = os.path.join( file_path, 'mean.fits' )

  for iteration,item in enumerate(bias_list):
    nimset = pyfits.getval(item,'nextend') // 3
    ncombine = pyfits.getval(item,'ncombine',ext=1)

    if (nimset > 1) | (ncombine <= 1):
      print('Input files have to be single imset files and have been CR-rejected')
      print('NIMSET: %d  NCOMBINE: %d'%(nimset,ncombine) )
      sys.exit(3)

    if (iteration==0):
      shutil.copy(item, sum_file )
    else:
      iraf.msarith( sum_file, '+', item, sum_file_tmp, verbose=0 )
      shutil.move( sum_file_tmp, sum_file )
      totalweight += ncombine

  # Then divide by the sum of the weighting factors.
  iraf.msarith( sum_file, '/', totalweight, mean_file, verbose = 0)  

  os.remove( sum_file )

  return mean_out

#---------------------------------------------------------------------------

def calibrate( input_file ):
  
    os.environ['oref'] = '/grp/hst/cdbs/oref/'

    output_blev = input_file.replace('.fits','_blev.fits')
    output_crj = input_file.replace('.fits','_crj.fits')
    # 
    # between the long file paths and not being able to find EPC files, 
    # need IRAF to run in the work dir
    # 
    #os.chdir( workdir )
    #
    # if cosmic-ray rejection has already been done on the input bias image,
    # skip all calstis-related calibration steps
    #

    fd = pyfits.open( input_file )
    nimset   = fd[0].header['nextend'] / 3
    nrptexp  = fd[0].header['nrptexp']
    crcorr   = fd[0].header['crcorr']
    blevcorr = fd[0].header['blevcorr']
    del fd

    if (nimset <= 1 and crcorr != "COMPLETE"):
        print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print("Bye now... better luck next time!")

    if (crcorr != "COMPLETE"):
       
       if (nrptexp != nimset):
            pyfits.setval(input_file,'NRPTEXP',value=nimset)
            pyfits.setval(input_file,'CRSPLIT',value=1)

       pyfits.setval(input_file, 'CRCORR', value='PERFORM')
       pyfits.setval(input_file, 'APERTURE', value='50CCD')
       pyfits.setval(input_file, 'APER_FOV', value='50x50')
       if (blevcorr != 'COMPLETE') :
           print('Performing BLEVCORR')
           pyfits.setval(input_file, 'BLEVCORR', value='PERFORM')
           iraf.basic2d(input_file, output_blev,
                        outblev = '', dqicorr = 'perform', atodcorr = 'omit',
                        blevcorr = 'perform', doppcorr = 'omit', lorscorr = 'omit',
                        glincorr = 'omit', lflgcorr = 'omit', biascorr = 'omit',
                        darkcorr = 'omit', flatcorr = 'omit', shadcorr = 'omit',
                        photcorr = 'omit', statflag = no, verb=no)
       else:
           print('Blevcorr alread Performed')
           shutil.copy(input_file,output_blev)

       print('Performing OCRREJECT')
       iraf.ocrreject(input=output_blev, output=output_crj, verb=no)

    elif (crcorr == "COMPLETE"):
        print "CR rejection already done"
        os.rename(input_file, output_crj )
  
    pyfits.setval(output_crj, 'FILENAME',
                  value=output_crj)

    return output_crj

#---------------------------------------------------------------------------

def make_basebias ( bias_list, refbias_name ):
  maxiter = 40            # Maximum number of iterations for imstat'))
  lower = INDEF           # Initial lower limit for imstat'))
  upper = INDEF           # Initial upper limit for imstat'))
  verbose = 0             # Show results of imstat iterations?'))
  PYprint = 0             # Print final results of imstat iterations?'))
  bias_path = os.path.split( bias_list[0] )[0]

  crj_list = [ calibrate(item) for item in bias_list ]

  mean_bias = average_biases( crj_list )
  bias_median = os.path.join( bias_path, 'median.fits')
  #
  #***************************************************************************
  # Median filter resulting superbias by a window of 15 x 3 pixels, and
  # subtract it from the superbias to produce a "residual" image containing
  # hot columns and such
  
  iraf.median( mean_bias + '[1]', bias_median, xwindow = 15, ywindow = 3, verb=yes)

  iraf.iterstat( mean_bias + '[1]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=yes)

  iraf.iterstat(bias_median + '[1]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=yes)

  diffmean = float(iraf.iterstat.mean) - float(iraf.iterstat.mean)
  med_hdu = pyfits.open( bias_median,mode='update' )
  med_hdu[ ('sci',1) ].data += diffmean
  med_hdu.flush()
  med_hdu.close()
  del med_hdu

  median_image = pyfits.getdata( bias_median, ext=('sci',1) )
  mean_image = pyfits.getdata( mean_bias, ext=('sci',1) )

  bias_residual = mean_image - median_image

  #
  #***************************************************************************
  # FIRST STEP TO REMOVE HOT COLUMNS:
  # Average all rows together and stretch resulting image to match input image
  #
  resi_cols = bias_path + os.sep + "resi_cols.fits"
  opusutil.RemoveIfThere(resi_cols)
  resi_cols2d = bias_path + os.sep + "resi_cols2d.fits"
  opusutil.RemoveIfThere(resi_cols2d)
  
  fd = pyfits.open( mean_image )
  xsize   = fd[1].header['naxis1']
  ysize   = fd[1].header['naxis2']
  xbin    = fd[0].header['binaxis1']
  ybin    = fd[0].header['binaxis2']
  gain    = fd[0].header['ccdgain']
  del fd
  iraf.blkavg(bias_resi, resi_cols, 1, ysize)
  iraf.blkrep(resi_cols, resi_cols2d, 1, ysize)
  
  print 'here'
  sys.exit()
  ####### Just about here

  #
  #***************************************************************************
  # Replace image values in hot columns by those in median filtered bias image
  # For now, the threshold above which a column is called "hot" is defined as
  # 3*sigma above the mean (using iterative statistics) in "resi_cols.fits".
  #
  #   replval = 0.08 / (gain*xbin)
  iraf.iterstat(resi_cols, nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
  print 'thresh mean,sigma = ',iraf.iterstat.mean,' ',iraf.iterstat.sigma
  replval = float(iraf.iterstat.mean) + 3.0 * (float(iraf.iterstat.sigma))
  iraf.imcalc(tmpsuper + '_sci.fits[0],'+bias_median+'[0],'+resi_cols2d+'[0]',
              tmpbias + '_tmp_col.fits',
              'if im3 .ge. ' + str(replval) + ' then im2 else im1', verb=no)
  
  #
  #***************************************************************************
  # SECOND STEP TO REMOVE HOT COLUMNS:
  # Average only the lower 20% of all rows together and stretch resulting
  # image to match input image.
  #
  opusutil.RemoveIfThere(tmpsuper + '_sci.fits')
  opusutil.RemoveIfThere(resi_cols)
  opusutil.RemoveIfThere(resi_cols2d)
  
  #
  # Now only use lower 20% of the Y range to check for hot columns
  #
  hotrange = int(math.floor(float(ysize) * 0.2 + 0.5))
  print 'hotrange =',hotrange
  iraf.blkavg(bias_resi+'[*,1:' + str(hotrange) + ']',
  		resi_cols, 1, hotrange)
  iraf.blkrep(resi_cols, resi_cols2d, 1, ysize)
  
  #
  #***************************************************************************
  # Replace image values in hot columns by those in median filtered bias image
  # For now, the threshold above which a column is called "hot" is defined as
  # 3*sigma above the mean (using iterative statistics) in "resi_cols.fits".
  #
  #   replval = 0.08 / (gain*xbin)
  iraf.iterstat(resi_cols, nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
  opusutil.PrintMsg('I','iraf.iterstat.mean,sigma = '+str(iraf.iterstat.mean)+' '+str(iraf.iterstat.sigma))
  replval = float(iraf.iterstat.mean) + 3.0 * (float(iraf.iterstat.sigma))

  iraf.imcalc(tmpbias + '_tmp_col.fits[0],'+bias_median+'[0],'+resi_cols2d+'[0]',
  		tmpbias + '_tmp2_col.fits', 
                'if im3 .ge. ' + str(replval) + ' then im2 else im1', verb=no)

  print('Columns hotter than '+str(replval)+' replaced by median value')
  #
  #***************************************************************************
  # Replace image values in residual single hot pixels (defined as those having
  # values greater than (mean + 5 sigma of Poisson noise) by those in
  # median-filtered bias image. This represents the science extension of the
  # final output reference superbias.
  #
  iraf.imarith(tmpbias + '_tmp2_col.fits[0]', '-', bias_median,
  		bias_resi2, verb=no)
  
  mn = -1.0
  sig = -1.0
  npx = -1
  med = -1.0
  mod = -1.0
  min = -1.0
  max = -1.0
  img = bias_resi2
  Pipe1 = iraf.imstat(img, fields = 'mean,stddev,npix,midpt,mode,min,max', 
                      lower = lower, upper = upper, PYfor=no, Stdout=1)
  parts = Pipe1[0].split( )
  mn    = float ( parts[0] )     ## string to float
  sig   = float ( parts[1] )     ## string to float
  npx   = int   ( parts[2] )     ## string to int
  med   = float ( parts[3] )     ## string to float
  mod   = float ( parts[4] )     ## string to float
  min   = float ( parts[5] )     ## string to float
  max   = float ( parts[6] )     ## string to float
  del Pipe1
  m = 1
  while (m <= maxiter):
      if (verbose):
  	opusutil.PrintMsg('I',str(m)+' '+img+': mean='+str(mn)+' rms='+str(sig))
  	opusutil.PrintMsg('I','   npix='+str(npx)+' median='+str(med)+' mode='+str(mod))
  	opusutil.PrintMsg('I','   min='+str(min)+'max='+str(max)) 
      ll = float(mn) - (5.0 * float(sig))
      ul = float(mn) + (5.0 * float(sig))
      if (lower != INDEF and ll < lower):
  	ll = lower
      if (upper != INDEF and ul > upper):
  	ul = upper
      nx = -1
      Pipe1 = iraf.imstat(img, fields = 'mean,stddev,npix,midpt,mode,min,max',
  				lower = ll, upper = ul, PYfor=no, Stdout=1)
      parts = Pipe1[0].split()
      mn    = float ( parts[0] )
      sig   = float ( parts[1] )
      nx    = int   ( parts[2] )
      med   = float ( parts[3] )
      mod   = float ( parts[4] )
      min   = float ( parts[5] )
      max   = float ( parts[6] )

      del Pipe1
      if (nx == npx):
  	break
      npx = nx
      m = m + 1

  if (PYprint and not verbose):
  	opusutil.PrintMsg('I','Median-subtracted bias: mean='+str(mn)+' rms='+str(sig))
  	opusutil.PrintMsg('I','   npix='+str(npx)+' median='+str(med)+' min='+str(min)+'max='+str(max))
  
  fivesig = float(mn) + (5.0 * float(sig))
  iraf.imcalc(tmpbias + '_tmp2_col.fits[0],'+bias_median+'[0],'+bias_resi2+'[0]',
  		tmpbias + '.fits',
  		'if im3 .ge. ' + str(fivesig) + ' then im2 else im1', verb=no)
  
  #
  #***************************************************************************
  # Build reference superbias file
  #
  # Copy "NULL reference bias" to current directory, taking into account
  #  the x- and y-binning (xsize and ysize) (the "NULL" ref. bias is one
  #  without a HISTORY keyword)
  #
  ref_template = os.path.expandvars('$oref/ref_null_bia' + str(xbin) + 'x' + str(ybin) +
                                    '.fits')
  shutil.copyfile(ref_template, thebasefile + '.fits')
 
  iraf.imcopy(tmpbias + '.fits[0]', thebasefile + '[sci,1][*,*]', verb=no)
  iraf.imcopy(tmpsuper + '.fits[err,1]', thebasefile + '[err,1][*,*]', verb=no)
  iraf.imcopy(tmpsuper + '.fits[dq,1]', thebasefile + '[dq,1][*,*]', verb=no)

  pyfits.setval(thebasefile, 'FILENAME', value=thebasefile + '.fits')
  pyfits.setval(thebasefile, 'FILETYPE', value='CCD BIAS IMAGE')
  pyfits.setval(thebasefile, 'CCDGAIN', value=gain)
  pyfits.setval(thebasefile, 'BINAXIS1', value=xbin)
  pyfits.setval(thebasefile, 'BINAXIS2', value=ybin)
  pyfits.setval(thebasefile, 'USEAFTER', value=' ')
  pyfits.setval(thebasefile, 'PEDIGREE', 'INFLIGHT')
  pyfits.setval(thebasefile, 'DESCRIP', value='Superbias created by R. de los Santos from proposals 7948/7949/8409/8439')
  pyfits.setval(thebasefile, 'NEXTEND', value='3')
  pyfits.setval(thebasefile, 'COMMENT', value='Reference file created by the STIS BIAS file pipeline')

  pyfits.setval(thebasefile, 'NCOMBINE', value=totweight, ext=1)
  

#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_basebias( glob.glob(sys.argv[1]), sys.argv[2] )
