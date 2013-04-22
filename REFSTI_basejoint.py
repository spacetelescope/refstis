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
from iraf import stsdas,hst_calib,nicmos,stis,imgtools,mstools,ttools
from pyraf.irafglobals import *
import os
import sys
import math
import glob
import shutil
import tempfile
import pyfits
import numpy as np
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

  return mean_file,totalweight

#---------------------------------------------------------------------------

def calibrate( input_file ):
    os.environ['oref'] = '/grp/hst/cdbs/oref/'
    print 'Calibrating %s'%(input_file)
    output_blev = input_file.replace('.fits','_blev.fits')
    REFSTI_functions.RemoveIfThere( output_blev )
    output_crj = input_file.replace('.fits','_crj.fits')
    REFSTI_functions.RemoveIfThere( output_crj )

    # between the long file paths and not being able to find EPC files, 
    # need IRAF to run in the work dir
    # 
    #os.chdir( workdir )
    #
    # if cosmic-ray rejection has already been done on the input bias image,


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
        return None

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
                        photcorr = 'omit', statflag = no, verb=no, Stdout='dev$null')
       else:
           print('Blevcorr alread Performed')
           shutil.copy(input_file,output_blev)

       print('Performing OCRREJECT')
       iraf.ocrreject(input=output_blev, output=output_crj, verb=no, Stdout='dev$null')

    elif (crcorr == "COMPLETE"):
        print "CR rejection already done"
        os.rename(input_file, output_crj )
  
    pyfits.setval(output_crj, 'FILENAME',
                  value=output_crj)

    os.remove( output_blev )

    return output_crj

#---------------------------------------------------------------------------

def make_basebias ( bias_list, refbias_name ):
  maxiter = 40            # Maximum number of iterations for imstat'))
  lower = INDEF           # Initial lower limit for imstat'))
  upper = INDEF           # Initial upper limit for imstat'))
  verbose = 0             # Show results of imstat iterations?'))
  PYprint = 0             # Print final results of imstat iterations?'))
  bias_path = os.path.split( bias_list[0] )[0]

  print 'Processing individual files'
  crj_list = [ calibrate(item) for item in bias_list ]
  crj_list = [ item for item in crj_list if item != None ]

  mean_bias,totalweight = average_biases( crj_list )
  bias_median = os.path.join( bias_path, 'median.fits')
  #
  #***************************************************************************
  # Median filter resulting superbias by a window of 15 x 3 pixels, and
  # subtract it from the superbias to produce a "residual" image containing
  # hot columns and such
  print 'Median filtering'
  REFSTI_functions.RemoveIfThere( bias_median )
  iraf.median( mean_bias + '[1]', bias_median, xwindow = 15, ywindow = 3, verb=yes)

  iraf.iterstat( mean_bias + '[1]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)
  iraf.iterstat( bias_median + '[0]', nsigrej = 3., maxiter = 40, PYprint=no, verbose=no)

  diffmean = float(iraf.iterstat.mean) - float(iraf.iterstat.mean)
  med_hdu = pyfits.open( bias_median,mode='update' )
  med_hdu[0].data += diffmean
  med_hdu.flush()
  med_hdu.close()
  del med_hdu

  #median_image = pyfits.getdata( bias_median, ext=0 )
  #mean_image = pyfits.getdata( mean_bias, ext=('sci',1) )
  #bias_residual = mean_image - median_image

  bias_residual = os.path.join( bias_path, 'residual.fits' )
  REFSTI_functions.RemoveIfThere( bias_residual )
  iraf.imarith( mean_bias + '[1]', '-', bias_median + '[0]', bias_residual, verb=no)

  # FIRST STEP TO REMOVE HOT COLUMNS:
  # Average all rows together and stretch resulting image to match input image
 
  resi_cols = os.path.join( bias_path, "resi_cols.fits")
  REFSTI_functions.RemoveIfThere(resi_cols)
  resi_cols2d = os.path.join( bias_path, "resi_cols2d.fits" )
  REFSTI_functions.RemoveIfThere(resi_cols2d)
  
  fd = pyfits.open( mean_bias )
  xsize   = fd[1].header['naxis1']
  ysize   = fd[1].header['naxis2']
  xbin    = fd[0].header['binaxis1']
  ybin    = fd[0].header['binaxis2']
  gain    = fd[0].header['ccdgain']
  del fd
  iraf.blkavg(bias_residual, resi_cols, 1, ysize)
  iraf.blkrep(resi_cols, resi_cols2d, 1, ysize)

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
  tmp_bias = os.path.join( bias_path, 'tmp_col.fits' )
  REFSTI_functions.RemoveIfThere( tmp_bias )
  #iraf.imcalc(mean_bias + '[1],'+bias_median+'[0],'+resi_cols2d+'[0]', tmp_bias,
  #            'if im3 .ge. ' + str(replval) + ' then im2 else im1', verb=no)

  hdu = pyfits.open( mean_bias )
  index = np.where( pyfits.getdata( resi_cols2d,ext=0 ) >= replval )
  hdu[ ('sci',1) ].data[index] = pyfits.getdata( bias_median, ext=0)[index]
  hdu.writeto(tmp_bias)
  
  #
  #***************************************************************************
  # SECOND STEP TO REMOVE HOT COLUMNS:
  # Average only the lower 20% of all rows together and stretch resulting
  # image to match input image.
  
  REFSTI_functions.RemoveIfThere(resi_cols)
  REFSTI_functions.RemoveIfThere(resi_cols2d)

  #
  # Now only use lower 20% of the Y range to check for hot columns
  #
  hotrange = int(math.floor(float(ysize) * 0.2 + 0.5))
  print 'hotrange =',hotrange
  iraf.blkavg(bias_residual+'[*,1:' + str(hotrange) + ']',
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
  print('iraf.iterstat.mean,sigma = '+str(iraf.iterstat.mean)+' '+str(iraf.iterstat.sigma))
  replval = float(iraf.iterstat.mean) + 3.0 * (float(iraf.iterstat.sigma))

  tmp_bias2 = os.path.join( bias_path, 'tmp_col2.fits' )
  REFSTI_functions.RemoveIfThere( tmp_bias2 )

  #iraf.imcalc(tmp_bias+'[1],'+bias_median+'[0],'+resi_cols2d+'[0]', tmp_bias2, 
  #            'if im3 .ge. ' + str(replval) + ' then im2 else im1', verb=yes)

  hdu = pyfits.open( tmp_bias )
  index = np.where( pyfits.getdata( resi_cols2d,ext=0 ) >= replval )
  hdu[ ('sci',1) ].data[index] = pyfits.getdata( bias_median, ext=0)[index]
  hdu.writeto(tmp_bias2)

  print('Columns hotter than '+str(replval)+' replaced by median value')
  #
  #***************************************************************************
  # Replace image values in residual single hot pixels (defined as those having
  # values greater than (mean + 5 sigma of Poisson noise) by those in
  # median-filtered bias image. This represents the science extension of the
  # final output reference superbias.
  #
  bias_residual_2 = os.path.join( bias_path, 'residual2.fits' )
  REFSTI_functions.RemoveIfThere( bias_residual_2 )
  print 'Imarith',tmp_bias2 + '[1]', '-', bias_median+'[0]'
  iraf.imarith(tmp_bias2 + '[1]', '-', bias_median+'[0]',
  		bias_residual_2, verb=yes)

  iter_count,mn,sig,npx,med,mod,min,max = REFSTI_functions.iterate( bias_residual_2+'[0]' )

  if (PYprint and not verbose):
  	print('Median-subtracted bias: mean='+str(mn)+' rms='+str(sig))
  	ptin('   npix='+str(npx)+' median='+str(med)+' min='+str(min)+'max='+str(max))
  
  fivesig = float(mn) + (5.0 * float(sig))
  
  #iraf.imcalc(tmp_bias + '[1],'+bias_median+'[0],'+bias_residual_2+'[0]',tmp_bias,
  # 		'if im3 .ge. ' + str(fivesig) + ' then im2 else im1', verb=no)
  
  hdu = pyfits.open( tmp_bias )
  index = np.where( pyfits.getdata( bias_residual_2,ext=0 ) >= replval )
  hdu[ ('sci',1) ].data[index] = pyfits.getdata( bias_median, ext=0)[index]
  hdu.writeto(tmp_bias,clobber=True)


  #
  #***************************************************************************
  # Build reference superbias file
  #
  # Copy "NULL reference bias" to current directory, taking into account
  #  the x- and y-binning (xsize and ysize) (the "NULL" ref. bias is one
  #  without a HISTORY keyword)
  #
  #ref_template = os.path.expandvars('$oref/ref_null_bia' + str(xbin) + 'x' + str(ybin) + '.fits')
  #print ref_template
  #shutil.copyfile(ref_template, thebasefile + '.fits')
 
  out_ref = os.path.join( bias_path, refbias_name)
  shutil.copy( mean_bias, out_ref )
  ref_hdu = pyfits.open( out_ref,mode='update')
  ref_hdu[1].data = pyfits.getdata( tmp_bias, ext=0 )
  ref_hdu.flush()
  ref_hdu.close()
  del ref_hdu

  pyfits.setval( out_ref, 'FILENAME', value=out_ref)
  pyfits.setval( out_ref, 'FILETYPE', value='CCD BIAS IMAGE')
  pyfits.setval( out_ref, 'CCDGAIN', value=gain)
  pyfits.setval( out_ref, 'BINAXIS1', value=xbin)
  pyfits.setval( out_ref, 'BINAXIS2', value=ybin)
  pyfits.setval( out_ref, 'USEAFTER', value=' ')
  pyfits.setval( out_ref, 'PEDIGREE', value='INFLIGHT')
  pyfits.setval( out_ref, 'DESCRIP', value='Superbias created by R. de los Santos from proposals 7948/7949/8409/8439')
  pyfits.setval( out_ref, 'NEXTEND', value='3')
  pyfits.setval( out_ref, 'COMMENT', value='Reference file created by the STIS BIAS file pipeline')
  pyfits.setval( out_ref, 'NCOMBINE', value=totalweight, ext=1)
  

#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_basebias( glob.glob(sys.argv[1]), sys.argv[2] )
