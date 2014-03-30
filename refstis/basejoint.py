"""
Script to produce a monthly base bias for STIS CCD

 Description:
        Python translation of a python translation
  of an original IRAF cl script to create superbias 
  reference file from a list of bias frames.
 
  METHOD:
   The input image (with multiple extensions) is overscan-subtracted and
   cosmic-ray-rejected using the CALSTIS algorithms within STSDAS. The
   cosmic-ray-rejected image is divided by the number of imsets present
   in the input image (since ocrreject adds up the individual imsets).
   After that, the superbias is median filtered using a window of 15 x 1
   pixels. The median-filtered is subtracted from the superbias to produce a
   "residual" image containing hot columns and such. This "residual" image
   is averaged along rows and replicated back to the original image size, so
   that hot columns clearly show up. After that, the image values in
   hot columns and pixels of the original superbias image (defined as those
   pixels having values greater than (mean + 5 sigma of Poisson noise) are
   replaced by those in the median-filtered bias image.
   Plots are made of the row- and column-averaged superbias, with plotting
   scales appropriate to the gain and binning settings of the superbias.
 

"""

import numpy as np
import os
import sys
import shutil
from astropy.io import fits as pyfits

import support
import functions

from pyraf import iraf
from iraf import stsdas, hst_calib, stis
from pyraf.irafglobals import *


#---------------------------------------------------------------------------

def average_biases( bias_list ):
    """
    Create a weighted sum of the individual input files.
    First make sure all individual input files have been ocrrejected.

    returns the filename of the averaged file

    """
    
    assert len(bias_list), 'Bias list is empty'

    file_path, file_name = os.path.split( bias_list[0] )
    mean_file = os.path.join( file_path, 'mean.fits' )

    for iteration, item in enumerate(bias_list):
        ofile = pyfits.open(item)
        hdr0 = ofile[0].header
        hdr1 = ofile[1].header
        nimset = hdr0['nextend'] // 3
        ncombine = hdr1['ncombine']

        #If input files have more than one imset or have not been cr-rejected, exit
        if (nimset > 1) | (ncombine <= 1):
            print('Input files have to be single imset files and have been CR-rejected')
            print('NIMSET: %d  NCOMBINE: %d'%(nimset, ncombine) )
            sys.exit(3)

        #Otherwise, add image to running sum
        if (iteration == 0):
            sum_arr = ofile[1].data
            err_arr = (ofile[2].data) ** 2
            dq_arr = ofile[3].data
            totalweight = ncombine
            totaltime = hdr0['texptime']
        else:
            sum_arr += ofile[1].data
            err_arr += (ofile[2].data) ** 2
            dq_arr = dq_arr | ofile[3].data
            totalweight += ncombine
            totaltime += hdr0['texptime']

    # Then divide by the sum of the weighting factors.
    mean_arr = sum_arr / totalweight
    mean_err_arr = np.sqrt(err_arr / (totalweight ** 2))
    #Update exptime and number of orbits

    hdr0['texptime'] = totaltime
    hdr1['ncombine'] = totalweight
    hdr1['exptime'] = totaltime

    hdu0 = pyfits.PrimaryHDU(header = hdr0)
    hdu1 = pyfits.ImageHDU(mean_arr, header = hdr1)
    hdu2 = pyfits.ImageHDU(mean_err_arr, header =  ofile[2].header)
    hdu3 = pyfits.ImageHDU(dq_arr, header =  ofile[3].header)

    hdulist = pyfits.HDUList([hdu0, hdu1, hdu2, hdu3])
    hdulist.writeto(mean_file)

    return mean_file, totalweight

#---------------------------------------------------------------------------

def calibrate( input_file ):
    """ calibrate input file

    """

    if not 'oref' in os.environ:
        os.environ['oref'] = '/grp/hst/cdbs/oref/'

    print 'Calibrating %s' % (input_file)
    output_blev = input_file.replace('.fits','_blev.fits')
    functions.RemoveIfThere( output_blev )
    output_crj = input_file.replace('.fits','_crj.fits')
    functions.RemoveIfThere( output_crj )

    hdu = pyfits.open( input_file )
    nimset = hdu[0].header['nextend'] / 3
    nrptexp = hdu[0].header['nrptexp']
    crcorr = hdu[0].header['crcorr']
    blevcorr = hdu[0].header['blevcorr']
    hdu.close()
    del hdu

    if (nimset <= 1 and crcorr != "COMPLETE"):
        print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print("Bye now... better luck next time!")
        return None

    if (crcorr != "COMPLETE"):
       
        if (nrptexp != nimset):
            pyfits.setval(input_file, 'NRPTEXP', value=nimset)
            pyfits.setval(input_file, 'CRSPLIT', value=1)

        pyfits.setval(input_file, 'CRCORR', value='PERFORM')
        pyfits.setval(input_file, 'APERTURE', value='50CCD')
        pyfits.setval(input_file, 'APER_FOV', value='50x50')
        if (blevcorr != 'COMPLETE') :
            #print('Performing BLEVCORR')
            pyfits.setval(input_file, 'BLEVCORR', value='PERFORM')
            iraf.basic2d(input_file, output_blev,
                         outblev = '', dqicorr = 'perform', atodcorr = 'omit',
                         blevcorr = 'perform', doppcorr = 'omit', lorscorr = 'omit',
                         glincorr = 'omit', lflgcorr = 'omit', biascorr = 'omit',
                         darkcorr = 'omit', flatcorr = 'omit', shadcorr = 'omit',
                         photcorr = 'omit', statflag = no, verb=no, Stdout='dev$null')
        else:
            #print('Blevcorr alread Performed')
            shutil.copy(input_file, output_blev)

        #print('Performing OCRREJECT')
        iraf.ocrreject(input=output_blev, output=output_crj, verb=no, Stdout='dev$null')

    elif (crcorr == "COMPLETE"):
        print "CR rejection already done"
        os.rename(input_file, output_crj )
  
    pyfits.setval(output_crj, 'FILENAME', value=os.path.split(output_crj)[1] )

    os.remove( output_blev )

    return output_crj

#---------------------------------------------------------------------------

def replace_hot_cols( mean_bias, median_image, residual_image, yfrac=None ):
    """ Replace hot columns in the mean_bias as identified from the 
    residual image with values from the bias_median

    'hot' is 3* sigma

    mean_bias will be updated in place

    """

    print 'Replacing hot column'

    residual_columns_2d = functions.make_resicols_image( residual_image, yfrac=yfrac )
    
    resi_cols_median, resi_cols_mean, resi_cols_std = support.sigma_clip( residual_columns_2d[0], sigma=3, iterations=40 )
    print 'thresh mean,sigma = {} {}'.format( resi_cols_mean, resi_cols_std )
    replval = resi_cols_mean + 3.0 * resi_cols_std
    index = np.where( residual_columns_2d >= replval )

    hdu = pyfits.open( mean_bias, mode='update' )
    hdu[ ('sci', 1) ].data[index] = median_image[index]
    hdu.flush()
    hdu.close()

#---------------------------------------------------------------------------

def replace_hot_pix( mean_bias, median_image ):
    """ Replace image values in residual single hot pixels (defined as those having
    values greater than (mean + 5 sigma of Poisson noise) by those in
    median-filtered bias image. This represents the science extension of the
    final output reference superbias.
    

    mean_bias will be updated in place.

    """
    print 'Replacing hot pixels'
    residual_image = pyfits.getdata( mean_bias, ext=('sci', 1) ) - median_image
    resi_median, resi_mean, resi_std = support.sigma_clip( residual_image )
    fivesig = resi_mean + (5.0 * resi_std)

    index = np.where( residual_image >= fivesig )

    hdu = pyfits.open( mean_bias, mode='update' )
    hdu[ ('sci', 1) ].data[index] = median_image[index]
    hdu.flush()
    hdu.close()

#---------------------------------------------------------------------------

def make_basebias( input_list, refbias_name='basebias.fits' ):
    """ Make the basebias for an anneal month 


    1- Calbrate each bias in the list
    2- Average together the biases
    3- Replace pixels and colums with median values
    4- Set header keywords

    """

    print '#-------------------------------#'
    print '#        Running basejoint      #'
    print '#-------------------------------#'
    print 'Output to %s' % refbias_name

    print 'Processing individual files'
    crj_list = [ calibrate(item) for item in input_list ]
    crj_list = [ item for item in crj_list if item != None ]
 
    mean_bias, totalweight = average_biases( crj_list )

    print 'Replacing hot columns and pixels by median-smoothed values'
    residual_image, median_image = functions.make_residual( mean_bias )

    replace_hot_cols( mean_bias, median_image, residual_image )
    ### for some reason this is done again, but only using the lower 20% of rows
    replace_hot_cols( mean_bias, median_image, residual_image, yfrac=(0, 20) )

    shutil.copy( mean_bias, refbias_name )

    pyfits.setval( refbias_name, 'NCOMBINE', value=totalweight, ext=1)
    functions.update_header_from_input( refbias_name, input_list )

    print 'Cleaning up...'
    functions.RemoveIfThere( mean_bias )
    for item in crj_list:
        functions.RemoveIfThere( item )

    print 'basejoint done'

#------------------------------------------------------------------------------------

