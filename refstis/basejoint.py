"""Procedure to make a monthyl bias for the STIS CCD.

Python translation of a python translation
of an original IRAF cl script to create superbias
reference file from a list of bias frames.

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

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import os
import shutil
import sys
import stistools

from . import functions

#-------------------------------------------------------------------------------

def average_biases(bias_list):
    """Create a weighted sum of the individual input files.

    First make sure all individual input files have been ocrrejected.

    Parameters
    ----------
    bias_list : list
        list of the input biases

    Returns
    -------
    mean_file : str
        name of the averaged filename
    totalweight : float
        sum of the NCOMBINE header keywords from the data

    """

    assert len(bias_list), 'Bias list is empty'

    file_path, file_name = os.path.split(bias_list[0])
    mean_file = os.path.join(file_path, 'mean.fits')

    for iteration, item in enumerate(bias_list):
        with fits.open(item) as hdu:
            hdr0 = hdu[0].header
            hdr1 = hdu[1].header
            nimset = hdr0['nextend'] // 3
            ncombine = hdr1['ncombine']

            #-- If input files have more than one imset or
            #-- have not been cr-rejected, exit
            if (nimset > 1) | (ncombine <= 1):
                print('Input files have to be single imset files and have been CR-rejected')
                print('NIMSET: %d  NCOMBINE: %d' % (nimset, ncombine))
                sys.exit(3)

            #Otherwise, add image to running sum
            if (iteration == 0):
                sum_arr = hdu[1].data
                err_arr = (hdu[2].data) ** 2
                dq_arr = hdu[3].data
                totalweight = ncombine
                totaltime = hdr0['texptime']
            else:
                sum_arr += hdu[1].data
                err_arr += (hdu[2].data) ** 2
                dq_arr |= hdu[3].data
                totalweight += ncombine
                totaltime += hdr0['texptime']

    # Then divide by the sum of the weighting factors.
    mean_arr = sum_arr / float(totalweight)
    mean_err_arr = np.sqrt(err_arr / (totalweight ** 2))
    #Update exptime and number of orbits

    hdr0['texptime'] = totaltime
    hdr1['ncombine'] = totalweight
    hdr1['exptime'] = totaltime

    out_hdu0 = fits.PrimaryHDU(header=hdr0)
    out_hdu1 = fits.ImageHDU(mean_arr, header=hdr1)
    out_hdu2 = fits.ImageHDU(mean_err_arr, header=hdu[2].header)
    out_hdu3 = fits.ImageHDU(dq_arr, header=hdu[3].header)

    hdulist = fits.HDUList([out_hdu0, out_hdu1, out_hdu2, out_hdu3])
    hdulist.writeto(mean_file, output_verify='exception')

    return mean_file, totalweight

#-------------------------------------------------------------------------------

def calibrate(input_file):
    """ calibrate input file

    """

    if not 'oref' in os.environ:
        os.environ['oref'] = '/grp/hst/cdbs/oref/'

    print('Calibrating %s' % (input_file))
    output_blev = input_file.replace('.fits', '_blev.fits')
    functions.RemoveIfThere(output_blev)
    output_crj = input_file.replace('.fits', '_crj.fits')
    functions.RemoveIfThere(output_crj)

    with fits.open(input_file, mode='update') as hdu:
        nimset = hdu[0].header['nextend'] / 3
        nrptexp = hdu[0].header['nrptexp']
        crcorr = hdu[0].header['crcorr']
        blevcorr = hdu[0].header['blevcorr']

        if (nimset <= 1 and crcorr != "COMPLETE"):
            print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
            print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
            print("Bye now... better luck next time!")
            return None

        if (crcorr != "COMPLETE"):

            if (nrptexp != nimset):
                hdu[0].header['NRPTEXP'] = nimset
                hdu[0].header['CRSPLIT'] = 1

            hdu[0].header['CRCORR'] = 'PERFORM'
            hdu[0].header['APERTURE'] = '50CCD'
            hdu[0].header['APER_FOV'] = '50x50'

            if (blevcorr != 'COMPLETE'):
                #print('Performing BLEVCORR')
                hdu[0].header['BLEVCORR'] = 'PERFORM'
                stistools.basic2d.basic2d(input=input_file,
                                          output=output_blev,
                                          outblev='',
                                          dqicorr='perform',
                                          atodcorr='omit',
                                          blevcorr='perform',
                                          doppcorr='omit',
                                          lorscorr='omit',
                                          glincorr='omit',
                                          lflgcorr='omit',
                                          biascorr='omit',
                                          darkcorr='omit',
                                          flatcorr='omit',
                                          shadcorr='omit',
                                          photcorr='omit',
                                          statflag=False,
                                          verbose=False,
                                          trailer="/dev/null")
            else:
                #print('Blevcorr alread Performed')
                shutil.copy(input_file, output_blev)

            #print('Performing OCRREJECT')
            stistools.ocrreject.ocrreject(input=output_blev,
                                          output=output_crj,
                                          verbose=False,
                                          trailer="/dev/null")

        elif (crcorr == "COMPLETE"):
            print("CR rejection already done")
            os.rename(input_file, output_crj)

    fits.setval(output_crj, 'FILENAME', value=os.path.split(output_crj)[-1])

    os.remove(output_blev)

    return output_crj

#-------------------------------------------------------------------------------

def replace_hot_cols(mean_bias, median_image, residual_image, yfrac=1):
    """ Replace hot columns in the mean_bias as identified from the
    residual image with values from the bias_median

    'hot' is 3* sigma

    mean_bias will be updated in place

    """

    print('Replacing hot column')
    residual_columns_2d = functions.make_resicols_image(residual_image,
                                                        yfrac=yfrac)

    resi_cols_mean, resi_cols_median, resi_cols_std = sigma_clipped_stats(residual_columns_2d[0], sigma=3, maxiters=40)

    print('thresh mean,sigma = {} {}'.format(resi_cols_mean, resi_cols_std))
    replval = resi_cols_mean + 3.0 * resi_cols_std
    index = np.where(residual_columns_2d >= replval)

    with fits.open(mean_bias, mode='update') as hdu:
        hdu[('sci', 1)].data[index] = median_image[index]

#-------------------------------------------------------------------------------

def replace_hot_pix(mean_bias, median_image):
    """ Replace image values in residual single hot pixels

    defined as those having
    values greater than (mean + 5 sigma of Poisson noise) by those in
    median-filtered bias image. This represents the science extension of the
    final output reference superbias.

    mean_bias will be updated in place.

    Parameters
    ----------
    mean_bias : str
        name of the mean bias bias
    median_image : np.ndarray
        2d median image of the bias

    """

    print('Replacing hot pixels')
    residual_image = fits.getdata(mean_bias, ext=('sci', 1)) - median_image
    resi_mean, resi_median, resi_std = sigma_clipped_stats(residual_image,
                                                           sigma=5,
                                                           maxiters=40)

    fivesig = resi_mean + (5.0 * resi_std)
    print("  hot is > {}".format(fivesig))
    index = np.where(residual_image >= fivesig)

    with fits.open(mean_bias, mode='update') as hdu:
        hdu[('sci', 1)].data[index] = median_image[index]

#-------------------------------------------------------------------------------

def make_basebias(input_list, refbias_name='basebias.fits'):
    """ Make the basebias for an anneal month


    1- Calbrate each bias in the list
    2- Average together the biases
    3- Replace pixels and colums with median values
    4- Set header keywords

    Parameters
    ----------
    input_list : list
        list of input bias files.
    refbias_name : str
        name of the output reference file.

    """

    print('#-------------------------------#')
    print('#        Running basejoint      #')
    print('#-------------------------------#')
    print('Output to %s' % refbias_name)

    print('Processing individual files')
    crj_list = [calibrate(item) for item in input_list]
    crj_list = [item for item in crj_list if item != None]

    mean_bias, totalweight = average_biases(crj_list)

    print('Replacing hot columns and pixels by median-smoothed values')
    residual_image, median_image = functions.make_residual(mean_bias)

    replace_hot_cols(mean_bias, median_image, residual_image)
    #-- then again, but only using the lower 20% of rows
    replace_hot_cols(mean_bias, median_image, residual_image, yfrac=.2)

    replace_hot_pix(mean_bias, median_image)

    shutil.copy(mean_bias, refbias_name)

    functions.update_header_from_input(refbias_name, input_list)
    fits.setval(refbias_name, 'NCOMBINE', value=totalweight, ext=1)
    fits.setval(refbias_name, 'TASKNAME', ext=0, value='BASEJOIN')

    print('Cleaning up...')
    functions.RemoveIfThere(mean_bias)
    for item in crj_list:
        functions.RemoveIfThere(item)

    print('basejoint done for {}'.format(refbias_name))

#-------------------------------------------------------------------------------
