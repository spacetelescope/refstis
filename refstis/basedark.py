"""

Functions to produce a monthly basedark for the STIS Darks and Biases
reference file pipeline

"""

from astropy.io import fits
import numpy as np
import shutil

import functions
import support

#--------------------------------------------------------------------------

def find_hotpix(filename):
    """ find pixels hotter that median + 5*sigma and update DQ array with 
    DQ == 16

    Updates file in place.
    
    Parameters:
    -----------
    filename, str
        filename of the input biasfile

    """

    with fits.open( filename, mode='update' ) as hdu:
        data_median, data_mean, data_std = support.sigma_clip(hdu[ ('sci', 1) ].data, 
                                                              sigma=3, 
                                                              iterations=40)

        five_sigma = data_median + 5 * data_std
        index = np.where((hdu[ ('SCI', 1) ].data > five_sigma) & 
                         (hdu[('SCI', 1)].data > data_mean + 0.1))

        hdu[('DQ', 1)].data[index] = 16

#--------------------------------------------------------------------------

def make_basedark(input_list, refdark_name='basedark.fits', bias_file=None):
    """Make a monthly baseline dark from the input list

    1 - If not already done, perform bias subtraction (I believe this should include blevcorr and dqicorr)
    2 - If after switch to side-2 electronics, perform temperature scaling
    3- Join all imsets from input list into single file
    4- combine and cr-reject
    5- normalize to e/s by dividing by (exptime/gain)
    6- update DQ array with hot pixel information

    Parameters:
    -----------
    input_list : list
        list of input dark files
    refdark_name : str
        name of the output reference file
    bias_file : str
        bias file to be used in calibration

    """

    print '#-------------------------------#'
    print '#        Running basedark       #'
    print '#-------------------------------#'
    print 'output to: %s' % refdark_name
    print 'with biasfile %s' % bias_file

    #-- bias subtract data if not already done
    for i, filename in enumerate(input_list):
        filename = functions.bias_subtract_data(filename)
        input_list[i] = filename
        #Side 1 operations ended on May 16, 2001. 
        #Side 2 operations started on July 10, 2001,
        #52091.0 corresponds to July 1, 2001
        texpstrt = fits.getval(filename, 'texpstrt', 0)
        if texpstrt > 52091.0:
            functions.apply_dark_correction(filename, texpstrt)

    joined_filename = refdark_name.replace('.fits', '_joined.fits') 
    crj_filename = joined_filename.replace('.fits', '_crj.fits')

    #if not bias_file: raise IOError('No biasfile specified, this task needs one to run')

    print 'Joining images'
    functions.msjoin(input_list, joined_filename)

    print 'CRREJECT'
    crdone = functions.bd_crreject(joined_filename)
    if not crdone:
        functions.bd_calstis(joined_filename, bias_file)

    functions.normalize_crj(crj_filename)

    shutil.copy(crj_filename, refdark_name)

    find_hotpix(refdark_name)

    functions.update_header_from_input(refdark_name, input_list)

    print 'Cleaning...'
    functions.RemoveIfThere(crj_filename)
    functions.RemoveIfThere(joined_filename)

    print 'basedark done for {}'.format(refdark_name)
    
#--------------------------------------------------------------------------

    
    
