"""
Functions to create weekly superdarks for the STIS darks and bias reference
file pipeline.

"""

import shutil
from astropy.io import fits as pyfits
from scipy.signal import medfilt
import numpy as np

import support
import functions

#-------------------------------------------------------------------------------

def create_superdark(crj_filename, basedark):
    """

    files will be updated in place

    # Add 'only baseline dark current' image to 'only hot pixels' image. 
    # This creates the science portion of the forthcoming reference dark. 

    """

    crj_hdu = pyfits.open( crj_filename, mode='update' )

    ## Perform iterative statistics on this normalized superdark 
    #imstat(img,fields="mean,stddev,npix,midpt,mode,min,max",lower=ll,upper=ul,for-) | scan(mn,sig,nx,med,mod,min,max)
    data_median, data_mean, data_std = support.sigma_clip( crj_hdu[ ('sci', 1) ].data , sigma=3, iterations=40 )

    #	p_fivesig = med + (5*sig)
    p_five_sigma = data_median + (5*data_std)
    print 'hot pixels are defined as above: ', p_five_sigma
    basedark_hdu = pyfits.open( basedark )
    base_median, base_mean, base_std = support.sigma_clip( basedark_hdu[ ('sci', 1) ].data, sigma=3, iterations=40 )

    zerodark = crj_hdu[ ('sci', 1) ].data - base_median
    only_hotpix = np.where( crj_hdu[ ('sci', 1) ].data >= p_five_sigma, zerodark, 0.0 )

    basedark_med = medfilt( basedark_hdu[('sci', 1)].data, (5, 5) )
    only_dark = np.where( basedark_hdu[ ('sci', 1) ].data >= p_five_sigma, basedark_med, basedark_hdu[ ('sci', 1) ].data )

    crj_hdu[ ('sci', 1) ].data = only_dark + only_hotpix

    #divide error by the sqrt of number of combined images
    #imcalc (s3//"[err]", "refERR.fits", "im1/sqrt("//ncombine//")", \
    #	pixtype="old", nullval=0., verbose-)
    #!!!!Very early data may not have ncombine
    #This doesn't appear to be used in the final xstis analysis
    #crj_hdu[('err', 1)].data /= np.sqrt(crj_hdu[1].header['ncombine'])

    #- update DQ extension 
    crj_hdu[ ('dq', 1) ].data = np.where( only_hotpix >= p_five_sigma, 16, crj_hdu[('dq', 1)].data)

    #- Update Error
    #imcalc (thebasedark//"_err.fits[0],only_hotpix.fits[0],"//theoutfile//"_err.fits[0]", \
	#"ERR_new.fits", "if im2 .eq. 0.0 then im1 else im3", \
	#pixtype="old", nullval=0., verbose-)
    

    crj_hdu[ ('err', 1) ].data = np.where(only_hotpix == 0, basedark_hdu[ ('err', 1) ].data, crj_hdu[('err', 1)].data)
    pdb.set_trace()
    crj_hdu.flush()
    crj_hdu.close()

#-------------------------------------------------------------------------------

def make_weekdark(input_list, refdark_name, thebasedark, thebiasfile = None):
    """
    1- If not already done, run basic2d with blevcorr, biascorr, and dqicorr 
        set to perform 
    2- Apply temperature correction to the data
    3- split all raw images into their imsets
    4- join imsets together into a single file
    5- combine and cr-reject
    6- normalize to e/s by dividing by (exptime/gain)
    7- do hot pixel things

    # Update ERR extension of new superdark by assigning the ERR values of the
    # basedark except for the new hot pixels that are updated from the weekly
    # superdark, for which the error extension of the weekly superdark is taken.

    """

    print '#-------------------------------#'
    print '#        Running weekdark       #'
    print '#-------------------------------#'
    if not thebiasfile:
        thebiasfile = pyfits.getval(input_list[0], 'biasfile', 0)

    print 'Making weekdark %s' % (refdark_name)
    print 'With : %s' % (thebiasfile)
    print '     : %s' % (thebasedark)

    joined_out = refdark_name.replace('.fits', '_joined.fits' )
    for i, filename in enumerate(input_list):
        filename = functions.bias_subtract_data(filename)
        input_list[i] = filename
        #Side 1 operations ended on May 16, 2001. 
        #Side 2 operations started on July 10, 2001, 
        #52091.0 corresponds to July 1, 2001
        if pyfits.getval(filename, 'texpstrt', 0) > 52091.0:
            functions.apply_dark_correction(filename, pyfits.getval(filename, 'texpstrt', 0))
    
    print 'Joining images to %s' % joined_out
    functions.msjoin( input_list, joined_out)

    crdone = functions.bd_crreject( joined_out )
    print "## crdone is ", crdone
    if (not crdone):
        functions.bd_calstis(joined_out, thebiasfile)
    crj_filename = joined_out.replace('.fits', '_crj.fits')
    shutil.copy( crj_filename, refdark_name )
    functions.normalize_crj( refdark_name )

    create_superdark( refdark_name, thebasedark )

    functions.update_header_from_input( refdark_name, input_list )

    print 'Cleaning up...'
    functions.RemoveIfThere( crj_filename )
    functions.RemoveIfThere( joined_out )

    print 'Weekdark done'

#-------------------------------------------------------------------------------