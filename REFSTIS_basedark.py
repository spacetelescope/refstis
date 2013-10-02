import REFSTIS_functions
import numpy as np
import pyfits
import glob
import sys
import shutil

#--------------------------------------------------------------------------

def hotpix( filename ):
    pass

#--------------------------------------------------------------------------

def update_header( filename,xbin,ybin ):
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

def make_basedark( input_list, refdark_name='basedark.fits', bias_file=None, oref=None):
    """
    1- split all raw images into their imsets
    2- join imsets together into a single file
    3- combine and cr-reject
    4- normalize to e/s by dividing by (exptime/gain)
    5- do hot pixel things

    """

    print '#-------------------------------#'
    print '#        Running basedark       #'
    print '#-------------------------------#'

    from pyraf import iraf
    from iraf import stsdas,toolbox,imgtools,mstools
    import os
    import shutil
    import REFSTIS_functions

    if not bias_file:
        raise IOError( 'No biasfile specified, this task needs one to run' )


    if (not oref) and (not 'oref' in os.environ):
        oref = '/grp/hst/cdbs/oref/'

    rootname_set = set( [item[:9] for item in input_list] )

    refdark_name = refdark_name
    refdark_path = os.path.split( refdark_name )[0] or './'

    imset_count = REFSTIS_functions.split_images( input_list )

    print 'Joining images'
    msjoin_list = [ item for item in 
                    glob.glob( os.path.join(refdark_path,'*raw??.fits') )  if os.path.split(item)[1][:9] in rootname_set]
    joined_out = refdark_name.replace('.fits', '_joined.fits') 

    REFSTIS_functions.msjoin( msjoin_list, joined_out  )

    # ocrreject
    print 'CRREJECT'
    iraf.chdir( refdark_path )
    crdone = REFSTIS_functions.bd_crreject( joined_out )
    if (not crdone):
        REFSTIS_functions.bd_calstis( joined_out, bias_file )

    # divide cr-rejected
    print 'Dividing'
    crj_filename = joined_out.replace('.fits', '_crj.fits')
    exptime = pyfits.getval( crj_filename, 'TEXPTIME', ext=0 )
    gain = pyfits.getval( crj_filename, 'ATODGAIN', ext=0 )
    xbin = pyfits.getval( crj_filename, 'BINAXIS1', ext=0 )
    ybin = pyfits.getval( crj_filename, 'BINAXIS2', ext=0 )
    
    normalize_factor = float(exptime) / gain # ensure floating point

    norm_filename = crj_filename.replace('_crj.fits','_norm.fits')
    iraf.msarith( crj_filename, '/', normalize_factor, norm_filename ,verbose=0)  

    pyfits.setval( norm_filename, 'TEXPTIME', value=1 )

    # hotpixel stuff
    iter_count,median,sigma,npx,med,mod,min,max = REFSTIS_functions.iterate( norm_filename )
    five_sigma = median + 5*sigma
    
    shutil.copy( norm_filename, refdark_name )

    norm_hdu = pyfits.open( norm_filename,mode='update' )
    index = np.where( norm_hdu[ ('SCI',1) ].data >= five_sigma + .1)[0]
    norm_hdu[ ('DQ',1) ].data[index] = 16
    norm_hdu.flush()
    norm_hdu.close()

    print 'Cleaning'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( norm_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )
    for item in msjoin_list:
        REFSTIS_functions.RemoveIfThere( item )

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

if __name__ == "__main__":
    make_basedark( glob.glob(sys.argv[1]), sys.argv[2], sys.argv[3] )

    
    
