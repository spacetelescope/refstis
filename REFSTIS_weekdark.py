from pyraf import iraf
from pyraf.irafglobals import *
import os
import sys
import glob
import shutil
import pyfits
from scipy.signal import medfilt
import REFSTIS_functions
import numpy as np

#-------------------------------------------------------------------------------

def split_and_join( input_list, refname ):
  # not done with this...tired.
  '''Split imsets and join them back together. 
  
  return the joined file name
  '''
  refname = refname.replace('.fits','')
  joined_out = refname + '_joined.fits'

  rootname_list = [ os.path.split(item)[1][:9] for item in input_list ]
  imset_count = REFSTIS_functions.split_images( input_list ) 
  
  print 'Joining images'
  msjoin_list = ','.join( [ item for item in glob.glob('*raw??.fits') if item[:9] in rootname_list] )

  print msjoin_list
  msjoin_file = open('msjoin_%s.txt'%(refname),'w')
  msjoin_file.write( '\n'.join(msjoin_list.split(',')) )
  msjoin_file.close()
  
  iraf.msjoin( inimg='@msjoin.txt', outimg=joined_out, Stderr='dev$null')

  return joined_out,imset_count

#-------------------------------------------------------------------------------

def make_weekdark( input_list, refdark_name, thebiasfile, thebasedark ):
    """
    1- split all raw images into their imsets
    2- join imsets together into a single file
    3- combine and cr-reject
    4- normalize to e/s by dividing by (exptime/gain)
    5- do hot pixel things

    """
    iraf.stsdas()
    iraf.imgtools()
    iraf.ttools()
    iraf.mstools()
    iraf.hst_calib()
    iraf.stis()


    print '#-------------------------------#'
    print '#        Running weekdark       #'
    print '#-------------------------------#'
    print 'Making weekdark %s' % (refdark_name)
    print 'With : %s' % (thebiasfile)
    print '     : %s' % (thebasedark)

    refdark_path = os.path.split( refdark_name )[0] or './'

    rootname_set = set( [ os.path.split(item)[1][:9] for item in input_list] )

    n_imsets = REFSTIS_functions.split_images( input_list ) 

    joined_out = refdark_name.replace('.fits', '_joined.fits' )
    print 'Joining images to %s' % joined_out
    msjoin_list = [ item for item in 
                    glob.glob( os.path.join(refdark_path,'*raw??.fits') )  if os.path.split(item)[1][:9] in rootname_set]
    REFSTIS_functions.msjoin( msjoin_list, joined_out)


    # test for the need to perform cosmic-ray rejection, and do it
    crdone = REFSTIS_functions.bd_crreject( joined_out )
    print "## crdone is ", crdone
    if (not crdone):
        REFSTIS_functions.bd_calstis(joined_out, thebiasfile)

    # divide cr-rejected
    crj_filename = joined_out.replace('.fits','_crj.fits')
    exptime = pyfits.getval( crj_filename, 'TEXPTIME', ext=0 )
    gain = pyfits.getval( crj_filename, 'ATODGAIN', ext=0 )
    xbin = pyfits.getval( crj_filename, 'BINAXIS1', ext=0 )
    ybin = pyfits.getval( crj_filename, 'BINAXIS2', ext=0 )

    normalize_factor = float(exptime) / gain # ensure floating point

    norm_filename = crj_filename.replace('_crj.fits','_norm.fits')
    iraf.msarith( crj_filename, '/', normalize_factor, norm_filename ,verbose=0)  

    pyfits.setval( norm_filename, 'TEXPTIME', value=1 )


    #
    #***************************************************************************
    # Perform iterative statistics on the normalized superdark
    #   (i.e., neglecting hot pixels in the process)
    #
    iter_count,median,sigma,npx,med,mod,min,max = REFSTIS_functions.iterate( norm_filename )
    five_sigma = median + 5*sigma

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
    iter_count,base_median,base_sigma,npx,basemed,mod,min,max = REFSTIS_functions.iterate( norm_filename )
    five_sigma = base_median + 5*base_sigma

    print "## Create median-filtered version of super-de-buper dark "
    basedark_hdu = pyfits.open( thebasedark )
    basedark_med = medfilt( basedark_hdu[('sci',1)].data, (5, 5) )

    norm_hdu = pyfits.open( norm_filename )

    zerodark = norm_hdu[ ('sci',1) ].data - basemed
    only_hotpix = np.where( basedark_hdu[ ('sci',1) ].data >= five_sigma, zerodark, 0 )

    print "## Create 'only baseline dark current' image from these two "
    only_dark = np.where( basedark_hdu[ ('sci',1) ].data >= five_sigma, basedark_med, basedark_hdu[ ('sci',1) ].data )

    print "## Add 'only baseline dark current' image to 'only hot pixels' image. "
    print "## This creates the science portion of the forthcoming reference dark. "
    superdark = only_dark + only_hotpix
    norm_hdu[ ('sci',1) ].data = superdark

    #
    #***************************************************************************
    # Update DQ extension of normalized dark by assigning the value 16 to the
    # hot pixels, and put the result in temporary image.
    #
    print "## Use imcalc to update DQ  extension of normalized dark"
    index = np.where( only_hotpix >= five_sigma )[0]
    norm_hdu[ ('dq',1) ].data[index] = 16

    #***************************************************************************
    # Update ERR extension of new superdark by assigning the ERR values of the
    # basedark except for the new hot pixels that are updated from the weekly
    # superdark, for which the error extension of the weekly superdark is taken.
    # Put the result in temporary ERR image.
    #
    print "## Use imcalc to update ERR extension of new superdark "
    ### This is obviously wrong
    index = np.where( only_hotpix == 0 )
    norm_hdu[ ('err',1) ].data[index] = basedark_hdu[ ('err', 1) ].data[index]
    index = np.where( only_hotpix != 0 )
    norm_hdu[ ('err',1) ].data[index] = basedark_hdu[ ('err', 1) ].data[index]


    pyfits.setval(norm_filename, "FILENAME", value=norm_filename)
    pyfits.setval(norm_filename, "FILETYPE", value="DARK IMAGE")
    pyfits.setval(norm_filename, "CCDAMP", value="ANY")
    pyfits.setval(norm_filename, "CCDGAIN", value="-1")
    pyfits.setval(norm_filename, "BINAXIS1", value=xbin)
    pyfits.setval(norm_filename, "BINAXIS2", value=ybin)
    pyfits.setval(norm_filename, "USEAFTER", value=" ")
    pyfits.setval(norm_filename, "PEDIGREE", value="INFLIGHT")
    pyfits.setval(norm_filename, "DESCRIP",
                  value="Weekly superdark created by J. Ely")
    pyfits.setval(norm_filename, "NEXTEND", value=3)
    pyfits.setval(norm_filename, "COMMENT", 
                  value="created by the STIS weekdark task in the reference file pipeline")

    shutil.copy( norm_filename, refdark_name )

    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( norm_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )
    for item in msjoin_list:
        REFSTIS_functions.RemoveIfThere( item ) 

    print '#-------------------------------#'
    print '#        Finished weekdark      #'
    print '#-------------------------------#'

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    make_weekdark( glob.glob(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4] )
