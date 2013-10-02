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

    rootname_set = set( [item[:9] for item in input_list] )

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
    basedark_med = medfilt( pyfits.getdata(thebasedark,ext=('sci',1)),(5,5) )

    norm_hdu = pyfits.open( norm_filename )

    zerodark = norm_hdu[ ('sci',1) ].data - basemed
    only_hot = np.where( theoutfile >= five_sigma, zerodark, 0 )

    print "## Create 'only baseline dark current' image from these two "
    only_dark = np.where( thebasedark >= five_sigma, basedark_med, thebasedark )

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
    index = np.where( only_hotpix >= p_fivesig )[0]
    norm_hdu[ ('dq',1) ].data[index] = 16

    #***************************************************************************
    # Update ERR extension of new superdark by assigning the ERR values of the
    # basedark except for the new hot pixels that are updated from the weekly
    # superdark, for which the error extension of the weekly superdark is taken.
    # Put the result in temporary ERR image.
    #
    print "## Use imcalc to update ERR extension of new superdark "
    index = np.where( only_hotpix == 0 )[0]
    norm_hdu[ ('err',1) ].data[index] = thebasedark_error[index]
    index = np.where( only_hotpix != 0 )[0]
    norm_hdu[ ('err',1) ].data[index] = theoutfile_error[index]


    pyfits.setval(thereffile+".fits[0]", "FILENAME", thereffile)
    pyfits.setval(thereffile+".fits[0]", "FILETYPE", "DARK IMAGE")
    pyfits.setval(thereffile+".fits[0]", "CCDAMP", "ANY")
    pyfits.setval(thereffile+".fits[0]", "CCDGAIN", "-1")
    pyfits.setval(thereffile+".fits[0]", "BINAXIS1", xbin)
    pyfits.setval(thereffile+".fits[0]", "BINAXIS2", ybin)
    pyfits.setval(thereffile+".fits[0]", "USEAFTER", " ")
    pyfits.setval(thereffile+".fits[0]", "PEDIGREE", "INFLIGHT")
    pyfits.setval(thereffile+".fits[0]", "DESCRIP",
                   "Weekly superdark created by J. Ely")
    pyfits.setval(thereffile+".fits[0]", "NEXTEND", 3)
    pyfits.setval(thereffile+".fits[0]", "COMMENT", 
                  "created by the STIS weekdark task in the reference file pipeline")

    print '#-------------------------------#'
    print '#        Finished weekdark      #'
    print '#-------------------------------#'

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    make_weekdark( glob.glob(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4] )
