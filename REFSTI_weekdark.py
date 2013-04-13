from pyraf import iraf
from pyraf.irafglobals import *
import os
import sys
import glob
import shutil
import tempfile
import pyfits
import stiref
import opusutil

nm = os.path.basename(sys.argv[0])

# bound "method" expected to be frequently used
MakeNameSafe = stiref.makeNameSafe

#---------------------------------------------------------------------------
def make_weekdark(infiles, basedark, biasfile, weekoutfile, outfile, reffile):
  do_cal = yes            # Perform bias subtraction and cr-reject?
  maxiter = 40            # Maximum number of iterations for imstat'))
  lower = INDEF           # Initial lower limit for imstat'))
  upper = INDEF           # Initial upper limit for imstat'))
  verbose = 0             # Show results of imstat iterations?'))

  #
  #*********************************************************************
  # Load necessary packages
  #
  iraf.stsdas()
  iraf.imgtools()
  iraf.ttools()
  iraf.mstools()
  iraf.hst_calib()
  iraf.stis()
  
  #
  #*********************************************************************
  # Expand input file list into temporary files
  #
  (dummy, tmp) = tempfile.mkstemp(suffix='_weekdark', dir=workdir)
  inlist = MakeNameSafe(tmp)
  iraf.sections(c_infiles, option = 'fullname', Stdout=inlist)
  nfiles = int(iraf.sections.nimages)
  imglist = inlist

  print 'Splitting images'
  imset_count = functions.split_images( imglist ) 
  
  print 'Joining images'
  msjoin_list = ','.join( [ item for item in glob.glob('*raw??.fits') ] )# if item[:9] in bias_list] )
  print msjoin_list
  joined_out = refdark_name+ '_joined' +'.fits' 
  print joined_out
  
  msjoin_file = open('msjoin.txt','w')
  msjoin_file.write( '\n'.join(msjoin_list.split(',')) )
  msjoin_file.close()
  
  iraf.msjoin( inimg='@msjoin.txt', outimg=joined_out, Stderr='dev$null')

  # test for the need to perform cosmic-ray rejection, and do it
  crdone = stiref.bd_crreject(joinedfile)
  print "## crdone is ", crdone
  if (not crdone):
      functions.bd_calstis(joinedfile, thebiasfile)

  # divide cr-rejected
  crj_filename = refdark_name + '_crj.fits'
  exptime = pyfits.getval( crj_filename, 'TEXPTIME', ext=0 )
  gain = pyfits.getval( crj_filename, 'ATODGAIN', ext=0 )
  xbin = pyfits.getval( crj_filename, 'XBIN', ext=0 )
  ybin = pyfits.getval( crj_filename, 'YBIN', ext=0 )
  
  normalize_factor = float(exptime)/gain # ensure floating point
  
  norm_filename = crj_filename.replace('_crj.fits','_norm.fits')
  iraf.msarith( crj_filename, '/', normalize_factor, norm_filename ,verbose=0)  

  pyfits.setval( norm_filename, 'TEXPTIME', value=1 )


  #
  #***************************************************************************
  # Perform iterative statistics on the normalized superdark
  #   (i.e., neglecting hot pixels in the process)
  #
  iter_count,median,sigma,npx,med,mod,min,max = functions.iterate( norm_filename )
  five_sigma = median + 5*sigma

  # save hot pixel level and the name of the baseline dark
  # for use in updating history
  out_fd = open(weekoutfile, 'w')
  out_fd.write(str(p_fivesig)+' '+basedark)
  out_fd.close()

  #
  #***************************************************************************
  # Perform iterative statistics on the baseline superdark
  #   (i.e., neglecting hot pixels in the process)
  #
  print "## Perform iterative statistics on the baseline superdark (thebasedark)"
  iter_count,base_median,base_sigma,npx,med,mod,min,max = functions.iterate( norm_filename )
  five_sigma = base_median + 5*base_sigma

  print "## Create zeodark", zerodark, "by subtct'g basemedian", basemed
  print "##           from", theoutfile+"_sci.fits"
  iraf.imarith(theoutfile+"_sci.fits", "-", basemed, zerodark)
  print "## imcalc theoutfile(+?):", theoutfile+"_sci.fits[0],"+zerodark+"[0]", ", only_hotpix:", only_hotpix
  iraf.imcalc (theoutfile+"_sci.fits[0],"+zerodark+"[0]", only_hotpix,
        "if im1 .ge. "+str(p_fivesig)+" then im2 else 0.0",
        pixtype="old", nullval=0., verbose=no)
  #
  #***************************************************************************
  # Create median-filtered version of super-de-buper dark
  #
  print "## Create median-filtered version of super-de-buper dark "
  print "## run iraf.median on thebasedark", thebasedark+".fits[sci]", "using:", basedrk_med
  iraf.median (thebasedark+".fits[sci]", basedrk_med, xwindow=5,
               ywindow=5, verbose=no)

  #
  #***************************************************************************
  # Create "only baseline dark current" image from these two
  #
  print "## Create 'only baseline dark current' image from these two "
  iraf.imcalc (thebasedark+"_sci.fits[0],"+basedrk_med+"[0]",only_dark,
        "if im1 .ge. "+str(p_fivesig)+" then im2 else im1",
        pixtype="old", nullval=0., verbose=no)
  #
  #***************************************************************************
  # Add "only baseline dark current" image to "only hot pixels" image.
  # This creates the science portion of the forthcoming reference dark.
  #
  print "## Add 'only baseline dark current' image to 'only hot pixels' image. "
  print "## This creates the science portion of the forthcoming reference dark. "
  print "## imarith(only_dark", only_dark+"[0]"
  print "##     + only_hotpix", only_hotpix+"[0]"
  print "##     =   superdark", superdark, ")"
  iraf.imarith(only_dark+"[0]","+",only_hotpix+"[0]",superdark)

  #
  #***************************************************************************
  # Update DQ extension of normalized dark by assigning the value 16 to the
  # hot pixels, and put the result in temporary image.
  #
  print "## Use imcalc to update DQ  extension of normalized dark"
  print "## imcalc input theoutfile", theoutfile+"_dq.fits[0],"+only_hotpix+"[0]", "by equation to output", refDQ
  iraf.imcalc (theoutfile+"_dq.fits[0],"+only_hotpix+"[0]", refDQ,
        "if im2 .ge. "+str(p_fivesig)+" then 16 else im1", 
        pixtype="old", nullval=0., verbose=no)
  #***************************************************************************
  # Update ERR extension of new superdark by assigning the ERR values of the
  # basedark except for the new hot pixels that are updated from the weekly
  # superdark, for which the error extension of the weekly superdark is taken.
  # Put the result in temporary ERR image.
  #
  print "## Use imcalc to update ERR extension of new superdark "
  print "## imcalc input thebasedark", thebasedark+"_err.fits[0],"+only_hotpix+"[0],"+theoutfile+"_err.fits[0]", 
  print "## by equation to output ERR_new", ERR_new
  iraf.imcalc(thebasedark+"_err.fits[0],"+only_hotpix+"[0],"+theoutfile+"_err.fits[0]",
        ERR_new, "if im2 .eq. 0.0 then im1 else im3", 
        pixtype="old", nullval=0., verbose=no)


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

  #
  #***************************************************************************
  # Copy new superdark into image extension of reference dark; copy ERR
  # extension from normalized weekly superdark, and DQ extension from
  # temporary file created previously.
  #
  print "## Copy new superdark into image extension of reference dark; "
  print "## copy ERR extension from normalized weekly superdark, "
  print "## and DQ extension from temporary file created previously. "
  print "## iraf.imcopy (superdark", superdark,    "to thereffile", thereffile+".fits[sci,1][*,*]", ", verbose=no)"
  print "## iraf.imcopy (ERR_new",  ERR_new+"[0]"
  print "##        to thereffile", thereffile+".fits[err,1][*,*]", ", verbose=no)"
  print "## iraf.imcopy (refDQ",    refDQ+"[0]",   
  print "##      to thereffile", thereffile+".fits[dq,1][*,*]",  ", verbose=no)"
  iraf.imcopy (superdark, thereffile+".fits[sci,1][*,*]", verbose=no)
  iraf.imcopy (ERR_new+"[0]", thereffile+".fits[err,1][*,*]", verbose=no)
  iraf.imcopy (refDQ+"[0]",  thereffile+".fits[dq,1][*,*]", verbose=no)



if __name__ == "__main__":
    
