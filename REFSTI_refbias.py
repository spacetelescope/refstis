from pyraf import iraf
from pyraf.irafglobals import *
from REFSTI_functions import split_images,crreject
import pyfits
import sys
import glob
import shutil
import os

def update_header( filename ):
    pass

def make_refbias( bias_list, refbias_name ):
    """ Basic premise for making a refbias:
    1- split all raw images into their raw imsets
    2- join imsets together into one large file
    3- combine and cosmic ray screen joined imset
    4- divide by number of imsets
    5- set header keywords
    """
    from pyraf import iraf
    from iraf import stsdas,toolbox,imgtools,mstools
    import os

    refbias_name.replace('.fits','') # ensure its just a rootname

    print 'Splitting images'
    split_images( bias_list )
    
    print 'Joining images'
    msjoin_list = ','.join( [ item for item in glob.glob('*raw??.fits') ] )# if item[:9] in bias_list] )
    print msjoin_list
    joined_out = refbias_name+ '_joined' +'.fits' 
    print joined_out
    
    msjoin_file = open('msjoin.txt','w')
    msjoin_file.write( '\n'.join(msjoin_list.split(',')) )
    msjoin_file.close()
    
    iraf.msjoin( inimg='@msjoin.txt', outimg=joined_out, Stderr='dev$null')
    
    crj_filename = crreject( joined_out )
    out_name = refbias_name + '.fits' 
    shutil.copy( crj_filename, out_name)

    pyfits.setval( out_name, 'FILENAME', value=out_name )
    pyfits.setval( out_name, 'FILETYPE', value='CCD BIAS IMAGE')
    pyfits.setval( out_name, 'CCDGAIN', value= 1 ) ### FIX
    pyfits.setval( out_name, 'BINAXIS1', value= 1 ) ### FIX
    pyfits.setval( out_name, 'BINAXIS2', value= 1 ) ### FIX
    pyfits.setval( out_name, 'USEAFTER', value= ' ' ) ### FIX
    pyfits.setval( out_name, 'PEDIGREE', value='INFLIGHT' ) 
    pyfits.setval( out_name, 'DESCRIP', value= 'Superbias created by J. Ely') ### FIX
    pyfits.setval( out_name, 'NEXTEND', value=3 ) 
    pyfits.setval( out_name, 'COMMENT', value='Replace this sometime' ) 
    
#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_refbias( glob.glob(sys.argv[1]), sys.argv[2] )
