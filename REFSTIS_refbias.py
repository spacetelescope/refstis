from pyraf import iraf
from pyraf.irafglobals import *
import REFSTIS_functions
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

    print '#-------------------------------#'
    print '#        Running refbias        #'
    print '#-------------------------------#'


    print 'Making refbias %s'%(refbias_name)

    refbias_name.replace('.fits','') # ensure its just a rootname
    refbias_path = os.path.split( refbias_name )[0]

    gain = REFSTIS_functions.get_keyword( bias_list, 'CCDGAIN', 0)
    xbin = REFSTIS_functions.get_keyword( bias_list, 'BINAXIS1', 0)
    ybin = REFSTIS_functions.get_keyword( bias_list, 'BINAXIS2', 0)

    REFSTIS_functions.split_images( bias_list )
    
    print 'Joining images'
    msjoin_list = ','.join( [ item for item in 
                              glob.glob( os.path.join(refbias_path,'*raw??.fits') ) ] )# if item[:9] in bias_list] )
    n_imsets = len(msjoin_list)
    joined_out = refbias_name+ '_joined.fits' 
    print joined_out
    
    msjoin_list_name = os.path.join( refbias_path,'msjoin.txt')
    msjoin_file = open( msjoin_list_name,'w')
    msjoin_file.write( '\n'.join(msjoin_list.split(',')) )
    msjoin_file.close()

    iraf.chdir( refbias_path )
    iraf.msjoin( inimg='@%s'%(msjoin_list_name), outimg=joined_out)#, Stderr='dev$null')

    #pdb.set_trace()
    crj_filename = REFSTIS_functions.crreject( joined_out )
    out_name = refbias_name + '.fits' 
    REFSTIS_functions.RemoveIfThere( out_name )
    shutil.copy( crj_filename, out_name)

    pyfits.setval( out_name, 'FILENAME', value=out_name )
    pyfits.setval( out_name, 'FILETYPE', value='CCD BIAS IMAGE')
    pyfits.setval( out_name, 'CCDGAIN', value= gain )
    pyfits.setval( out_name, 'BINAXIS1', value= xbin )
    pyfits.setval( out_name, 'BINAXIS2', value= ybin ) 
    pyfits.setval( out_name, 'USEAFTER', value= ' ' ) ### FIX
    pyfits.setval( out_name, 'PEDIGREE', value='INFLIGHT' ) 
    pyfits.setval( out_name, 'DESCRIP', value= 'refbias created by J. Ely') ### FIX
    pyfits.setval( out_name, 'NEXTEND', value= 3 ) 
    pyfits.setval( out_name, 'COMMENT', value='Replace this sometime' ) 

    for item in msjoin_list.split(','):
        os.remove(item)
        
    REFSTIS_functions.RemoveIfThere( msjoin_list_name )
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )
    for item in msjoin_list.split(','):
        REFSTIS_functions.RemoveIfThere( item ) 
    
#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_refbias( glob.glob(sys.argv[1]), sys.argv[2] )
