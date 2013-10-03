"""

Script used to create monthly reference biases for the STIS Darks and Biases reference
file pipeline

"""

import REFSTIS_functions
import pyfits
import sys
import glob
import shutil
import os

#------------------------------------------------------------------------------------

def update_header( filename ):
    pass

#------------------------------------------------------------------------------------

def make_refbias( input_list, refbias_name='basebias.fits' ):
    """ Basic premise for making a refbias:
    1- split all raw images into their raw imsets
    2- join imsets together into one large file
    3- combine and cosmic ray screen joined imset
    4- divide by number of imsets
    5- set header keywords

    """

    import os
    import REFSTIS_functions

    print '#-------------------------------#'
    print '#        Running refbias        #'
    print '#-------------------------------#'
    print 'Making refbias %s' % (refbias_name)

    refbias_path = os.path.split( refbias_name )[0]

    rootname_set = set( [ os.path.split(item)[1][:9] for item in input_list] )

    gain = REFSTIS_functions.get_keyword( input_list, 'CCDGAIN', 0)
    xbin = REFSTIS_functions.get_keyword( input_list, 'BINAXIS1', 0)
    ybin = REFSTIS_functions.get_keyword( input_list, 'BINAXIS2', 0)

    n_imsets = REFSTIS_functions.split_images( input_list )
    
    joined_out = refbias_name.replace('.fits', '_joined.fits' )
    print 'Joining images to %s' % joined_out
    msjoin_list = [ item for item in 
                    glob.glob( os.path.join(refbias_path,'*raw??.fits') )  if os.path.split(item)[1][:9] in rootname_set]
    REFSTIS_functions.msjoin( msjoin_list, joined_out)


    print 'Checking for cosmic ray rejection'
    crj_filename = REFSTIS_functions.crreject( joined_out )
    print 'CR-Rejected output will be %s' % crj_filename 

    REFSTIS_functions.RemoveIfThere( refbias_name )
    shutil.copy( crj_filename, refbias_name)

    pyfits.setval( refbias_name, 'FILENAME', value=refbias_name )
    pyfits.setval( refbias_name, 'FILETYPE', value='CCD BIAS IMAGE')
    pyfits.setval( refbias_name, 'CCDGAIN', value= gain )
    pyfits.setval( refbias_name, 'BINAXIS1', value= xbin )
    pyfits.setval( refbias_name, 'BINAXIS2', value= ybin ) 
    pyfits.setval( refbias_name, 'USEAFTER', value= ' ' ) ### FIX
    pyfits.setval( refbias_name, 'PEDIGREE', value='INFLIGHT' ) 
    pyfits.setval( refbias_name, 'DESCRIP', value= 'refbias created by J. Ely') ### FIX
    pyfits.setval( refbias_name, 'NEXTEND', value= 3 ) 
    pyfits.setval( refbias_name, 'COMMENT', value='Replace this sometime' ) 
 
    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )
    for item in msjoin_list:
        REFSTIS_functions.RemoveIfThere( item ) 
    
    print '#-------------------------------#'
    print '#        Finished refbias       #'
    print '#-------------------------------#'


#------------------------------------------------------------------------------------
    
if __name__ == "__main__":
    make_refbias( glob.glob(sys.argv[1]), sys.argv[2] )
