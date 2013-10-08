"""

Script used to create monthly reference biases for the STIS Darks and Biases 
reference file pipeline

"""

import REFSTIS_functions
import shutil

#-------------------------------------------------------------------------------

def make_refbias( input_list, refbias_name='basebias.fits' ):
    """ Basic premise for making a refbias:
    1- join imsets from each datset together into one large file
    2- combine and cosmic ray screen joined imset
    3- divide by number of imsets
    4- set header keywords

    """

    print '#-------------------------------#'
    print '#        Running refbias        #'
    print '#-------------------------------#'
    print 'Making refbias %s' % (refbias_name)
    joined_out = refbias_name.replace('.fits', '_joined.fits' )

    print 'Joining images to %s' % joined_out
    REFSTIS_functions.msjoin( input_list, joined_out)

    print 'Checking for cosmic ray rejection'
    crj_filename = REFSTIS_functions.crreject( joined_out )

    shutil.copy( crj_filename, refbias_name)
 
    REFSTIS_functions.update_header_from_input( refbias_name, input_list )

    print 'Cleaning up...'
    REFSTIS_functions.RemoveIfThere( crj_filename )
    REFSTIS_functions.RemoveIfThere( joined_out )
    
    print 'refbias done'

#-------------------------------------------------------------------------------
