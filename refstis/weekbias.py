"""Functions to create a weekly bias for the STIS instrument.

"""

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import shutil

from . import functions
from .basejoint import replace_hot_cols

#-------------------------------------------------------------------------------

def make_weekbias(input_list, refbias_name, basebias):
    """ Make 'weekly' bias from list of input bias files

    1. join imsets from each datset together into one large file
    2. combine and cosmic ray screen joined imset
    3. find hot colums
    4. add hot colums in to basebias sci data as output science data
    5. update error, dq, and headers

    .. note::

      For the SCI extensions, the baseline bias rate is taken from the aptly named
      basebias because anything that is "hot" is assumed to be potentially
      transient and is thus taken from the weekly biases.

      Update ERR extension of new superbias by assigning the ERR values of the
      baseline superbias except for the new hot pixels that are updated from
      the weekly superbias, for which the error extension of the weekly
      superbias is taken. Put the result in temporary ERR image.

    Parameters
    ----------
    input_list : list
        list of STIS bias files
    refbias_name : str
        filename of the output reference file
    basebias : str
        filename of the monthly basebias

    """

    print('#-------------------------------#')
    print('#        Running weekbias       #')
    print('#-------------------------------#')
    print('Output to %s' % (refbias_name))
    print('using {}'.format(basebias))

    joined_out = refbias_name.replace('.fits', '_joined.fits')
    functions.msjoin(input_list, joined_out)

    crj_filename = functions.crreject(joined_out)
    residual_image, median_image = functions.make_residual(crj_filename, (3, 15))

    resi_columns_2d = functions.make_resicols_image(residual_image, yfrac=.25)

    resi_mean, resi_median, resi_std = sigma_clipped_stats(resi_columns_2d[0],
                                                           sigma=3,
                                                           maxiters=20)
    replval = resi_mean + 5.0 * resi_std
    only_hotcols = np.where(resi_columns_2d >= replval, residual_image, 0)

    with fits.open(crj_filename, mode='update') as hdu:
        #-- update science extension
        baseline_sci = fits.getdata(basebias, ext=('sci', 1))
        hdu[('sci', 1)].data = baseline_sci + only_hotcols

        #-- update DQ extension
        hot_index = np.where(only_hotcols > 0)
        hdu[('dq', 1)].data[hot_index] = 16

        #- update ERR
        baseline_err = fits.getdata(basebias, ext=('err', 1))
        no_hot_index = np.where(only_hotcols == 0)
        hdu[('err', 1)].data[no_hot_index] = baseline_err[no_hot_index]

    shutil.copy(crj_filename, refbias_name)
    functions.update_header_from_input(refbias_name, input_list)
    fits.setval(refbias_name, 'TASKNAME', ext=0, value='WEEKBIAS')

    print('Cleaning up...')
    functions.RemoveIfThere(crj_filename)
    functions.RemoveIfThere(joined_out)

    print('weekbias done for {}'.format(refbias_name))

#-------------------------------------------------------------------------------
