import os
import shutil
import numpy as np
from astropy.io import fits
from refstis import functions


def test_residuals_image():
    """Run a couple tests to check that residual columms images are
    being made correctly
    """

    res_image = functions.make_resicols_image(np.ones((10, 10)))
    assert (res_image == np.ones((10, 10))).all(), 'Error in simple array'

    sample_image = np.arange(20).repeat(20).reshape((20, 20))
    res_image = functions.make_resicols_image(sample_image)
    assert (res_image == np.ones((20, 20)) * 9.5).all(), 'Error in interesting array'


def test_refaver(tmp_path, basedark):
    """Tests refaver and call to msarith.
    """
    outfile = str(tmp_path / 'avg_dark.fits')

    input_files = []
    for i in range(1, 2 + 1):
        input_file = os.path.join(tmp_path, f'basedark_{i:.0f}.fits')
        input_files.append(input_file)
        shutil.copy(str(basedark), input_file)

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        functions.refaver(input_files, combined_name=outfile)
    finally:
        os.chdir(cwd)

    fits.info(outfile)

    with fits.open(outfile) as f:
        assert f[1].data[500, 500] == 0.
        assert f[1].data[592, 593] == 10.
