from refstis import functions
import numpy as np


def test_residuals_image():
    """Run a couple tests to check that residual columms images are
    being made correctly
    """

    res_image = functions.make_resicols_image(np.ones((10, 10)))
    assert (res_image == np.ones((10, 10))).all(), 'Error in simple array'

    sample_image = np.arange(20).repeat(20).reshape((20, 20))
    res_image = functions.make_resicols_image(sample_image)
    assert (res_image == np.ones((20, 20)) * 9.5).all(), 'Error in interesting array'
