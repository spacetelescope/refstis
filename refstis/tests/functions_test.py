import REFSTIS_functions
import numpy as np

#-------------------------------------------------------------------------------

def test_residuals_image():
    """ Run a couple tests to check that residual columms images are
    being made correctly
    
    """

    res_image = REFSTIS_functions.make_resicols_image( np.ones( (10,10) ) )
    assert (res_image == np.ones( (10,10) )).all() , 'Error in simple array'

    sample_image = np.arange( 20 ).repeat( 20 ).reshape( (20,20) )
    res_image = REFSTIS_functions.make_resicols_image( sample_image )
    assert (res_image == np.ones( (20,20) )*9.5).all(), 'Error in interesting array'

#-------------------------------------------------------------------------------

def test_figure_number_of_periods():
    pass

#-------------------------------------------------------------------------------

def test_figure_days_in_period():
    pass

#-------------------------------------------------------------------------------

def test_translate_date_string():
    pass

#-------------------------------------------------------------------------------
