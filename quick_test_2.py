import REFSTIS_basejoint
import glob

REFSTIS_basejoint.make_basebias( glob.glob('test/*raw.fits') ) 
