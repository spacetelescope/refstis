import pytest
import os
import tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits


@pytest.fixture
def basebias():
    return str(Path(__file__).parent / 'data' / 'basebias.fits')


@pytest.fixture
def basedark():
    return str(Path(__file__).parent / 'data' / 'basedark.fits')


@pytest.fixture
def rawbias_list(tmp_path):
    filename1 = tmp_path / "test1_raw.fits"
    filename2 = tmp_path / "test2_raw.fits"
    filename3 = tmp_path / "test3_raw.fits"

    header = {
        'TARGNAME': 'BIAS',
        'PROPOSID': 12345,
        'PROPTTL1': 'title',
        'TELESCOP': 'HST',
        'INSTRUME': 'STIS',
        'ROOTNAME': 'test.fits',
        'DETECTOR': 'CCD',
        'OPT_ELEM': 'MIRVIS',
        'APERTURE': 'F28X50LP',
        'PROPAPER': 'F28X50LP',
        'CRSPLIT': 1,
        'RA_TARG': 0.,
        'DEC_TARG': 0.,
        'EXTEND': True,
        'NEXTEND': 9,
        'GROUPS': False,
        'NRPTEXP': 3,
        'DQICORR': 'PERFORM',
        'CRCORR': 'PERFORM',
        'BLEVCORR': 'PERFORM',
        'BIASCORR': 'OMIT',
        'DARKCORR': 'OMIT',
        'FLATCORR': 'OMIT',
        'HELLCORR': 'OMIT',
        'ATODGAIN': 1.,
        'CCDGAIN': 1,
        'CCDOFFST': 3,
        'CCDAMP': 'D',
        'OBSTYPE': 'SPECTROSCOPIC',
        'OBSMODE': 'ACCUM',
        'CENTERA1': 532,
        'CENTERA2': 523,
        'SIZAXIS1': 1062,
        'SIZAXIS2': 1044,
        'BINAXIS1': 1,
        'BINAXIS2': 1,
        'NCOMBINE': 1,
        'BPIXTAB':  'oref$h1v11475o_bpx.fits',
        'DARKFILE': 'oref$a5i1349do_drk.fits',
        'PFLTFILE': 'oref$h4s1351lo_pfl.fits',
        'LFLTFILE': 'oref$jaj1058ho_lfl.fits',
        'PHOTTAB':  'oref$l7a15023o_pht.fits',
        'IMPHTTAB': 'oref$97a1641fo_imp.fits',
        'APERTAB':  'oref$y2r1559to_apt.fits',
        'CCDTAB':   'oref$16j1600do_ccd.fits',
        'BIASFILE': 'N/A',
        'CRREJTAB': 'oref$j3m1403io_crr.fits',
        'IDCTAB':   'oref$o8g1508do_idc.fits',
        'TDSTAB':   'oref$8712049eo_tds.fits',
        'TEXPTIME': 1.,
        'TEXPSTRT': 61158.25250479,
        'TEXPEND': 61158.252516364075,
        }

    sci_header = {
        'EXPTIME':    1.,
        'EXPSTART':   61158.25250479,
        'EXPEND':     61158.252516364075,
        'V_HELIO': 0.,
        }

    # ERR, DQ:
    empty_header = {
        'BITPIX':     16,
        'NPIX1':    1062,
        'NPIX2':    1044,
        'PIXVALUE':  0.0,
        }

    data1 = np.zeros((1044, 1062), dtype=np.uint16) + 1500
    data2 = np.zeros((1044, 1062), dtype=np.uint16) + 1501
    data3 = np.zeros((1044, 1062), dtype=np.uint16) + 1502

    data1[512:522, 512:522] += 10
    data2[512:522, 512:522] += 10
    data3[512:522, 512:522] += 10

    hdu = fits.HDUList([
        fits.PrimaryHDU(header=fits.Header(header)),
        fits.ImageHDU(data=data1,
            header=fits.Header(sci_header), name='SCI', ver=1),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='ERR', ver=1),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='DQ', ver=1),
        fits.ImageHDU(data=data2,
            header=fits.Header(sci_header), name='SCI', ver=2),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='ERR', ver=2),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='DQ', ver=2),
        fits.ImageHDU(data=data3,
            header=fits.Header(sci_header), name='SCI', ver=3),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='ERR', ver=3),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='DQ', ver=3),])

    hdu.writeto(filename1)
    hdu.writeto(filename2)
    hdu.writeto(filename3)

    return [str(filename1), str(filename2), str(filename3)]


@pytest.fixture
def rawdark_list(tmp_path):
    filename1 = tmp_path / "testdark1_raw.fits"
    filename2 = tmp_path / "testdark2_raw.fits"
    filename3 = tmp_path / "testdark3_raw.fits"

    header = {
        'TARGNAME': 'DARK',
        'PROPOSID': 12345,
        'PROPTTL1': 'title',
        'TELESCOP': 'HST',
        'INSTRUME': 'STIS',
        'ROOTNAME': 'test.fits',
        'DETECTOR': 'CCD',
        'OPT_ELEM': 'MIRVIS',
        'APERTURE': 'F28X50LP',
        'PROPAPER': 'F28X50LP',
        'CRSPLIT': 1,
        'RA_TARG': 0.,
        'DEC_TARG': 0.,
        'EXTEND': True,
        'NEXTEND': 3,
        'GROUPS': False,
        'NRPTEXP': 1,
        'DQICORR': 'PERFORM',
        'CRCORR': 'OMIT',
        'BLEVCORR': 'PERFORM',
        'BIASCORR': 'PERFORM',
        'DARKCORR': 'OMIT',
        'FLATCORR': 'OMIT',
        'HELLCORR': 'OMIT',
        'ATODGAIN': 1.,
        'CCDGAIN': 1,
        'CCDOFFST': 3,
        'CCDAMP': 'D',
        'OBSTYPE': 'SPECTROSCOPIC',
        'OBSMODE': 'ACCUM',
        'CENTERA1': 532,
        'CENTERA2': 523,
        'SIZAXIS1': 1062,
        'SIZAXIS2': 1044,
        'BINAXIS1': 1,
        'BINAXIS2': 1,
        'NCOMBINE': 1,
        'BPIXTAB':  'oref$h1v11475o_bpx.fits',
        'DARKFILE': 'ref$basedark.fits',
        'PFLTFILE': 'oref$h4s1351lo_pfl.fits',
        'LFLTFILE': 'oref$jaj1058ho_lfl.fits',
        'PHOTTAB':  'oref$l7a15023o_pht.fits',
        'IMPHTTAB': 'oref$97a1641fo_imp.fits',
        'APERTAB':  'oref$y2r1559to_apt.fits',
        'CCDTAB':   'oref$16j1600do_ccd.fits',
        'BIASFILE': 'ref$basebias.fits',
        'CRREJTAB': 'oref$j3m1403io_crr.fits',
        'IDCTAB':   'oref$o8g1508do_idc.fits',
        'TDSTAB':   'oref$8712049eo_tds.fits',
        'TEXPTIME': 1.,
        'TEXPSTRT': 61158.25250479,
        'TEXPEND': 61158.252516364075,
        }

    sci_header = {
        'EXPNAME': 'EXPNAME',
        'EXPTIME':    1.,
        'EXPSTART':   61158.25250479,
        'EXPEND':     61158.252516364075,
        'V_HELIO': 0.,
        'LTV1':     19.,
        'LTV2':     20.,
        'LTM1_1':    1.,
        'LTM2_2':    1.,
        'OCCDHTAV': 18.,
        }

    # ERR, DQ:
    empty_header = {
        'BITPIX':     16,
        'NPIX1':    1062,
        'NPIX2':    1044,
        'PIXVALUE':  0.0,
        'LTV1':     19.,
        'LTV2':     20.,
        'LTM1_1':    1.,
        'LTM2_2':    1.,
        }

    data1 = np.zeros((1044, 1062), dtype=np.uint16) + 1500

    data1[612:622, 612:622] += 10

    hdu = fits.HDUList([
        fits.PrimaryHDU(header=fits.Header(header)),
        fits.ImageHDU(data=data1,
            header=fits.Header(sci_header), name='SCI', ver=1),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='ERR', ver=1),
        fits.ImageHDU(data=None, header=fits.Header(empty_header), name='DQ', ver=1),])

    hdu.writeto(filename1)
    hdu.writeto(filename2)
    hdu.writeto(filename3)

    return [str(filename1), str(filename2), str(filename3)]
