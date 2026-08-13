import os
from astropy.io import fits
from refstis.basejoint import make_basebias
from refstis.refbias import make_refbias
from refstis.weekbias import make_weekbias
from refstis.basedark import make_basedark
from refstis.weekdark import make_weekdark


def test_basebias(rawbias_list, tmp_path):
    '''Test refstis.basejoint.make_basebias().
    '''
    outfile = str(tmp_path / 'basebias.fits')
    make_basebias(input_list=rawbias_list, refbias_name=str(outfile))

    fits.info(outfile)

    with fits.open(outfile) as f:
        assert f[1].data[400, 400] == 0.
        assert f[1].data[493, 493] == 10.


def test_weekbias(rawbias_list, tmp_path, basebias):
    '''Test refstis.weekbias.make_weekbias() with stored output of test_basebias().
    '''
    outfile = str(tmp_path / 'weekbias.fits')
    make_weekbias(input_list=rawbias_list, refbias_name=str(outfile), basebias=basebias)

    fits.info(outfile)

    with fits.open(outfile) as f:
        assert f[1].data[400, 400] == 0.
        assert f[1].data[492, 493] == 90.


def test_make_refbias(rawbias_list, tmp_path):
    '''Test refstis.refbias.make_refbias().
    '''
    outfile = str(tmp_path / 'refbias.fits')
    make_refbias(input_list=rawbias_list, refbias_name=outfile)

    fits.info(outfile)

    with fits.open(outfile) as f:
        assert f[1].data[400, 400] == 0.
        assert f[1].data[492, 493] == 90.


def test_basedark(rawdark_list, tmp_path, basebias):
    '''Test refstis.basedark.make_basedark() with stored output of test_basebias().
    '''
    outfile = str(tmp_path / 'basedark.fits')
    make_basedark(input_list=rawdark_list, refdark_name=outfile, bias_file=basebias)

    fits.info(outfile)

    with fits.open(outfile) as f:
        assert f[1].data[500, 500] == 0.
        assert f[1].data[592, 593] == 10.


def test_weekdark(rawdark_list, tmp_path, basebias, basedark):
    '''Test refstis.weekdark.make_weekdark() with stored outputs of test_basebias()
    and test_basedark().
    '''
    outfile = str(tmp_path / 'weekdark.fits')
    make_weekdark(input_list=rawdark_list,
                  refdark_name=outfile,
                  thebasedark=basedark,
                  thebiasfile=basebias)

    fits.info(outfile)

    with fits.open(outfile) as f:
        assert f[1].data[500, 500] == 0.
        assert f[1].data[592, 593] == 10.
