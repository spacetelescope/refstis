from astropy.io import fits as pyfits
from astropy.time import Time
import numpy as np
import scipy.stats as stats
from matplotlib import pyplot
import subprocess
import os
import glob
import shutil
import pdb

def make_testing_directories(cur_dir, product_testing):
    '''
    This function creates a directory called testing and another directory under
    testing with the name of the item to be tested (either bias or dark)
    If a directory already exists, it is not overwritten
    '''
    if not os.path.exists(os.path.join(cur_dir, 'test')):
        os.mkdir(os.path.join(cur_dir, 'test'))
    testing_dir = os.path.join(cur_dir, 'test', product_testing)
    if not os.path.exists(testing_dir):
        os.mkdir(os.path.join(cur_dir, 'test', product_testing))

    return testing_dir

def copy_file_to_testing_dir(data_dir, testing_dir, filename = '*raw.fits'):
    '''
    Copy all raw files from data_dir to testing_dir
    Consider in the future copying wav files
    '''
    flist = glob.glob(os.path.join(data_dir, filename))
    filenames = []
    for ifile in flist:
        shutil.copy(ifile, os.path.join(testing_dir, ifile.split('/')[-1]))
        filenames.append(ifile.split('/')[-1])
    return filenames

def calibrate(calstis_func, filename, reffile_keyword, reffile_name, commandline_options = ''):
    '''
    This code calls one of the Calstis modules with user defined commandline_options
    to calibrate a file. This allows you to avoid messing with the calibration switches
    in the header
    '''
    pyfits.setval(filename, reffile_keyword, ext = 0, value = reffile_name)
    subprocess_str = '%s %s %s' %(calstis_func, commandline_options, filename)
    subprocess.call(subprocess_str, shell = True)

def make_fractional_year(date):
    date = np.array(date)
    date = Time(date, format = 'mjd', scale = 'utc').yday
    date = [float(i.split(':')[0]) + float(i.split(':')[1])/ 365.25 for i in date]
    return date
    

def calc_mean_level(test_dir, filetype = 'flt'):
    '''
    This code plots the STD of each calibrated image in test_dir 
    input:
    test_dir: name of directory where data are located
    ax: axis object to plot onto
    '''
    flist = glob.glob(os.path.join(test_dir, '*%s.fits' %(filetype)))
    date = []
    stdev = []
    median_val = []
    for ifile in flist:
        date.append(pyfits.getval(ifile, 'texpstrt', 0))
        tbdata = pyfits.getdata(ifile, 1)
        stdev.append(np.std(tbdata))
        median_val.append(np.median(tbdata))
    date = np.array(make_fractional_year(date))
    median_val = np.array(median_val)
    stdev = np.array(stdev)

    return date, median_val, stdev

def examine_stdev(date_flt, stdev_flt, date_raw, stdev_raw, savefig = False, fig_filename = None):
    #Test that a filename has been specified if savefig is turned on
    assert (savefig and fig_filename) or ((not savefig) and (not fig_filename)), 'You must specify both savefig and fig_filename to save a figure'

    if savefig:
        fig2 = pyplot.figure(2)
        ax2 = fig2.add_subplot(1, 1, 1)
        ax2.plot(date_flt, stdev_flt, 'bo')
        ax2.plot(date_raw, stdev_raw, 'rx')
        ax2.legend(['Standard Deviation of FLT image', 'Standard Deviation of RAW image'], loc = 'best')
        ax2.set_xlabel('Date')
        ax2.set_ylabel('Standard Deviation')
        ax2.set_title('Standard Deviation of RAW and FLT Images')
        pyplot.savefig(fig_filename)
        pyplot.close()
    if  not np.all(date_raw - date_flt == 0):
        'raw and flt stdev correspond to different dates, STDEV not compared'
        return
    else:
        if not np.all(stdev_raw - stdev_flt >= 0):
            
            print '\nWARNING: Standard Deviation for at least one file was larger in the FLT than in the RAW\n'
    return

def delete_calibration_products(flist, testing_dir):
    '''
    This function deletes all flt files so CalSTIS will create new products.
    At some point we may want to consider
    extending this to other calibration products
    '''
    for ifile in flist:
        if os.path.exists(os.path.join(testing_dir, ifile.replace('raw', 'flt'))):
            os.remove(os.path.join(testing_dir, ifile.replace('raw', 'flt')))

def calibrate_bias_data(testing_dir, cur_dir, data_dir, reffile_dir, reffile_name, ax):

    flist = copy_file_to_testing_dir(data_dir, testing_dir)
    copy_file_to_testing_dir(reffile_dir, testing_dir, filename = reffile_name)
    os.chdir(testing_dir)
    delete_calibration_products(flist, testing_dir)
    for ifile in flist:
        calibrate('cs1.e', ifile, 'BIASFILE', reffile_name, '-blev -bias -dqi')


def test_bias(cur_dir, data_dir, reffile_dir, reffile_name, ax, calibrate = False):
    testing_dir = make_testing_directories(cur_dir, 'bias')
    if calibrate:
        calibrate_bias_data(testing_dir, cur_dir, data_dir, reffile_dir, reffile_name, ax)
    date, median_val, stdev = calc_mean_level(testing_dir)
    ax.errorbar(date, median_val, stdev, fmt = 'o')
    ax.set_xlabel('Date')
    ax.set_ylabel('Median Value of FLT image')
    pyplot.savefig('bias_mean.pdf')
    date_raw, median_val_raw, stdev_raw = calc_mean_level(testing_dir, filetype = 'raw')
    examine_stdev(date, stdev, date_raw, stdev_raw, savefig = True, fig_filename = 'bias_stdev.pdf')
    os.chdir(cur_dir)
    
if __name__ == "__main__":
    os.environ['oref'] = '/grp/hst/cdbs/oref/'
    fig = pyplot.figure(1)
    ax = fig.add_subplot(1,1,1)
    os.chdir('/user/bostroem/science/cte/2012_04')
    cur_dir= os.getcwd()
    test_bias(cur_dir, '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk01', '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk01', 'refbias__wk01.fits', ax, calibrate = True)
    test_bias(cur_dir, '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk02', '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk02', 'refbias__wk02.fits', ax, calibrate = True)
    test_bias(cur_dir, '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk03', '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk03', 'refbias__wk03.fits', ax, calibrate = True)
    test_bias(cur_dir, '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk04', '/user/bostroem/science/cte/2012_04/bias/biases/1-1x1/wk04', 'refbias__wk04.fits', ax, calibrate = True)

    raw_input('Pause before existing to view plot')
    
