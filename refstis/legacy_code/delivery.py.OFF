#!/usr/bin/env python

'''
Functions necessary to check the reference files before delivery
to CDBS

'''

from astropy.io import fits as pyfits
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import shutil
import sys
from datetime import date

from crds import certify
import stistools
from stistools.calstis import calstis

from .functions import send_email

#----------------------------------------------------------------

def plot_obset(folder):
    '''
    Check collapsed columns and rows against last month for irregularities
    '''

    plt.ioff()

    print('#----------#')
    print('Making Plots')
    print('#----------#')
    bias = []
    biwk = []
    dark = []
    for filename in glob.glob(os.path.join(folder, '*.fits')):
        name = os.path.split(filename)[-1]
        if 'bias' in name and '_wk' in name:
            bias.append(filename)
        if 'bias' in name and '_biwk' in name:
            biwk.append(filename)
        if 'dark' in name and '_wk' in name:
            dark.append(filename)

    plt.rcParams['figure.subplot.hspace'] = .35
    plt.figure(figsize=(14, 20))
    plt.suptitle('Bias: collapsed rows, colums, and means')
    for i, ifile in enumerate(bias):
        print(ifile)
        data = pyfits.getdata(ifile, 1)
        plt.subplot(4, 1, 1)
        plt.plot( np.sum(data, axis=0), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(500, 3000)
        plt.xlabel('X pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 2)
        plt.plot( np.sum(data, axis=1), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(500, 3000)
        plt.xlabel('Y pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 3)
        plt.plot(i + 1, data.mean(), markersize=10, marker='d', label=ifile)
        plt.xlim(0, 9)
        plt.xlabel('Dataset')
        plt.ylabel('Mean')
        plt.subplot(4, 1, 4)
        plt.plot(i + 1, data.std(), markersize=10, marker='d', label=ifile)
        plt.xlim(0, 9)
        plt.xlabel('Dataset')
        plt.ylabel('Std')
    plt.savefig(os.path.join(folder, 'biases.pdf'))
    plt.close()

    plt.figure(figsize=(14, 20))
    plt.suptitle('Dark: collapsed rows, colums, and means')
    for i, ifile in enumerate(dark):
        print(ifile)
        data = pyfits.getdata(ifile, 1)
        plt.subplot(4, 1, 1)
        plt.plot( np.sum(data, axis=0), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(0, 200)
        plt.xlabel('X pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 2)
        plt.plot( np.sum(data, axis=1), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(0, 200)
        plt.xlabel('Y pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 3)
        plt.plot(i + 1, data.mean(), markersize=10, marker='d', label=ifile)
        plt.xlim(0, 9)
        plt.xlabel('Dataset')
        plt.ylabel('Mean')
        plt.subplot(4, 1, 4)
        plt.plot(i + 1, data.std(), markersize=10, marker='d', label=ifile)
        plt.xlim(0, 9)
        plt.xlabel('Dataset')
        plt.ylabel('Std')
    plt.savefig(os.path.join(folder, 'darks.pdf'))
    plt.close()

    plt.figure(figsize=(14, 20))
    plt.suptitle('BiWeek bias: collapsed rows, colums, and means')
    for i, ifile in enumerate(biwk):
        print(ifile)
        data = pyfits.getdata(ifile, 1)
        plt.subplot(4, 1, 1)
        plt.plot( np.sum(data, axis=0), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(2000, 6000)
        plt.xlabel('X pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 2)
        plt.plot( np.sum(data, axis=1), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(1000, 5000)
        plt.xlabel('Y pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 3)
        plt.plot(i + 1, data.mean(), markersize=10, marker='d', label=ifile)
        plt.xlim(0, 5)
        plt.xlabel('Dataset')
        plt.ylabel('Mean')
        plt.subplot(4, 1, 4)
        plt.plot(i + 1, data.std(), markersize=10, marker='d', label=ifile)
        plt.xlim(0, 5)
        plt.xlabel('Dataset')
        plt.ylabel('Std')
    plt.savefig(os.path.join(folder, 'biweeks.pdf'))
    plt.close()

#----------------------------------------------------------------

def set_descrip(folder):
    """ Make sure the descriptions are useful.  Should be removed when using
        the new version of the pipeline.
    """
    print('Making headers pretty')
    print('WARNING: Make sure to take this out when using the new refstis')

    for item in glob.glob(os.path.join(folder, '*.fits')):
        hdu = pyfits.open(item)

        gain = hdu[0].header['CCDGAIN']
        useafter = hdu[0].header['USEAFTER'][:9]

        if '_drk' in item:
            descrip = 'Weekly Dark for STIS CCD data taken after %s'% useafter
        elif '_bia' in item and hdu[0].header['CCDGAIN'] == 1:
            descrip = 'Weekly Gain=%d Bias for STIS CCD data taken after %s' % (gain, useafter)
        elif '_bia' in item and hdu[0].header['CCDGAIN'] == 4:
                descrip = 'Bi-Weekly Gain=%d Bias for STIS CCD data taken after %s' % (gain, useafter)

        while len(descrip) < 67:
            descrip += '-'

        hdu[0].header['DESCRIP'] = descrip
        hdu.writeto(item, overwrite=True)

#----------------------------------------------------------------

def regress(folder):
    """ Run `*drk` and `*bia` files in folder through CalSTIS to check
    for errors in processing

    """

    start_dir = os.getcwd()

    print('#------------------#')
    print('Running regression for')
    print(folder)
    print('#------------------#')

    monitor_dir = '/grp/hst/stis/darks_biases'
    test_suite = os.path.join(monitor_dir, 'test_suite')
    test_dark = os.path.join(monitor_dir, 'test_dark')

    print((glob.glob(os.path.join(folder, '*bia.fits'))))
    reference_files = glob.glob(os.path.join(folder, '*bia.fits')) + \
                    glob.glob(os.path.join(folder, '*drk.fits'))
    print(reference_files)
    assert len(reference_files) >= 1, 'No reference files in folder'

    print('Copying files and removing old files')
    for testing_dir in [test_suite, test_dark]:
        for oldfile in glob.glob(os.path.join(testing_dir, '*_drk.fits')) + \
                glob.glob(os.path.join(testing_dir, '*_bia.fits')):
            print(('removing', oldfile))
            os.remove(os.path.join(testing_dir, oldfile))

        for newfile in reference_files:
            print(('moving', newfile, '-->', testing_dir))
            shutil.copy(newfile, testing_dir)


    #######################################
    # Run checks in the test_suite folder #
    #######################################


    os.chdir(test_suite)
    print((os.getcwd()))

    bias_biwk_refs = glob.glob('*bias*_bi*.fits')
    bias_biwk_refs.sort()
    biasrefs = glob.glob('*bias*_wk*.fits')
    biasrefs.sort()
    darkrefs = glob.glob('*dark*.fits')
    darkrefs.sort()

    raws = glob.glob('*raw.fits')

    print('Setting IMPHTTAB in datasets')

    for dark, bias in zip(darkrefs, biasrefs):
        remove_products()

        print('#-------------------------------------------#')
        print(('Running CalSTIS with %s %s ' % (dark, bias)))
        print('#-------------------------------------------#')

        for rawfile in raws:
            if os.path.exists( rawfile.replace('raw.fits', 'wav.fits') ):
                wavefile = rawfile.replace('raw.fits', 'wav.fits')
            else:
                wavefile = ''

            pyfits.setval(rawfile, 'DARKFILE', value=dark, ext=0)
            pyfits.setval(rawfile, 'BIASFILE', value=bias, ext=0)
            pyfits.setval(rawfile, 'IMPHTTAB', value='oref$x9r1607mo_imp.fits', ext=0)
            if wavefile:
                pyfits.setval(wavefile, 'DARKFILE', value=dark, ext=0)
                pyfits.setval(wavefile, 'BIASFILE', value=bias, ext=0)
                pyfits.setval(wavefile, 'IMPHTTAB', value='oref$x9r1607mo_imp.fits', ext=0)

            status = calstis(rawfile, wavecal=wavefile)

            if status: sys.exit('Calstis Error detected for %s' % (dark[5:9]))

    ######################################
    # Run checks in the test_dark folder #
    ######################################

    os.chdir(test_dark)
    print((os.getcwd()))

    raws = glob.glob('*raw.fits')

    print('Setting IMPHTTAB in datasets')
    for item in raws:
        pyfits.setval(item, 'IMPHTTAB', value='oref$x9r1607mo_imp.fits', ext=0)

    for bias in bias_biwk_refs:

        remove_products()

        print('#------------------------------------------#')
        print(('Running CalSTIS with %s' % (bias)))
        print('#------------------------------------------#')

        for rawfile in raws:

            pyfits.setval(rawfile, 'BIASFILE', value=bias, ext=0)
            pyfits.setval(rawfile, 'DARKFILE', value=darkrefs[0], ext=0)

            status = calstis(input=rawfile)

            if status: sys.exit('Calstis Error detected for %s' % (dark[5:9]))


    os.chdir(start_dir)

#----------------------------------------------------------------


def send_forms(folder):

    start_dir = os.getcwd()
    os.chdir(folder)

    today_obj = date.today()
    today = str(today_obj.month) + '/' + \
        str(today_obj.day) + '/' + str(today_obj.year)
    message = '1-Name of deliverer: Allyssa Riley\n'
    message += ' (other e-mail addresses) debes@stsci.edu,\n'
    message += '\n'
    message += ' 2-Date of delivery: ' + today + '\n'
    message += '\n'
    message += ' 3-Instrument: STIS\n'
    message += '\n'
    message += ' 4-Type of file (bias,pht,etc.): bia, drk\n'
    message += '\n'
    message += ' 5-Has HISTORY section in header [0] been updated to describe in detail\n'
    message += '   why it is being delivered and how the file was created? (yes/no): yes\n'
    message += '\n'
    message += ' 6-USEAFTER, PEDIGREE, DESCRIP, and COMMENT have been checked? yes\n'
    message += '\n'
    message += ' 6a-Was the DESCRIP keyword updated with a summary of why the file was updated or created? \n'
    message += '   (yes/no) yes\n'
    message += '\n'
    message += ' 6b-If the reference files are replacing previous versions, do the new USEAFTER dates \n'
    message += '    exactly match the old ones? N/A\n'
    message += '\n'
    message += ' 7-CDBS Verification complete? (fitsverify,certify,etc.): yes\n'
    message += '\n'
    message += ' 8-Should these files be ingested in the OPUS, DADS and CDBS databases? yes\n'
    message += '   (if not indicate it clearly which ones):\n'
    message += '\n'
    message += ' 8a-If files are synphot files, should they be delivered to ETC? N/A\n'
    message += '\n'
    message += ' 9-Files run through CALXXX or SYNPHOT in the IRAF version of STSDAS and the IRAF* \n'
    message += '   version used by the Archive pipeline? (yes/no): yes \n'
    message += '   List the versions used:  calstis v {} \n'.format(stistools.calstis.__version__ )
    message += '\n'
    message += ' 10-Does it replace an old reference file? (yes/no): no\n'
    message += '\n'
    message += ' 10a-If yes, which one? \n'
    message += '     (If the file being replaced is bad, and should not be used with any data, please\n'
    message += '      indicate this here.)\n'
    message += '\n'
    message += ' 11- What is the level of change of the file? (e.g. compared to old file it\n'
    message += '     could be: SEVERE, MODERATE, TRIVIAL, 1\%, 5\% etc.): SEVERE\n'
    message += '\n'
    message += ' 11a-If files are tables, please indicate exactly which rows have changed. Show output \n'
    message += '     of compare_table.pro. N/A\n'
    message += '\n'
    message += ' 12-Please indicate which modes (e.g. all the STIS, FUVMAMA, E140L modes) are affected by\n'
    message += '     the changes in the file.  All CCD modes are affected \n'
    message += '\n'
    message += ' 13-Description of how the files were "tested" for correctness: Used calstis v {} \n'.format(stistools.calstis.__version__ )
    message += ' to reduce a test suite of CCD data and reduced a test suite of dark images as if \n'
    message += ' they were science images. The reduced darks were significantly lower in median and \n'
    message += ' mean values. The CCD images were reduced to either flt, crj, x1d, x2d, sx1, and sx2 \n'
    message += ' as appropriate. \n'
    message += '\n'
    message += ' 14-Additional Considerations: Some of the useafter dates DO NOT match the first date \n'
    message += ' in the pedigree. This is fine as the pipeline that creates the superdarks and \n'
    message += ' superbiases pulls the dates and times from the anneal proposal.\n'
    message += '\n'
    message += ' 15-Reason for delivery: New weekly biases and darks were created for the new anneal \n'
    message += ' month and need to be delivered for GO observations. \n'
    message += '\n '

    for search_string in ('*drk.fits', 'bias_wk*.fits', 'bias_bi*.fits'):
        file_list = glob.glob(search_string)
        file_list.sort()
        USEAFTER = []
        for item in file_list:
            USEAFTER.append(pyfits.getval(item, 'USEAFTER')[:11])
        for i, name in enumerate(file_list):
            if name != file_list[-1]:
                add_str = name + ' is for dates ' + \
                    USEAFTER[i] + ' to ' + USEAFTER[i + 1] + ' \n '
            else:
                add_str = name + ' is for dates after ' + USEAFTER[i] + '\n '
            message += add_str

    message += '\n\n 16-Disk location and name of files:\n'
    message += os.getcwd() + '\n'
    os.system('ls -la *.fits > tmp.txt')
    tmp = open('tmp.txt', 'r')
    lines = tmp.readlines()
    for line in lines:
        message += line
    os.remove('tmp.txt')

    delivery_form = open('deliveryform.txt', 'w')
    delivery_form.write(message)

    send_email(subject='STIS Darks and Bias Delivery Form', message=message)

    os.chdir(start_dir)

#----------------------------------------------------------------

def remove_products():
    ext_list = ['*_crj*', '*_flt*', '*_sx1*', '*_sx2*', '*_x1d*', '*_x2d*', '*_tmp*']
    for ext in ext_list:
        file_list = glob.glob(ext)
        if file_list != []:
            print(('removing {} files'.format( ext )))
            for file in file_list:
                os.remove(file)

#----------------------------------------------------------------

def move(source, destination):
    if not os.path.exists(destination):
        os.mkdir(destination)

    for root, dirs, files in os.walk(source):
        for filename in files:
            if not filename.endswith('.fits'):
                continue

            full_path = os.path.join(root, filename)

            if 'biases/4-1x1/' in full_path and 'weekbias_' in filename:
                full_destination = os.path.join(destination, filename)
            elif 'biases/1-1x1/' in full_path and 'weekbias_' in filename and not 'grp' in filename:
                full_destination = os.path.join(destination, filename)
            elif 'darks' in full_path and 'weekdark_' in filename:
                full_destination = os.path.join(destination, filename)
            else:
                continue

            print(full_path, '-->', full_destination)
            if os.path.exists(full_destination):
                os.remove(full_destination)
            shutil.copy(full_path, full_destination)


#----------------------------------------------------------------

def run_crds_checks(folder):
    datasets = ' '.join(glob.glob(os.path.join(folder, '*.fits')))
    errors = certify.CertifyScript("crds.certify {}".format(datasets))()

    assert not errors, "{} error found running crds checks.".format(errors)

#----------------------------------------------------------------

def check_all(folder, delivery_dir):
    move(folder, delivery_dir)
    regress(delivery_dir)
    plot_obset(delivery_dir)
    run_crds_checks(delivery_dir)
    send_forms(delivery_dir) #AER 25 Aug 2016: Moved this here from before regress(delivery_dir)

    print('#-------------------------------------------#')
    print('Darks and Bias Monitor complete.  ')
    print('Please run certify and fitsverify')
    print('Please send delivery form to redcat@stsci.edu.')
    print('#-------------------------------------------#')

#----------------------------------------------------------------
