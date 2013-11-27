#!/usr/bin/env python

'''
Functions necessary to check the reference files before delivery
to CDBS

'''

from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import shutil
import sys

from datetime import date, datetime

#----------------------------------------------------------------


def plot_files():
    '''
    Check collapsed columns and rows against last month for irregularities
    '''

    print '#----------#'
    print 'Making Plots'
    print '#----------#'
    bias = []
    biwk = []
    dark = []
    for top in [last, new]:
        for root, dirs, files in os.walk(top):
            if 'delivery' not in root:
                for ifile in files:
                    if 'bias_wk' in ifile:
                        bias.append(root + '/' + ifile)
                    if 'bias_biwk' in ifile:
                        biwk.append(root + '/' + ifile)
                    if 'dark_wk' in ifile:
                        dark.append(root + '/' + ifile)

    plt.rcParams['figure.subplot.hspace'] = .35
    plt.figure(figsize=(14, 20))
    plt.suptitle('Bias: collapsed rows, colums, and means')
    for i, ifile in enumerate(bias):
        print ifile
        data = pyfits.getdata(ifile, 1)
        plt.subplot(4, 1, 1)
        plt.plot(collapse_array(data, 0), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(500, 3000)
        plt.xlabel('X pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 2)
        plt.plot(collapse_array(data, 1), label=ifile)
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
    raw_input('Bias all ok?')
    plt.close()

    plt.figure(figsize=(14, 20))
    plt.suptitle('Dark: collapsed rows, colums, and means')
    for i, ifile in enumerate(dark):
        print ifile
        data = pyfits.getdata(ifile, 1)
        plt.subplot(4, 1, 1)
        plt.plot(collapse_array(data, 0), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(0, 200)
        plt.xlabel('X pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 2)
        plt.plot(collapse_array(data, 1), label=ifile)
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
    raw_input('Dark all ok?')
    plt.close()

    plt.figure(figsize=(14, 20))
    plt.suptitle('BiWeek bias: collapsed rows, colums, and means')
    for i, ifile in enumerate(biwk):
        print ifile
        data = pyfits.getdata(ifile, 1)
        plt.subplot(4, 1, 1)
        plt.plot(collapse_array(data, 0), label=ifile)
        plt.xlim(0, 1024)
        plt.ylim(2000, 6000)
        plt.xlabel('X pixels')
        plt.ylabel('Counts')
        plt.subplot(4, 1, 2)
        plt.plot(collapse_array(data, 1), label=ifile)
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
    raw_input('Biweek bias all ok?')
    plt.close()

#----------------------------------------------------------------

def calibrate( folder ):
    """ Run *drk and *bia files in folder through CalSTIS to check
    for errors in processing

    """
    monitor_dir = '/grp/hst/stis/darks_biases'
    test_suite = os.path.join( monitor_dir, 'test_suite' )
    test_dark = os.path.join( monitor_dir, 'test_dark' )

    
    reference_files = glob.glob(os.path.join( folder, '*bia.fits' ) ) + \
        glob.glob(os.path.join( folder, '*drk.fits' ) )

    for infile in reference_files:
        (path, filename) = os.path.split(infile)

        for location in [test_suite, test_dark]:
            testing_file = os.path.join( location, filename )
            if os.path.exists( testing_file ):
                os.remove( testing_file )

            shutil.copy( infile, testing_file )

    #######################################
    # Run checks in the test_suite folder #
    #######################################

    os.chdir( test_suite )
    bias_biwk_refs = glob.glob('*bias_bi*.fits')
    bias_biwk_refs.sort()
    biasrefs = glob.glob('bias*_wk*.fits')
    biasrefs.sort()
    darkrefs = glob.glob('dark*.fits')
    darkrefs.sort()
    raws = glob.glob('*raw.fits')
    wavs = glob.glob('*wav.fits')
    for dark, bias in zip(darkrefs, biasrefs):
        remove()
        for rawfile in raws:
            pyfits.setval(rawfile, 'DARKFILE', value=dark, ext=0)
            pyfits.setval(rawfile, 'BIASFILE', value=bias, ext=0)
        for wavefile in wavs:
            pyfits.setval(wavefile, 'DARKFILE', value=dark, ext=0)
            pyfits.setval(wavefile, 'BIASFILE', value=bias, ext=0)

        print '#-------------------------------------------#'
        print 'Running CalSTIS with %s %s ' % (dark, bias)
        print '#-------------------------------------------#'

        for txt_file in (dark[5:9] + '_err.txt', dark[5:9] + '_stdout.txt'):
            if os.path.exists(txt_file):
                os.remove(txt_file)
        iraf.calstis('*raw.fits', Stdout=dark[5:9] + '_stdout.txt')

        if not check_txt(dark[5:9] + '_stdout.txt'):
            sys.exit('Calstis Error detected for %s' % (dark[5:9]))

    ######################################
    # Run checks in the test_dark folder #
    ######################################

    os.chdir( test_dark )
    raws = glob.glob('*raw.fits')
    wavs = glob.glob('*wav.fits')
    print bias_biwk_refs
    for bias in bias_biwk_refs:
        remove()
        for rawfile in raws:
            pyfits.setval(rawfile, 'BIASFILE', value=bias, ext=0)
            pyfits.setval(rawfile, 'DARKFILE', value=darkrefs[0], ext=0)
        for wavefile in wavs:
            pyfits.setval(wavefile, 'BIASFILE', value=bias, ext=0)
            pyfits.setval(wavefile, 'DARKFILE', value=biasrefs[0], ext=0)

        print '#------------------------------------------#'
        print 'Running CalSTIS with %s' % (bias)
        print '#------------------------------------------#'

        name = str(bias)
        for txt_file in (name[5:13] + '_err.txt', name[5:13] + '_stdout.txt'):
            if os.path.exists(txt_file):
                os.remove(txt_file)
        iraf.calstis('*raw.fits,*wav.fits', Stdout=name[5:13] + '_stdout.txt')

        if not check_txt(name[5:13] + '_stdout.txt'):
            sys.exit('Calstis Error detected for %s' % (name[5:13]))

    print '#-------------------------------------------#'
    print 'Darks and Bias Monitor complete.  '
    print 'Please run certify and fitsverify'
    print 'Please send delievery form to cdbs@stsci.edu.'
    print '#-------------------------------------------#'

#----------------------------------------------------------------


def send_forms():
    today_obj = date.today()
    today = str(today_obj.month) + '/' + \
        str(today_obj.day) + '/' + str(today_obj.year)
    message = '1-Name of deliverer: Justin Ely\n'
    message += ' (other e-mail addresses) proffitt@stsci.edu,aloisi@stsci.edu,debes@stsci.edu,\n'
    message += 'osten@stsci.edu,bohlin@stsci.edu\n'
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
    message += '   List the versions used:  calstis v 2.36 \n'
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
    message += ' 12-Description of how the files were "tested" for correctness: Used calstis v 2.36 \n'
    message += ' to reduce a test suite of CCD data and reduced a test suite of dark images as if \n'
    message += ' they were science images. The reduced darks were significantly lower in median and \n'
    message += ' mean values. The CCD images were reduced to either flt, crj, x1d, x2d, sx1, and sx2 \n'
    message += ' as appropriate. \n'
    message += '\n'
    message += ' 13-Additional Considerations: Some of the useafter dates DO NOT match the first date \n'
    message += ' in the pedigree. This is fine as the pipeline that creates the superdarks and \n'
    message += ' superbiases pulls the dates and times from the anneal proposal.\n'
    message += '\n'
    message += ' 14-Reason for delivery: New weekly biases and darks were created for the new anneal \n'
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

    message += '\n\n 15-Disk location and name of files:\n'
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

#----------------------------------------------------------------


def check_txt(ifile):
    '''
    Check text file for any errors or failures.
    '''
    lines = open(ifile)
    'calstis0  failed'
    for line in lines.readlines():
        if (('calstis' in line) & ('failed' in line)):
            return False
    return True

#----------------------------------------------------------------

def remove():
    ext_list = ['*_crj*', '*_flt*', '*_sx1*', '*_sx2*', '*_x1d*', '*_x2d*', '*_tmp*']
    for ext in ext_list:
        file_list = glob.glob(ext)
        if file_list != []:
            print 'removing %s files'%(ext)
            for file in file_list:
                os.remove(file)

#----------------------------------------------------------------

def run_cdbs_checks():
    print 'Round one'
    os.system('python %slength_descrip_CDBS.py' % (monitor_dir))
    print '\n\nRound two'
    os.system('python %slength_descrip_CDBS.py' % (monitor_dir))

    for command in ['certify *.fits', 'fitsverify_delivery *.fits']:
        os.system(command)

#----------------------------------------------------------------


def check_all():
    pass

#----------------------------------------------------------------

if __name__ == "__main__":
    calibrate( sys.argv[1] )
