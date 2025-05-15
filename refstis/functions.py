from astropy.io import fits as pyfits
from astropy.stats import sigma_clipped_stats
import numpy as np
import os
import stat
import shutil
from scipy.signal import medfilt
from scipy.ndimage import median_filter
from astropy.time import Time
import math
import sqlite3
import datetime

from stistools.calstis import calstis
from stistools.ocrreject import ocrreject
from stistools.basic2d import basic2d
#--------------------------------------------------------------------------------

def send_email(subject=None, message=None, from_addr=None, to_addr=None):
    '''
    Send am email via SMTP server.
    This will not prompt for login if you are already on the internal network.
    '''
    import os
    import getpass
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    users_email=getpass.getuser()+'@stsci.edu'

    if not subject:
        subject='Message from %s'%(__file__)
    if not message:
        message='You forgot to put a message into me'
    if not from_addr:
        from_addr=users_email
    if not to_addr:
        to_addr=users_email

    svr_addr='smtp.stsci.edu'
    msg = MIMEMultipart()
    msg['Subject']=subject
    msg['From']=from_addr
    msg['To']=to_addr
    msg.attach(MIMEText(message))
    s = smtplib.SMTP(svr_addr)
    s.sendmail(from_addr, to_addr, msg.as_string())
    s.quit()
    print('\nEmail sent to %s \n' %(from_addr))

#---------------------------------------------------------------------------

def mjd_to_greg(mjd):
   #This comes from http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/julian-date-form
   JD = mjd + 2400000.5
   JD = int(JD)
   L= JD+68569
   N= 4*L/146097
   L= L-(146097*N+3)/4
   I= 4000*(L+1)/1461001
   L= L-1461*I/4+31
   J= 80*L/2447
   K= L-2447*J/80
   L= J/11
   J= J+2-12*L
   I= 100*(N-49)+I+L
   Year = I
   Month = J
   Day = K
   month_to_day = {'0': 0,'1':31, '2':59, '3':90, '4':120, '5':151, '6':181, '7':212, '8':243, '9':273, '10':304, '11':334, '12':365}
   tot_day = (month_to_day[str(int(Month)-1)] + Day)
   day_in_year = 365.0
   if (Month >= 2) & (Year%4 == 0): #for leap year
      tot_day = tot_day + 1.0
      day_in_year = 365.0 + 1.0
   frac_day = tot_day / day_in_year
   fractional_year = Year + frac_day
   return (Year, Month, Day, fractional_year)

#------------------------------------------------------------------------

def update_header_from_input(filename, input_list):
    """ Updates header of output file using keywords from the input data

    If a header keyword is not consistent in this step, an error will be
    raised

    """
    targname = get_keyword(input_list, 'TARGNAME', 0)
    if targname == 'BIAS':
        filetype = 'CCD BIAS IMAGE'
    elif targname == 'DARK':
        filetype = 'DARK IMAGE'
    else:
        raise ValueError('targname %s not understood' % str(targname))

    gain = get_keyword(input_list, 'CCDGAIN', 0)
    if gain == 1:
        frequency = 'Weekly'
        N_period = 4
    elif gain == 4:
        frequency = 'Bi-Weekly'
        N_period = 2
    else:
        raise ValueError( 'Frequency %s not understood' % str(frequency))

    nimsets = count_imsets(input_list)

    proposals = list(set([pyfits.getval(item, 'PROPOSID') for item in input_list]))
    prop_titles = list(set([pyfits.getval(item, 'PROPTTL1') for item in input_list]))

    data_start_pedigree, data_end_pedigree, data_start_mjd, data_end_mjd = get_start_and_endtimes(input_list)
    #anneal_weeks = divide_anneal_month(data_start_mjd, data_end_mjd, '/grp/hst/stis/calibration/anneals/', N_period)
    #useafter = 'value'
    #for begin, end in anneal_weeks:
    #    if (begin < data_start_mjd) and (end > data_end_mjd):
    #        begin_time = Time(begin, format = 'mjd', scale = 'utc').iso  #anneal week start
    #        useafter = datetime.datetime.strptime(begin_time.split('.')[0], '%Y-%m-%d %H:%M:%S').strftime('%b %d %Y %X')

    month, day, time, year = Time(data_start_mjd, format='mjd').datetime.ctime().split()[1:]
    useafter = '{:s} {:02d} {:s} {:s}'.format(month, int(day), year, time)

    hdu_out = pyfits.HDUList(pyfits.PrimaryHDU())
    hdu_out[0].header['FILENAME'] = os.path.split(filename)[1]
    hdu_out[0].header['NEXTEND'] = 3
    hdu_out[0].header['TELESCOP'] = 'HST'
    hdu_out[0].header['INSTRUME'] = 'STIS'
    hdu_out[0].header['DETECTOR'] = ('CCD', 'detector in use: CCD')

    if targname == 'BIAS':
        hdu_out[0].header['CCDAMP'] = get_keyword(input_list, 'CCDAMP', 0)
        hdu_out[0].header['CCDGAIN'] = gain
        hdu_out[0].header['CCDOFFST'] = get_keyword(input_list, 'CCDOFFST', 0)
    elif targname == 'DARK':
        hdu_out[0].header['CCDAMP'] = 'ANY'
        hdu_out[0].header['CCDGAIN'] = -1
        hdu_out[0].header['CCDOFFST'] = -1
    else:
        raise ValueError('{} for targname not understood'.format(targname))

    hdu_out[0].header['OPT_ELEM'] = 'ANY'
    hdu_out[0].header['APERTURE'] = 'ANY'
    hdu_out[0].header['OBSTYPE'] = 'ANY'
    hdu_out[0].header['BINAXIS1'] = get_keyword(input_list, 'BINAXIS1', 0)
    hdu_out[0].header['BINAXIS2'] = get_keyword(input_list, 'BINAXIS2', 0)
    hdu_out[0].header['FILETYPE'] = filetype
    hdu_out[0].header['PEDIGREE'] = 'INFLIGHT %s %s' % (data_start_pedigree, data_end_pedigree)
    hdu_out[0].header['USEAFTER'] = useafter
    hdu_out[0].header['DESCRIP'] = "%s gain=%d %s for STIS CCD data taken after %s" % (frequency, gain, targname.lower(), useafter[:11])
    while len(hdu_out[0].header['DESCRIP']) < 67:
        hdu_out[0].header['DESCRIP'] = hdu_out[0].header['DESCRIP'] + '-'
    if len(hdu_out[0].header['DESCRIP']) > 67:
        raise ValueError('DESCRIP is too long! {}'.format(hdu_out[0].header['DESCRIP']))

    hdu_out[0].header.add_comment('Reference file created by %s' % __name__ )

    hdu_out[0].header.add_history('Super{} image, combination of {} input {} frames taken in'.format(targname.lower(),
                                                                                                     nimsets,
                                                                                                     targname.lower()))
    hdu_out[0].header.add_history('CCDGAIN={}, BINAXIS1={}, BINAXIS2={} mode.'.format(hdu_out[0].header['CCDGAIN'],
                                                                                      hdu_out[0].header['BINAXIS1'],
                                                                                      hdu_out[0].header['BINAXIS2']))

    hdu_out[0].header.add_history('All input frames were from Proposal(s) {}:'.format('/'.join(map(str, proposals))))
    for item in prop_titles:
        hdu_out[0].header.add_history(item)

    hdu_out[0].header.add_history('The following input files were used:')
    for item in input_list:
        hdu_out[0].header.add_history(os.path.split(item)[1])
    hdu_out[0].header.add_history('')

    if targname == 'BIAS':
        hdu_out[0].header.add_history('The data were split into sub-lists of less than 30 ')
        hdu_out[0].header.add_history('imsets each. The script "refbias" was run on the')
        hdu_out[0].header.add_history('individual sub-lists and then averaged using the script')
        hdu_out[0].header.add_history('"refaver". The "refbias" procedure works as follows.')
        hdu_out[0].header.add_history('After joining the files together into a multi-imset')
        hdu_out[0].header.add_history('file, overscan subtraction is performed for every')
        hdu_out[0].header.add_history('individual bias frame. The bias frames in the multi-')
        hdu_out[0].header.add_history('imset file are then combined using ocrreject, which')
        hdu_out[0].header.add_history('performs a three-step iterative cosmic-ray rejection')
        hdu_out[0].header.add_history('with rejection thresholds of 5, 4, and 3 sigma,')
        hdu_out[0].header.add_history('respectively.')
        hdu_out[0].header.add_history('')
        hdu_out[0].header.add_history('After cosmic ray rejection, the combined bias was')
        hdu_out[0].header.add_history('divided by the number of frames combined (ocrreject')
        hdu_out[0].header.add_history('adds images instead of averaging them), using the ')
        hdu_out[0].header.add_history('stsdas.toolbox.imgtools.mstools.msarith routine.')
        hdu_out[0].header.add_history('Subsequently, hot pixels (and columns) were identified ')
        hdu_out[0].header.add_history('as follows. A median-filtered version of the averaged ')
        hdu_out[0].header.add_history('bias is created (kernel = 15 x 3 pixels) and')
        hdu_out[0].header.add_history('subtracted from the averaged bias, leaving a')
        hdu_out[0].header.add_history('"residual" image containing hot pixels and columns.')
        hdu_out[0].header.add_history('The pixels hotter than 5 sigma of the effective RMS')
        hdu_out[0].header.add_history('noise in the residual bias image are then identified,')
        hdu_out[0].header.add_history('and flagged as such in the Data Quality (DQ)')
        hdu_out[0].header.add_history('extension of the output reference bias file.')

    elif targname == 'DARK':
        hdu_out[0].header['REF_TEMP'] = 18
        hdu_out[0].header['DRK_VS_T'] = 0.07

        hdu_out[0].header.add_history('This superdark image is a combination of 2 images.')
        hdu_out[0].header.add_history('The first is a "baseline dark" image which is a  ')
        hdu_out[0].header.add_history('(typically) monthly average of (cosmic-ray-rejected)')
        hdu_out[0].header.add_history('dark files. The second image is made locally within the')
        hdu_out[0].header.add_history('"weekdark" script.')
        hdu_out[0].header.add_history('These gain {} darks are first'.format(hdu_out[0].header['CCDGAIN']))
        hdu_out[0].header.add_history('corrected for temperature (C) using the factor')
        hdu_out[0].header.add_history('1.0 / (1.0 + DRK_V_TMP * (OCCDHTAV - S2_REF_TEMP)) and then ')
        hdu_out[0].header.add_history('combined together (cosmic rays are rejected in the')
        hdu_out[0].header.add_history('combination) using calstis,')
        hdu_out[0].header.add_history('and normalized to a dark time of 1 second. After that,')
        hdu_out[0].header.add_history('hot pixels in that normalized dark are updated into the ')
        hdu_out[0].header.add_history('baseline dark. These hot pixels have a value higher than ')
        hdu_out[0].header.add_history('(baseline dark current + 5 sigma of the new dark).')
        hdu_out[0].header.add_history('The pixels hotter than 0.1 electron/sec in this dark')
        hdu_out[0].header.add_history('are being assigned a DQ value of 16.')
        hdu_out[0].header.add_history('Previously hot pixels that have fully annealed out ')
        hdu_out[0].header.add_history('between the observing dates of the baseline and the')
        hdu_out[0].header.add_history('new dark are being assigned a value equal to that in')
        hdu_out[0].header.add_history('a median-filtered (kernel = 5x5 pixels) version of')
        hdu_out[0].header.add_history('the baseline dark.')

    with pyfits.open(filename) as ref:
        hdu_out.append(pyfits.ImageHDU(data=ref[1].data))
        hdu_out[1].header['EXTNAME'] = 'SCI'
        hdu_out[1].header['EXTVER'] = 1
        hdu_out[1].header['PCOUNT'] = 0
        hdu_out[1].header['GROUNT'] = 1

        hdu_out.append(pyfits.ImageHDU(data=ref[2].data))
        hdu_out[2].header['EXTNAME'] = 'ERR'
        hdu_out[2].header['EXTVER'] = 1
        hdu_out[2].header['PCOUNT'] = 0
        hdu_out[2].header['GROUNT'] = 1

        hdu_out.append(pyfits.ImageHDU(data=ref[3].data))
        hdu_out[3].header['EXTNAME'] = 'DQ'
        hdu_out[3].header['EXTVER'] = 1
        hdu_out[3].header['PCOUNT'] = 0
        hdu_out[3].header['GROUNT'] = 1

    hdu_out.writeto(filename, overwrite=True, output_verify='exception')

#------------------------------------------------------------------------

def get_start_and_endtimes(input_list):
    times = []
    for ifile in input_list:
        times.append(pyfits.getval(ifile, 'texpstrt', 0))
        times.append(pyfits.getval(ifile, 'texpend', 0))
    times.sort()
    start_mjd = times[0]
    end_mjd = times[-1]
    times = np.array([x.split()[0] for x in Time(times, format = 'mjd', scale = 'utc').iso])

    start_list = times[0].split('-')
    end_list = times[-1].split('-')

    #return strings in format mm/dd/yyyy
    start_str = '%s/%s/%s' %( start_list[2], start_list[1], start_list[0])
    end_str = '%s/%s/%s' %(end_list[2], end_list[1], end_list[0])
    return start_str, end_str, start_mjd, end_mjd

#------------------------------------------------------------------------

def make_resicols_image(residual_image, yfrac=1):
    print("Making residual column image")

    ystart = 0
    yend = min(1024, int(np.floor(yfrac * residual_image.shape[0] + .5)))

    residual_columns = np.mean(residual_image[ystart:yend], axis=0)
    print(ystart, '-->', yend)
    residual_columns_image = residual_columns * np.ones(residual_image.shape[1])[:, np.newaxis]

    return residual_columns_image

#------------------------------------------------------------------------

def make_residual(mean_bias, kern=(3, 15)):
    """Create residual image

    Median filter the median with a 15 x 3 box and subtract from the mean
    to produce the residual image.


    """
    mean_hdu = pyfits.open(mean_bias)
    mean_image = mean_hdu[('sci', 1)].data

    median_image = median_filter(mean_image, kern)

    medi_mean = sigma_clipped_stats(median_image, sigma=3, maxiters=40)[0]
    mean_mean = sigma_clipped_stats(mean_image, sigma=3, maxiters=40)[0]
    diffmean = mean_mean - medi_mean

    median_image += diffmean
    residual_image = mean_image - median_image

    return residual_image, median_image

#------------------------------------------------------------------------

def normalize_crj(filename):
    """ Normalize the input filename by exptim/gain and flush hdu

    """

    with pyfits.open(filename, mode='update') as hdu:
        exptime = hdu[0].header['TEXPTIME']
        gain = hdu[0].header['ATODGAIN']

        norm_factor = float(exptime)/gain
        print('Normalizing by ', norm_factor)
        hdu[('sci', 1)].data /= norm_factor
        hdu[('err', 1)].data /= abs(norm_factor)

        hdu[0].header['TEXPTIME'] = 1

#------------------------------------------------------------------------

def msjoin(imset_list, out_name='joined_out.fits'):
    """ Replicate msjoin functionality in pure python

    """

    hdu = pyfits.open( imset_list[0] )

    ext_count = 0
    n_offset = (len( hdu[1:] ) // 3) + 1
    for dataset in imset_list[1:]:
        add_hdu = pyfits.open( dataset )
        for extension in add_hdu[1:]:
            extension.header['EXTVER'] = (ext_count // 3) + n_offset
            hdu.append( extension )
            ext_count += 1

    hdu[0].header['NEXTEND'] = len( hdu ) - 1
    hdu.writeto(out_name, output_verify='exception', overwrite=True)

    if not os.path.exists(out_name):
        raise IOError('Error in refstis:functions:msjoin() -- output file not written!')

#------------------------------------------------------------------------

def crreject(input_file, workdir=None):
    if not 'oref' in os.environ:
        os.environ['oref'] = '/grp/hst/cdbs/oref/'

    path, name = os.path.split(input_file)
    name, ext = os.path.splitext(name)
    trailerfile = os.path.join(path, name+'_crreject_log.txt')

    output_blev = input_file.replace('.fits','_blev.fits')
    output_crj = input_file.replace('.fits','_crj.fits')

    with pyfits.open(input_file) as hdu:
        nimset = hdu[0].header['nextend'] / 3
        nrptexp = hdu[0].header['nrptexp']
        crcorr = hdu[0].header['crcorr']
        blevcorr = hdu[0].header['blevcorr']

    if (nimset <= 1 and crcorr != "COMPLETE"):
        print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print("Bye now... better luck next time!")
        raise ValueError('nimset <=1 and CRCORR not complete')

    if (crcorr != "COMPLETE"):
        if (nrptexp != nimset):
            pyfits.setval(input_file,'NRPTEXP',value=nimset)
            pyfits.setval(input_file,'CRSPLIT',value=1)

        pyfits.setval(input_file, 'CRCORR', value='PERFORM')
        #pyfits.setval(input_file, 'DQICORR', value='PERFORM')
        pyfits.setval(input_file, 'APERTURE', value='50CCD')
        pyfits.setval(input_file, 'APER_FOV', value='50x50')
        if (blevcorr != 'COMPLETE') :
            print('Performing BLEVCORR')
            pyfits.setval(input_file, 'BLEVCORR', value='PERFORM')
            status = basic2d(input_file,
                             output_blev,
                             outblev='',
                             dqicorr='perform',
                             blevcorr='perform',
                             doppcorr='omit',
                             lorscorr='omit',
                             glincorr='omit',
                             lflgcorr='omit',
                             biascorr='omit',
                             darkcorr='omit',
                             flatcorr='omit',
                             photcorr='omit',
                             statflag=False,
                             verbose=False,
                             trailer=trailerfile)
            if status != 0:
                try:
                    print()
                    with open(trailerfile) as tr:
                        for line in tr.readlines():
                            print('    {}'.format(line.strip()))
                finally:
                    raise Exception('BASIC2D failed to properly reduce {}'.format(input_file))
        else:
            print('Blevcorr already Performed')
            shutil.copy(input_file,output_blev)

        print('Performing OCRREJECT')
        status = ocrreject(input=output_blev,
                           output=output_crj,
                           verbose=False,
                           trailer=trailerfile)
        if status != 0:
            try:
                print()
                with open(trailerfile) as tr:
                    for line in tr.readlines():
                        print('    {}'.format(line.strip()))
            finally:
                raise Exception('OCRREJECT failed to properly reduce {}'.format(output_blev))

    elif (crcorr == "COMPLETE"):
        print("CR rejection already done")
        os.rename(input_file, output_crj)

    pyfits.setval(output_crj, 'FILENAME', value=output_crj)

    with pyfits.open(output_crj) as hdu:
        gain = hdu[0].header['atodgain']
        ccdgain = hdu[0].header['ccdgain']
        xsize = hdu[1].header['naxis1']
        ysize = hdu[1].header['naxis2']
        xbin = hdu[0].header['binaxis1']
        ybin = hdu[0].header['binaxis2']

        try:
            ncombine = hdu[0].header['ncombine']
        except:
            ncombine = hdu[1].header['ncombine']

    print(('Number of combined imsets is '+str(ncombine)+' while number of imsets is '+str(nimset ) ))
    print(('Dividing cosmic-ray-rejected image by '+str(ncombine)+'...'))
    out_div = output_crj.replace('.fits','_div.fits')

    ###this used to be a call to MSARITH, is anything else needed?
    ###modifying the error too (done), etc?
    hdu = pyfits.open(output_crj)
    hdu[('sci', 1)].data /= ncombine
    hdu[('err', 1)].data /= ncombine
    hdu.writeto(out_div, output_verify='exception', overwrite=True)

    os.remove(output_blev)
    os.remove(output_crj)

    return out_div

#------------------------------------------------------------------------

def count_imsets(file_list):
    """Count the total number of imsets in a file list.

    The total number of imsets is counted by dividing the total number
    of extensions (NEXTEND) by 3 (SCI, ERR, DQ).

    Parameters
    ----------
    file_list : list
        list of all files to be counted

    Returns
    -------
    total : int
        number of imsets

    """

    if not isinstance(file_list, list):
        file_list = [file_list]

    total = 0
    for item in file_list:
        total += pyfits.getval(item,'NEXTEND',ext=0) / 3

    return total

#------------------------------------------------------------------------

def get_keyword(file_list,keyword,ext=0):
    """ return the value from a header keyword over a list of files

    if the value is not consistent accross the input files, an assertion error
    will be raised

    """

    kw_set = set([pyfits.getval(item,keyword,ext=ext) for item in file_list])
    assert len(kw_set) == 1,' multiple values found for kw: % s'% (keyword)

    return list(kw_set)[0]

#------------------------------------------------------------------------
def get_anneal_month_dates(data_begin, data_end, database_path):
    '''
    This function uses the anneal database to get the dates of the anneal

    This is written under the assumption that data_begin and data_end fall
    in the same anneal period

    data_begin and data_end are the start and end dates of the data and should
    be in mjd
    '''
    assert data_begin > 50000, 'data_begin should be in mjd'
    assert data_end > 50000, 'data_end should be in mjd'

    db = sqlite3.connect(os.path.join(database_path, 'anneal_info.db'))
    c = db.cursor()
    c.execute("""SELECT DISTINCT start, end FROM anneals""")
    rows = [row for row in c]
    anneal_start_date = np.array([row[0] for row in rows])
    anneal_end_date = np.array([row[1] for row in rows])

    start_indx = np.where(data_begin - anneal_end_date > 0)
    anneal_period_start_indx = start_indx[0][-1] #want to start an anneal month at the end of the anneal
    end_indx = np.where(anneal_start_date - data_end > 0)
    anneal_period_end_indx = end_indx[0][0] #want to end an anneal month at the start of the anneal

    anneal_month_start = Time(anneal_end_date[anneal_period_start_indx], format = 'mjd', scale = 'utc')
    anneal_month_end = Time(anneal_start_date[anneal_period_end_indx], format = 'mjd', scale = 'utc')
    assert anneal_period_end_indx - anneal_period_start_indx == 1, \
        'data cannot cross anneals, data date range [%f - %f], anneal month [%f - %f]' %(data_begin, data_end, anneal_month_start.val, anneal_month_end.val)
    return anneal_month_start, anneal_month_end

#------------------------------------------------------------------------

def divide_anneal_month(data_begin, data_end, database_path, N_period):
    '''
    This function divides an anneal month into anneal weeks and returns
    tuples of the start and end dates of the anneal weeks
    '''

    anneal_month_start, anneal_month_end = get_anneal_month_dates(data_begin, data_end, database_path)
    total_num_days = anneal_month_end.val - anneal_month_start.val #these are in mjd so you get # of days
    base_num_days = math.floor(total_num_days / N_period)
    #For remaining days, add 1 to each week until all days are gone
    remaining_days = math.floor(total_num_days - base_num_days*N_period)
    remaining_time = total_num_days - base_num_days*N_period - remaining_days
    num_days_per_period = np.array([base_num_days for i in range(N_period)])
    num_days_per_period[-1] = num_days_per_period[-1] + remaining_time
    for i, day in enumerate(range(int(remaining_days))):
        num_days_per_period[i] += 1
    anneal_weeks = []
    wk_end = anneal_month_start.val
    for period in num_days_per_period:
        wk_start = wk_end
        wk_end = wk_start + period
        anneal_weeks.append((wk_start, wk_end))
    return anneal_weeks

#------------------------------------------------------------------------

def figure_number_of_periods(number_of_days, mode) :
    """ Determines the number of periods ('weeks') that the anneal
    'month' should be split into.

    Takes the number of days in the period and the mode ('WK' or 'BIWK')
    and returns the total number of periods.

    Parameters
    ----------
    number_of_days : int
        total number of days in the 'month'
    mode : str
        wk or biwk

    Returns
    -------
    number_of_periods : int
        the total number of periods to split the month into

    """

    #print "periods(", number_of_days,", ",mode,") called ..."
    nm = 'periods'
    msg  = "called w/ (number_of_days="+str(number_of_days)
    msg += ", mode="+mode+") "

    # set upper limits for period lengths
    MIN_DAYS_IN_WK = 6
    MAX_DAYS_IN_WK = 9
    MIN_DAYS_IN_BIWK = 11
    MAX_DAYS_IN_BIWK = 22

    DELTA_WK   =  9
    DELTA_BIWK = 18
    PERIOD_NUMBER_START = 1
    # Boundary condition catcher
    number_of_periods = -1

    # WK mode
    if mode == "WK":
        for n in range(2, number_of_days) :
            if  number_of_days < 12 :
                if (number_of_days >  MAX_DAYS_IN_WK or
                    number_of_days <  MIN_DAYS_IN_WK) :
                    # raises an ERROR
                    number_of_periods = -1
                    break
                else :
                    number_of_periods = PERIOD_NUMBER_START
                    break

            elif number_of_days <= DELTA_WK * n :
                number_of_periods = PERIOD_NUMBER_START + ( n - 1 )
                # For easier comparison to periods.py
                if number_of_days > 72 :
                    msg = "length of anneal month = " +str(number_of_days)+ "; For real? "

                    # Do we really need to set this to 0?  What's wrong with a bigger number?
                    msg  = "I give it "+str(number_of_periods)+" (WKs). "
                    msg += "Is this bad? "
                    msg += "Code would usually give such a large difference "
                    msg += "number_of_periods = 0, but why?  "

                break

    # BIWK mode
    elif mode == 'BIWK':
        for n in range(2, number_of_days) :
            if number_of_days < MAX_DAYS_IN_BIWK :
                # raises an ERROR if < MIN_DAYS_IN_BIWK1
                if number_of_days < MIN_DAYS_IN_BIWK : number_of_periods = -1
                else : number_of_periods = PERIOD_NUMBER_START
                break

            elif number_of_days <= DELTA_BIWK * n :
                number_of_periods = PERIOD_NUMBER_START + ( n - 1 )
                # For easier comparison to periods.py
                if number_of_days > 72 :
                    msg = "length of anneal month = " +str(number_of_days)+ "; For real? "
                    print(("E", msg, nm))
                    # Do we really need to set this to 0?  What's wrong with a bigger number?
                    msg  = "I give it "+str(number_of_periods)+" (BIWKs). "
                    msg += "Is this bad? "
                    msg += "Code would usually give such a large difference "
                    msg += "number_of_periods = 0, but why?  "
                    print(("W", msg, nm))
                break
    else:
        sys.exit('what mode did you put in?')

    # return value
    #print "return number_of_periods =",  number_of_periods
    return number_of_periods

#------------------------------------------------------------------------

def figure_days_in_period(N_periods, N_days, add_remainder=False):
    """Spreads out the extra days among the periods.

    Notes
    -----
    Extra days will be added to the periods beginning with the first, until there
    are no more, so the lengths may not be perfectly even but will not differ by more
    than a day.

    Parameters
    ----------
    N_periods : int
        total number of periods to be split into
    N_days : float, int
        total number of days to be divided

    Returns
    -------
    periods : list
        list of days/period: e.g. [8,8,8,7]

    Examples
    --------
    >>> figure_days_in_period(4, 28)
    [7, 7, 7, 7]

    >>> figure_days_in_period(4, 29)
    [8, 7, 7, 7]

    >>> figure_days_in_period(3, 21)
    [7, 7, 7]

    """

    remainder = N_days - int(N_days)
    N_days = int(N_days)

    base_length = N_days // N_periods
    N_extra_days = N_days - N_periods * base_length

    period_lengths = [base_length for item in range(N_periods)]

    for i in range(N_extra_days):
        index = i % N_periods
        period_lengths[index] += 1

    assert (sum(period_lengths)) == N_days, 'ERROR: extra days not spread around correctly'

    if add_remainder:
        #period_lenghts[-1] += remainder
        period_lengths[-1] += remainder

    return period_lengths

#------------------------------------------------------------------------

def translate_date_string(input_string):
    month_dict = {'Jan': 1,
                  'Feb': 2,
                  'Mar': 3,
                  'Apr': 4,
                  'May': 5,
                  'Jun': 6,
                  'Jul': 7,
                  'Aug': 8,
                  'Sep': 9,
                  'Oct': 10,
                  'Nov': 11,
                  'Dec': 12}


    date_list = input_string.split()
    time_list = date_list[3].split(':')

    month = month_dict[date_list[0]]
    day = int(date_list[1])
    year = int(date_list[2])

    hour = int(time_list[0])
    minute = int(time_list[1])
    second = int(time_list[2])

    #if time_list[3].endswith('PM'):
    #    hour += 12

    a = (14-month)/12
    y = year + 4800 - a
    m = month + 12*a - 3

    JDN = day + (153 * m + 2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
    JD = JDN + (hour - 12)/24.0 + minute/1440.0 + second/86400.0
    MJD = JD - 2400000.5
    return MJD

#------------------------------------------------------------------------

def bd_crreject(joinedfile):
    """ Check if cosmic-ray rejection has been performed on input file

    if cosmic-ray rejection has already been done on the input bias image,
    skip all calstis-related calibration steps


    """

    print(joinedfile)

    fd = pyfits.open(joinedfile)
    nimset   = fd[0].header['nextend'] / 3
    nrptexp  = fd[0].header['nrptexp']
    crcorr   = fd[0].header['crcorr']
    crdone = 0

    if (crcorr == "COMPLETE") :
        crdone = 1
        print('OK, CR rejection already done')
        os.rename(joinedfile, joinedfile.replace('_joined', '_crj') )
    else:
        print(('crcorr found = ' + crcorr))

    if (nimset <= 1 and not crdone):
        print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print("Bye now... better luck next time!")
        raise ValueError( 'Something bad happened' )

    if not crdone:
        print('FYI: CR rejection not already done')
        print(('Keyword NRPTEXP = ' + str(nrptexp) + ' while nr. of imsets = ' + str(nimset)))
        if (nrptexp != nimset):
            pyfits.setval( joinedfile,'NRPTEXP',value=nimset)
            pyfits.setval( joinedfile,'CRSPLIT',value=1)

            print(('>>>> Updated keyword NRPTEXP to '+str(nimset) ))
            print('    (and set keyword CRSPLIT to 1)' )
            print(('     in ' + joinedfile ))

    return crdone

#------------------------------------------------------------------------

def bd_calstis(joinedfile, thebiasfile=None):
    """ Run CalSTIS on the joined file

    Header keywords will be set for ocrreject to work correctly and not
    flag regions outside the original aperture:
    APERTURE --> 50CCD
    APER_FOV --> '50x50'
    DARKCORR --> 'OMIT'
    FLATCORR --> 'OMIT'

    Parameters
    ----------
    joinedfile : str
        join of multiple input darks
    thebiasfile : str, bool
        the biasfile to be subtracted by basic2d

    """

    with pyfits.open(joinedfile, 'update') as hdu:
        hdu[0].header['CRCORR'] = 'PERFORM'
        hdu[0].header['APERTURE'] = '50CCD'
        hdu[0].header['APER_FOV'] = '50x50'
        hdu[0].header['DARKCORR'] = 'OMIT'
        hdu[0].header['FLATCORR'] = 'OMIT'

        if thebiasfile:
            hdu[0].header['BIASFILE'] = thebiasfile

    crj_file = joinedfile.replace('.fits', '_crj.fits')

    path, name = os.path.split(joinedfile)
    name, ext = os.path.splitext(name)
    trailerfile = os.path.join(path, name+'_bd_calstis_log.txt')

    if os.path.exists(crj_file):
        print('Deleting old file: %s' % crj_file)
        os.remove(crj_file)

    if os.path.exists(trailerfile):
        print('Deleting old file: %s' % trailerfile)
        os.remove(trailerfile)

    print('Running CalSTIS on %s' % joinedfile)
    print('to create: %s' % crj_file)
    status = calstis(joinedfile,
                     wavecal="",
                     outroot="",
                     savetmp=False,
                     verbose=False,
                     trailer=trailerfile)
    if status != 0:
        try:
            print()
            with open(trailerfile) as tr:
                for line in tr.readlines():
                    print('    {}'.format(line.strip()))
        finally:
            raise Exception('CalSTIS failed to properly reduce {}'.format(joinedfile))

    pyfits.setval(crj_file, 'FILENAME', value=os.path.split(crj_file)[1])

#------------------------------------------------------------------------

def RemoveIfThere(item):
    """Remove a file only if it already exists

    Parameters
    ----------
    item : str
        file to be removed

    Examples
    --------
    >>> RemoveIfThere('/path/to/file.txt')

    """

    if os.path.exists(item):
        os.remove(item)

#------------------------------------------------------------------------

def refaver(reffiles, combined_name):
    """Average two reference files together using itools msarith.

    Parameters
    ----------
    reffiles : list
        List of reference files to be averaged together
    combined_name : str
        Output name of the combined file

    """
    from .msarith import msarith

    print('#-----------------------#')
    print('combining datasets')
    print(reffiles)
    print('into')
    print(combined_name)
    print('#-----------------------#')


    all_paths = {os.path.split(item)[0] for item in reffiles}
    assert len(all_paths) == 1, "More than one path found"

    initial_dir = os.getcwd()
    if list(all_paths)[0] != '':
        os.chdir(list(all_paths)[0])

    all_subfiles = []
    for subfile in reffiles:
        subfile = os.path.split(subfile)[-1]
        outfile = subfile.replace('.fits', '_aver.fits')
        print("Running (itools) msarith / 2 on {}".format(subfile))
        msarith(subfile, '/', 2, outfile, verbose=True)
        all_subfiles.append(outfile)

    assert len(all_subfiles) == 2, 'Length of subfiles doesnt equal 2: {}'.format(all_subfiles)

    if not combined_name.endswith('.fits'):
        combined_name = combined_name + '.fits'

    #-- remove path from output name
    combined_name = os.path.split(combined_name)[-1]

    msarith(all_subfiles[0], '+', all_subfiles[1], combined_name, verbose=True)

    for filename in all_subfiles:
        os.remove(filename)

    #-- move back to beginning location
    os.chdir(initial_dir)

#------------------------------------------------------------------------

def apply_dark_correction(filename, expstart):
    """Perform temperature scaling to input dark file

    All science extensions in the input filename will be scaled to the reference
    temperatue of 18.0 c.

    Parameters
    ----------
    filename : str
        full path to input FITS file
    expstart : str
        start time in MJD of the dataset

    """

    dark_v_temp = 0.07
    s2ref_temp = 18.0
    with pyfits.open(filename, mode = 'update') as ofile:
        if 'tempcorr' not in ofile[0].header:
            nextend = ofile[0].header['nextend']

            for ext in np.arange(1, nextend, 3):
                occdhtav = ofile[ext].header['OCCDHTAV']
                factor = 1.0 / (1.0 + dark_v_temp * (float(occdhtav) - s2ref_temp))
                ofile[ext].data = ofile[ext].data * factor
                print('{}, ext {}: Scaling data by '.format(filename, ext), factor, ' for temperature: ', occdhtav)
                ofile[ext+1].data = np.sqrt((ofile[ext+1].data)**2 * (factor**2)) #Modify the error array
                ofile[ext].header.add_history('File scaled for Side-2 temperature uncertainty by data * (1.0 + %f * (%f - %f)) following description is STIS TIR 2004-01' %(dark_v_temp, occdhtav, s2ref_temp))

            ofile[0].header['tempcorr'] = 'COMPLETE'
        else:
            print('TEMPCORR = %s, no temperature correction applied to %s' %(ofile[0].header['tempcorr'], filename))

#-------------------------------------------------------------------------------

def bias_subtract_data(filename, biasfile, outdir=''):
    """Perform bias subtraction on input dataset

    basic2d from calstis will be run on the input dataset with the
    steps DQICORR, BLEVCORR, and BIASCORR set to perform.

    Parameters
    ----------
    filename : str
        full path to input FITS file
    biasfile : str
        full path to the bias FITS file to be subtracted
    outdir : str, optional
        if specified, the directory to send the bias-subtracted output file

    Returns
    -------
    filename : str
        full_path to the bias subtracted file

    """

    with pyfits.open(filename) as hdu:
        if (hdu[0].header['BLEVCORR'] == 'COMPLETE') or (hdu[0].header['BIASCORR'] == 'COMPLETE'):
            print("BIAS correction already done for {}".format(filename))
            return filename

        # Did we run the pixel-based CTE correction on the input file?
        # If so, this changes the extension of the output file.
        cte_corrected = hdu[0].header.get('PCTECORR', 'OMIT').strip() == 'COMPLETE'

    path, name = os.path.split(filename)
    name, ext = os.path.splitext(name)
    if outdir:
        path = outdir
    trailerfile = os.path.join(path, name + '_bias_subtract_log.txt')

    biasfile = make_path_safe(biasfile)

    output_filename = name.rsplit('_raw', 1)[0] + ('_flc' if cte_corrected else '_flt') + '.fits'
    output_filename = os.path.join(path, output_filename)
    if os.access(output_filename, os.F_OK | os.W_OK):
        os.remove(output_filename)

    # Check to see if input directory allows write access.  Copy data to output directory if not.
    if (not os.access(filename, os.W_OK)) and (os.path.dirname(filename) != os.path.dirname(output_filename)):
        shutil.copy(filename, os.path.dirname(output_filename))
        filename = os.path.join(os.path.dirname(output_filename), os.path.basename(filename))
        # Inherit destination directory read/write permissions (not execute) + user write:
        all_read_write = stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH | stat.S_IWOTH
        os.chmod(filename, (os.stat(os.path.dirname(filename)).st_mode | stat.S_IWUSR) & all_read_write)
    elif not os.access(filename, os.W_OK):
        raise IOError('Directory with input data does not allow write access. ' + \
                      f'Please specify an explicit outdir.  {filename}')

    pyfits.setval(filename, 'BIASFILE', ext=0, value=biasfile, comment='')
    status = basic2d(filename,
                     output=output_filename,
                     dqicorr='perform',
                     blevcorr='perform',
                     biascorr='perform',
                     doppcorr='omit',
                     lorscorr='omit',
                     glincorr='omit',
                     lflgcorr='omit',
                     darkcorr='omit',
                     flatcorr='omit',
                     photcorr='omit',
                     verbose=False,
                     trailer=trailerfile)
    if status != 0:
        try:
            print()
            with open(trailerfile) as tr:
                for line in tr.readlines():
                    print('    {}'.format(line.strip()))
        finally:
            raise Exception(f'BASIC2D failed to properly reduce {filename}')

    return output_filename

#-------------------------------------------------------------------------------

def make_path_safe(filename):
    """Make a full path to file safe for use in FITS headers.

    For full paths that are less than 67 characters, the filename is simply
    returned.  When the full path is equal to or greater than 67, then the path
    is inserted into the environment as 'refdir' and the filename is returned as
    'refdir$filename'.  This will prevent the filename from being split accross
    multiple FITS keywords, and is a convention understood by many
    tasks/pipelines.

    Parameters
    ----------
    filename : str
        Full path + name to the file.

    Returns
    -------
    filename : str
        Safe filename that can fit in a single FITS keyword.

    Examples
    --------
    >>> make_path_safe('/short/path/reference_file.fits')
    /short/path/reference_file.fits

    >>> make_path_safe('/really/really/really/really/really/really/really/really/reference_file.fits')
    refdir$reference_file.fits

    """
    if len(filename) < 67:
        return filename

    path, filename = os.path.split(filename)
    #-- Calstis wants the '/' at the end
    if not path.endswith('/'):
        path += '/'

    os.environ['refdir'] = path

    filename = 'refdir$'+filename

    return filename

#-------------------------------------------------------------------------------
