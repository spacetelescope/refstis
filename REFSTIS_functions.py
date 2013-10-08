import pyfits
import numpy as np
import os
from scipy.signal import medfilt
from astropy.time import Time
import pdb
import support

#-------------------------------------------------------------------------------

def iterclip( in_array, sigma, maxiter ):
    skpix = indata.reshape( indata.size, )
 
    ct = indata.size
    iter = 0; c1 = 1.0 ; c2 = 0.0
 
    while (c1 >= c2) and (iter < maxiter):
        lastct = ct
        medval = numpy.median(skpix)
        sig = numpy.std(skpix)
        wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]
 
        c1 = abs(ct - lastct)
        c2 = converge_num * lastct
        iter += 1

#---------------------------------------------------------------------------

def update_header_from_input( filename, input_list ):
    """  Updates header of output file using keywords from the input data
    
    If a header keyword is not consistent in this step, an error will be 
    raised

    """

    gain = get_keyword( input_list, 'CCDGAIN', 0)
    xbin = get_keyword( input_list, 'BINAXIS1', 0)
    ybin = get_keyword( input_list, 'BINAXIS2', 0)
    amp = get_keyword( input_list, 'CCDAMP', 0)
    targname = get_keyword( input_list, 'TARGNAME', 0)
    
    if targname == 'BIAS':
        filetype = 'CCD BIAS IMAGE'
    elif targname == 'DARK':
        filetype = 'DARK IMAGE'
    else:
        raise ValueError( 'targname %s not understood' % str( targname ) )

    if gain == 1:
        frequency = 'Weekly'
    elif gain == 4:
        frequency = 'Bi-Weekly'
    else:
        raise ValueError( 'Frequency %s not understood' % str( frequency ) )
    data_start, data_end = get_start_and_endtimes(input_list)
    hdu = pyfits.open( filename, mode='update' )
    hdu[0].header['FILENAME'] = os.path.split( filename )[1]
    hdu[0].header['CCDGAIN'] = gain
    hdu[0].header['BINAXIS1'] = xbin
    hdu[0].header['BINAXIS2'] = ybin  
    hdu[0].header['NEXTEND'] = 3
    hdu[0].header['PEDIGREE'] = 'INFLIGHT %s %s' %(data_start, data_end)
    hdu[0].header['USEAFTER'] = 'value'
    hdu[0].header['DESCRIP'] = "%s gain=%d %s for STIS CCD data taken after XXX" % (frequency, gain, targname.lower)
    hdu[0].header.add_comment('Reference file created by %s' % __name__ )

    while len( hdu[0].header['DESCRIP'] ) < 67:
        hdu[0].header['DESCRIP'] = hdu[0].header['DESCRIP'] + '-'

    hdu[0].header.add_history('Used input files')
    for item in input_list:
        hdu[0].header.add_history(item)

    hdu.flush()
    hdu.close()

def get_start_and_endtimes(input_list):
    times = []
    for ifile in input_list:
        times.append(pyfits.getval(ifile, 'texpstrt', 0))
        times.append(pyfits.getval(ifile, 'texpend', 0))
    times.sort()
    times = Time(times, format = 'mjd', scale = 'utc', out_subfmt = 'date').iso
    start_list = times[0].split('-')
    end_list = times[-1].split('-')
    
    #return strings in format mm/dd/yyyy
    start_str = '%s/%s/%s' %(start_list[1], start_list[0], start_list[2])
    end_str = '%s/%s/%s' %(end_list[1], end_list[0], end_list[2])
    return start_str, end_str

    

#---------------------------------------------------------------------------

def make_resicols_image( residual_image, yfrac=None ):    

    import numpy as np

    if yfrac:
        ystart = int( (yfrac[0] / 100.0) * residual_image.shape[0] )
        yend = int( (yfrac[1] / 100.0) * residual_image.shape[0] )
        residual_columns = np.mean( residual_image[ ystart:yend ], axis=0 )
    else:
        residual_columns = np.mean( residual_image, axis=0 )

    residual_columns_image = residual_columns.repeat( residual_image.shape[0] ).reshape( residual_image.shape )

    return residual_columns_image

#-------------------------------------------------------------------------------

def make_residual( mean_bias ):
    """Create residual image

    Median filter the median with a 15 x 3 box and subtract from the mean
    to produce the residual image.


    """
    import pyfits
    
    mean_hdu = pyfits.open( mean_bias )
    mean_image = mean_hdu[ ('sci', 1) ].data

    median_image = medfilt( mean_image, (3, 15) )

    diffmean = support.sigma_clip( mean_image, sigma=3, iterations=40)[1] - support.sigma_clip( median_image, sigma=3, iterations=40 )[1]
    residual_image = median_image - diffmean

    return residual_image, median_image

#-------------------------------------------------------------------------------

def normalize_crj( filename ):
    """ Normalize the input filename by exptim/gain and flush hdu 

    """

    import pyfits

    hdu = pyfits.open( filename, mode='update' )

    exptime = hdu[0].header[ 'TEXPTIME' ]
    gain = hdu[0].header[ 'ATODGAIN' ]

    hdu[ ('sci', 1) ].data /= (float(exptime) / gain)

    hdu[0].header['TEXPTIME'] = 1

    hdu.flush()
    hdu.close()

#-------------------------------------------------------------------------------

def msjoin( imset_list, out_name='joined_out.fits' ):
    """ Replicate msjoin functionality in pure python

    """

    import pyfits

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
    hdu.writeto( out_name )

#-------------------------------------------------------------------------------

def crreject( input_file, workdir=None) :
    import pyfits
    import os
    from stistools.basic2d import basic2d
    from stistools.ocrreject import ocrreject

    if not 'oref' in os.environ:
        os.environ['oref'] = '/grp/hst/cdbs/oref/'

    output_blev = input_file.replace('.fits','_blev.fits')
    output_crj = input_file.replace('.fits','_crj.fits')
    # 
    # between the long file paths and not being able to find EPC files, 
    # need IRAF to run in the work dir
    # 
    #os.chdir( workdir )
    #
    # if cosmic-ray rejection has already been done on the input bias image,
    # skip all calstis-related calibration steps
    #

    fd = pyfits.open( input_file )
    nimset   = fd[0].header['nextend'] / 3
    nrptexp  = fd[0].header['nrptexp']
    crcorr   = fd[0].header['crcorr']
    blevcorr = fd[0].header['blevcorr']
    del fd

    if (nimset <= 1 and crcorr != "COMPLETE"):
        print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print("Bye now... better luck next time!")
        raise ValueError( 'nimset <=1 and CRCORR not complete' )

    if (crcorr != "COMPLETE"):
        if (nrptexp != nimset):
            pyfits.setval(input_file,'NRPTEXP',value=nimset)
            pyfits.setval(input_file,'CRSPLIT',value=1)

        pyfits.setval(input_file, 'CRCORR', value='PERFORM')
        pyfits.setval(input_file, 'APERTURE', value='50CCD')
        pyfits.setval(input_file, 'APER_FOV', value='50x50')
        if (blevcorr != 'COMPLETE') :
            print('Performing BLEVCORR')
            pyfits.setval(input_file, 'BLEVCORR', value='PERFORM')
            basic2d(input_file, output_blev,
                         outblev = '', dqicorr = 'perform', 
                         blevcorr = 'perform', doppcorr = 'omit', lorscorr = 'omit',
                         glincorr = 'omit', lflgcorr = 'omit', biascorr = 'omit',
                         darkcorr = 'omit', flatcorr = 'omit', 
                         photcorr = 'omit', statflag = False, verbose = False)
        else:
            print('Blevcorr alread Performed')
            shutil.copy(input_file,output_blev)

        print('Performing OCRREJECT')
        ocrreject(input=output_blev, output=output_crj, verbose=False)

    elif (crcorr == "COMPLETE"):
        print "CR rejection already done"
        os.rename(input_file, output_crj )
  
    pyfits.setval(output_crj, 'FILENAME',
                  value=output_crj)
  
    fd = pyfits.open(output_crj)
    gain    = fd[0].header['atodgain']
    ccdgain = fd[0].header['ccdgain']
    xsize   = fd[1].header['naxis1']
    ysize   = fd[1].header['naxis2']
    xbin    = fd[0].header['binaxis1']
    ybin    = fd[0].header['binaxis2']
    try:
        ncombine = fd[0].heade['ncombine']
    except:
        ncombine = fd[1].header['ncombine']
    fd.close()
    del fd
  
    print('Number of combined imsets is '+str(ncombine)+' while number of imsets is '+str(nimset ) )
    print('Dividing cosmic-ray-rejected image by '+str(ncombine)+'...')
    out_div = output_crj.replace('.fits','_div.fits')

    #this used to be a call to MSARITH, is anything else needed?
    #modifying the error too, etc?
    hdu = pyfits.open( output_crj )
    hdu[1].data /= ncombine
    hdu.writeto( out_div )

    os.remove( output_blev )
    os.remove( output_crj )

    return out_div

#-------------------------------------------------------------------------------
                     
def count_imsets( file_list ):
    """ Count the total number of imsets in a file list by dividing the number
    of extensions (NEXTEND) by 3

    """

    import pyfits
    total = 0
    for item in file_list:
        total += pyfits.getval(item,'NEXTEND',ext=0) / 3
    return total

#-------------------------------------------------------------------------------

def get_keyword( file_list,keyword,ext=0):
    """ return the value from a header keyword over a list of files
    
    if the value is not consistent accross the input files, an assertion error
    will be raised
    
    """

    import pyfits

    kw_set = set( [pyfits.getval(item,keyword,ext=ext) for item in file_list] )
    assert len(kw_set) == 1,' multiple values found for kw: % s'% (keyword)

    return list(kw_set)[0]

#------------------------------------------------------------------------------

def figure_number_of_periods(number_of_days, mode) :
    """ Determines the number of periods ('weeks') that the anneal
    'month' should be split into.  

    Takes the number of days in the period and the mode ('WK' or 'BIWK')
    and returns the total number of periods. 
 
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
    if ( mode == "WK" ) :
        for n in xrange( 2, number_of_days) :
            if  number_of_days < 12 :
                if (number_of_days >  MAX_DAYS_IN_WK or
                    number_of_days <  MIN_DAYS_IN_WK ) :
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
    elif ( mode == 'BIWK') : 
        for n in xrange( 2, number_of_days) :
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
                    print("E", msg, nm)
                    # Do we really need to set this to 0?  What's wrong with a bigger number?
                    msg  = "I give it "+str(number_of_periods)+" (BIWKs). "
                    msg += "Is this bad? "
                    msg += "Code would usually give such a large difference "
                    msg += "number_of_periods = 0, but why?  "
                    print("W", msg, nm)
                break
    else:
        sys.exit('what mode did you put in?')

    # return value
    #print "return number_of_periods =",  number_of_periods
    return number_of_periods

#-------------------------------------------------------------------------------

def figure_days_in_period(N_periods, N_days):
    """Spreads out the extra days among the periods.

    Extra days will be added to the periods beginning with the first, until there
    are no more, so the lengths may not be perfectly even but will not differ by more
    than a day.

    Returns:
       list of days/period: e.g. [8,8,8,7]
    
    """

    base_length = N_days // N_periods
    N_extra_days = N_days - N_periods * base_length

    period_lengths = [ base_length for item in xrange(N_periods) ]

    for i in range(N_extra_days):
        index = i % N_periods
        period_lengths[index] += 1

    assert (sum(period_lengths)) == N_days, 'ERROR: extra days not spread around correctly'

    return period_lengths

#-------------------------------------------------------------------------------

def translate_date_string(input_string):
    month_dict = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6,
                  'Jul':7, 'Aug':8, 'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12}


    date_list = input_string.split()
    time_list = date_list[3].split(':')

    month = month_dict[ date_list[0] ]
    day = int( date_list[1] )
    year = int( date_list[2] )
    
    hour = int( time_list[0] )
    minute = int( time_list[1] )
    second = int( time_list[2] )

    if time_list[3].endswith('PM'): hour += 12
 
    a = (14-month)/12
    y = year + 4800 - a
    m = month + 12*a - 3

    JDN = day + (153 * m + 2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
    JD = JDN + (hour - 12)/24.0 + minute/1440.0 + second/86400.0
    MJD = JD - 2400000.5
    return MJD

#---------------------------------------------------------------------------

def bd_crreject(joinedfile):
    """ Check if cosmic-ray rejection has been performed on input file
    
    if cosmic-ray rejection has already been done on the input bias image,
    skip all calstis-related calibration steps
    

    """

    import pyfits
    import os

    print joinedfile

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
        print('crcorr found = ' + crcorr)

    if (nimset <= 1 and not crdone):
        print("Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print("This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print("Bye now... better luck next time!")
        raise ValueError( 'Something bad happened' )

    if not crdone:
        print('FYI: CR rejection not done')
        print('Keyword NRPTEXP = ' + str(nrptexp) + ' while nr. of imsets = ' + str(nimset))
        if (nrptexp != nimset):
            pyfits.setval( joinedfile,'NRPTEXP',value=nimset)
            pyfits.setval( joinedfile,'CRSPLIT',value=1)
            
            print('>>>> Updated keyword NRPTEXP to '+str(nimset) )
            print('    (and set keyword CRSPLIT to 1)' )
            print('     in ' + joinedfile )

    return crdone

#--------------------------------------------------------------------------

def bd_calstis(joinedfile, thebiasfile=None ) :
    """ Run CalSTIS on the joined file

    Header keywords will be set for ocrreject to work correctly and not 
    flag regions outside the original aperture:
    APERTURE --> 50CCD
    APER_FOV --> '50x50'

    """

    import pyfits
    from stistools.calstis import calstis
    import os
    import shutil

    pyfits.setval(joinedfile, 'CRCORR', value='PERFORM')
    pyfits.setval(joinedfile, 'APERTURE', value='50CCD')
    pyfits.setval(joinedfile, 'APER_FOV', value='50x50')
    pyfits.setval(joinedfile, 'DARKCORR', value='OMIT')
    pyfits.setval(joinedfile, 'FLATCORR', value='OMIT')

    ### This was causing a floating point error in CalSTIS
    ### Perhaps the biasfile i was using was bad?
    #if thebiasfile:
    #   pyfits.setval(joinedfile, 'BIASFILE', value=thebiasfile)

    crj_file = joinedfile.replace('.fits', '_crj.fits')

    print 'Running CalSTIS on %s' % joinedfile 
    print 'to create: %s' % crj_file
    calstis(joinedfile, wavecal="", outroot="",
                 savetmp=no, verbose=yes)
    
    pyfits.setval(crj_file, 'FILENAME', value=os.path.split(crj_file)[1] )

#--------------------------------------------------------------------------

def RemoveIfThere(item):
    """  Simple copy from OPUSUtils to remove a file if it is there

    """

    import os
    if os.path.exists(item):
        os.remove(item)

#--------------------------------------------------------------------------

def refaver( reffiles, combined_name ):
    """ Average two reference files toghether

    Arithmetic for the combination is done using msarith

    """

    from pyraf import iraf
    from iraf import mstools
    import os

    if not combined_name.endswith('.fits'):
        combined_name = combined_name + '.fits'

    print '#-----------------------#'
    print 'combining datasets'
    print reffiles
    print 'into'
    print combined_name
    print '#-----------------------#'

    all_subfiles = []
    for subfile in reffiles:
        outfile = subfile.replace('.fits', '_aver.fits')
        iraf.msarith( subfile, '/', 2, outfile, verbose=1 )
        all_subfiles.append( outfile )

    assert len(all_subfiles) == 2, 'Length of subfiles doesnt equal 2'

    iraf.msarith( all_subfiles[0], '+', all_subfiles[1], combined_name, verbose=1)

    for filename in all_subfiles:
        os.remove( filename )

#--------------------------------------------------------------------------
