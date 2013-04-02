
def crreject( input_file, workdir=None) :
    import pyfits
    import os
    from pyraf import iraf
    from iraf import stsdas,hst_calib,stis
    from pyraf.irafglobals import *

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
        print('E',"Sorry, your input image seems to have only 1 imset, but it isn't cr-rejected.")
        print('E',"This task can only handle 'raw' or 'flt images with the NEXTEND keyword equal to 3*N (N > 1).")
        print('E',"Bye now... better luck next time!")


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
           iraf.basic2d(input_file, output_blev,
                        outblev = '', dqicorr = 'perform', atodcorr = 'omit',
                        blevcorr = 'perform', doppcorr = 'omit', lorscorr = 'omit',
                        glincorr = 'omit', lflgcorr = 'omit', biascorr = 'omit',
                        darkcorr = 'omit', flatcorr = 'omit', shadcorr = 'omit',
                        photcorr = 'omit', statflag = no, verb=no)
       else:
           print('Blevcorr alread Performed')
           shutil.copy(input_file,output_blev)

       print('Performing OCRREJECT')
       iraf.ocrreject(input=output_blev, output=output_crj, verb=no)

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
    print('HEADER INFO: CCD gain : '+str(gain)+' electrons/ADU' )
    print('             BINAXIS1 : '+str(xbin ) )
    print('             BINAXIS2 : '+str(ybin ) )
    print('Dividing cosmic-ray-rejected image by '+str(ncombine)+'...')
    out_div = output_crj.replace('.fits','_div.fits')
    iraf.msarith( output_crj, '/', ncombine, out_div, verbose = 0)
    
    return out_div
    #return tmpsuper, xbin, ybin, ccdgain, gain, xsize, ysize, ncombine

#------------------------------------------------------------------------------------

def split_images( imglist,outname='' ):
    from pyraf import iraf
    from iraf import stsdas,toolbox,imgtools,mstools
    print 'Splitting images'
    for dataset in imglist:
        print dataset
        iraf.mssplit(inimg=dataset, outimg=outname, extension = '*', retain='no', 
                     Stderr='dev$null'
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
    else : 
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
                    PrintMsg("E", msg, nm)
                    # Do we really need to set this to 0?  What's wrong with a bigger number?
                    msg  = "I give it "+str(number_of_periods)+" (BIWKs). "
                    msg += "Is this bad? "
                    msg += "Code would usually give such a large difference "
                    msg += "number_of_periods = 0, but why?  "
                    PrintMsg("W", msg, nm)
                break


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
    base_length = N_days//N_periods
    N_extra_days = N_days - N_periods * base_length

    period_lengths = [ base_length for item in xrange(N_periods) ]

    for i in range(N_extra_days):
        index = i%N_periods
        period_lengths[index] += 1

    assert (sum(period_lengths)) == N_days,'ERROR: extra days not spread around correctly'

    return period_lengths

#-------------------------------------------------------------------------------

def translate_date_string(input_string):
    month_dict = {'Jan':0, 'Feb':1, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6,
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

    #print 'Converting %d/%d/%d %d:%d:%d to MJD'%(month,day,year,hour,minute,second)

    a = (14-month)/12
    y = year + 4800 - a
    m = month + 12*a - 3

    JDN = day + (153*m + 2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
    JD = JDN + (hour-12)/24.0 + minute/1440.0 + second/86400.0
    MJD = JD - 2400000.5

    return MJD

#-------------------------------------------------------------------------------
