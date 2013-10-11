#!/usr/bin/env python

"""Script to fill database with cos anneal month start and end times.  
"""

from __future__ import division

import sqlite3
import glob
import os
import sys
import pyfits
import argparse
import textwrap
import support
import time
import shutil
import re
import numpy as np
from support import SybaseInterface
from support import createXmlFile, submitXmlFile

from REFSTIS_functions import figure_days_in_period, figure_number_of_periods, translate_date_string
import REFSTIS_pop_db
import REFSTIS_basedark
import REFSTIS_weekdark
import REFSTIS_refbias
import REFSTIS_weekbias
import REFSTIS_basejoint
import REFSTIS_functions

products_directory = '/user/ely/STIS/refstis/darks_biases/'
retrieve_directory = '/user/ely/STIS/refstis/requested/'

#dark_proposals = [7600, 7601, 8408, 8437, 8837, 8864, 8901, 8902, 9605, 9606, 
#                  10017, 10018, 11844, 11845, 12401, 12402, 12741, 12742]
#bias_proposals = [7600, 7601, 8409, 8439, 8838, 8865, 8903, 8904, 9607, 9608, 
#                  10019, 10020, 11846, 11847, 12403, 12404, 12743, 12744]

dark_proposals = [ 11844, 11845, 12400, 12401, 12741, 12742, 13131, 13132]
bias_proposals = [ 11846, 11847, 12402, 12403, 12743, 12744, 13133, 13134]


#-------------------------------------------------------------------------------

def get_new_periods():
    print '#-------------------#'
    print 'Reading from database'
    print '#-------------------#\n'
    db = sqlite3.connect( "anneal_info.db" )
    c = db.cursor()
    table = 'anneals'

    c.execute("""SELECT * FROM %s """ % (table))

    all_info = [row for row in c]

    table_id_all = [row[0] for row in all_info]
    proposal_id_all = [row[1] for row in all_info]
    visit_id_all = [ int(row[2]) for row in all_info ]
    anneal_start_all = [row[3] for row in all_info]
    anneal_end_all = [row[4] for row in all_info]

    dirs_to_process = []

    for i in range( len(table_id_all) )[::-1]:
        if i == len(table_id_all) - 1: continue
        ref_begin = anneal_end_all[i]
        ref_end = anneal_start_all[i + 1]
        proposal = proposal_id_all[i + 1]  # defined from proposal of the next anneal
        visit = visit_id_all[i + 1]
        year, month, day, dec_year = support.mjd_to_greg(ref_begin)
        end_year, end_month, end_day, dec_year = support.mjd_to_greg(ref_end)       
 
        if visit < 10:
            visit = '0'+str(visit)
        else:
            visit = str(visit)

        print '#--------------------------------#'
        print 'Searching for new observations for'
        print 'Period: %d_%d_%s'%(year, proposal, visit)  
        print 'MJD %5.5f %5.5f'%(ref_begin, ref_end)
        print month, day, year, ' to ', end_month, end_day, end_year
        print '#--------------------------------#'

        products_folder = os.path.join( products_directory, '%d_%d_%s' % (year, proposal, visit) )
        dirs_to_process.append( products_folder )

        if not os.path.exists( products_folder ): 
            os.mkdir( products_folder )
        
        already_retrieved = []
        for root, dirs, files in os.walk( products_folder ):
            for filename in files:
                if filename.endswith('_raw.fits'):
                    already_retrieved.append( filename[:9].upper() )

        new_obs = get_new_obs('DARK', ref_begin, ref_end) + get_new_obs('BIAS', ref_begin, ref_end)
        obs_to_get = [ obs for obs in new_obs if not obs in already_retrieved ]

        if not len( obs_to_get ): 
            print 'No new obs to get, skipping this period\n\n'
            continue
        else: 
            print 'Found new observations for this period'
            print obs_to_get, '\n\n'

        response = collect_new( obs_to_get )
        move_obs( obs_to_get, products_folder) 
        separate_obs( products_folder, ref_begin, ref_end )

    return dirs_to_process

#-------------------------------------------------------------------------------

def split_files(  all_files ):
    all_info = [ ( pyfits.getval(filename,'EXPSTART', 1), filename ) for filename in all_files ]
    all_info.sort()
    all_files = [ line[1] for line in all_info ]

    halfway = len(all_files)//2

    super_list = [ all_files[:halfway],
                   all_files[halfway:] ]

    return super_list

#-------------------------------------------------------------------------------

def pull_out_subfolders( root_folder ):
    """ Walk through input folder and use regular expressions to return lists of
    the folders for each gain and for each week

    parameters
    ----------
    root_folder
        string, the folder to walk through

    returns
    -------
    gain_folders
        list, containing folders for each gain
        
    week_folders
        list, containing folders for each week

    """

    gain_folders = []
    week_folders = []
    for root, dirs, files in os.walk( root_folder ):
        tail = os.path.split( root )[-1]

        if 'wk' in tail:     
            week_folders.append( root )
        
        # ex: 1-1x1 or 4-1x1
        if re.search( '([0-4]-[0-4]x[0-4])', tail):
            gain_folders.append( root )

    return gain_folders, week_folders

#-------------------------------------------------------------------------------

def pull_info( foldername ):
    """ Pull proposal and week number from folder name 

    A valid proposal is a string of 5 numbers from 0-9
    A valid week is a string of 'wk' + 2 numbers ranging from 0-9.
    A valid biweek is the string 'biwk' + 2 numbers from 0-9.
    
    parameters
    ----------
    foldername
        string, name of folder to search in

    returns
    -------
    proposal
        string, proposal number as string

    week
        string, week of the anneal (wk01, etc)
    """
    try:
        proposal = re.search('(_[0-9]{5}_)', foldername ).group().strip('_')
    except:
        proposal = ''
    try:
        week = re.search('([bi]*wk0[0-9])', foldername ).group()
    except:
        week = ''
    return proposal, week

#-------------------------------------------------------------------------------

def make_ref_files( root_folder, clean=False ):
    """ Make all refrence files for a given folder

    This functions is very specific to the REFSTIS pipeline, and requires files
    and folders to have certain naming conventions.  

    """

    print '#-----------------------------#'
    print '#  Making all ref files for   #'
    print  root_folder
    print '#-----------------------------#'

    if not os.path.exists( root_folder ): raise IOError( 'Root folder does not exist' )

    if clean:  clean_directory( root_folder )

    bias_threshold = { (1, 1, 1):98,  (1, 1, 2):25,  (1, 2, 1):25,  (1, 2, 2):7, 
                       (1, 4, 1):7,  (1, 4, 2):4,  (4, 1, 1):1 }

    gain_folders, week_folders = pull_out_subfolders( root_folder )
        

    ######################
    # make the base biases
    ######################
    if os.path.exists(os.path.join(root_folder, 'biases')):
        for folder in gain_folders:
            all_dir = os.path.join( folder, 'all' )
            if not os.path.exists( all_dir ):  os.mkdir( all_dir )

            for root, dirs, files in os.walk( folder ):
                if root.endswith('all'): continue
                for filename in files:
                    if filename.endswith('_raw.fits'):
                        shutil.copy( os.path.join( root, filename), all_dir )

            all_files = glob.glob( os.path.join( all_dir, '*_raw.fits') )
            basebias_name = os.path.join( all_dir, 'basebias.fits' )
            if not os.path.exists( basebias_name ):
                REFSTIS_basejoint.make_basebias( all_files , basebias_name )
            else:
                print 'Basebias already created, skipping'
    else:
        print 'no folder %s exists, not making a basebias' %(os.path.join(root_folder, 'biases'))
    ######################
    # make the base darks
    ######################
    if os.path.exists(os.path.join(root_folder, 'darks')):
        dark_folder = os.path.join( root_folder, 'darks' )
        all_dir = os.path.join( dark_folder, 'all' )
        if not os.path.exists( all_dir ):  os.mkdir( all_dir )

        for root, dirs, files in os.walk( dark_folder ):
            if root.endswith('all'): continue
            for filename in files:
                if filename.endswith('_raw.fits'):
                    shutil.copy( os.path.join( root, filename), all_dir )

        all_files = glob.glob( os.path.join( all_dir, '*_raw.fits') )

        basebias_name = os.path.join( root_folder, 'biases/1-1x1/all/', 'basebias.fits' )
        basedark_name = os.path.join( all_dir, 'basedark.fits' )
        if not os.path.exists( basedark_name ):
            REFSTIS_basedark.make_basedark( all_files , basedark_name, basebias_name )
        else:
            print 'Basedark already created, skipping'
    else:
        print 'no folder %s exists, not making a basedark' %(os.path.join(root_folder, 'darks'))
    
    ####################
    # make the weekly biases and darks
    ####################

    for folder in week_folders:
        REFBIAS = False
        WEEKBIAS = False
 
        BASEDARK = False
        WEEKDARK = False
        print 'Processing %s'%(folder)
        
        proposal, wk = pull_info( folder )

        raw_files = glob.glob( os.path.join( folder, '*raw.fits') )
        n_imsets = REFSTIS_functions.count_imsets( raw_files )

        #if n_imsets > 140: 
        #    sys.exit('error, too many imsets found: %d'%(n_imsets) )
        
        gain = REFSTIS_functions.get_keyword( raw_files, 'CCDGAIN', 0)
        xbin = REFSTIS_functions.get_keyword( raw_files, 'BINAXIS1', 0)
        ybin = REFSTIS_functions.get_keyword( raw_files, 'BINAXIS2', 0)
        
        if re.search('/biases/', folder):
            filetype = 'bias'
            REFBIAS = True

            if n_imsets < bias_threshold[ (gain, xbin, ybin) ]:
                WEEKBIAS = True

        elif re.search('/darks/', folder):
            filetype = 'dark'
            BASEDARK = True
            WEEKDARK = True

        else:
            print 'ERROR', folder
            sys.exit()


        print 'Making REFFILE for ', filetype
        print '%d files found with %d imsets'%(len(raw_files), n_imsets)

        if REFBIAS: 
            refbias_name = os.path.join( folder, 'refbias_%s_%s.fits'%(proposal, wk) )
            if os.path.exists( refbias_name ):
                print 'Refbias already created, skipping'
            else:
                if n_imsets > 120:
                    super_list = split_files( raw_files )
                    all_subnames = []
                    for i, sub_list in enumerate( super_list ):
                        subname = refbias_name.replace('.fits', '_grp0'+str(i+1)+'.fits')
                        print 'Making sub-file for datasets'
                        print sub_list
                        REFSTIS_refbias.make_refbias( sub_list, subname )
                        all_subnames.append( subname )
                    REFSTIS_functions.refaver( all_subnames, refbias_name )
                else:
                    REFSTIS_refbias.make_refbias( raw_files, refbias_name )

        if WEEKBIAS:
            weekbias_name = os.path.join( folder, 'weekbias_%s_%s.fits'%(proposal, wk) )
            if os.path.exists( weekbias_name ):
                print 'Weekbias already created, skipping'
            else:
                REFSTIS_refbias.make_refbias( raw_files, weekbias_name, basebias_name )
 
        if WEEKDARK:
            weekdark_name = os.path.join( folder, 'weekdark_%s_%s.fits'%(proposal, wk) )
            if os.path.exists( weekdark_name ):
                print 'Weekdark already created, skipping'
            else:
                weekbias_name = os.path.join( root_folder, 'biases/1-1x1', wk, 'refbias_%s_%s.fits'%(proposal, wk) ) ### probably need to be final file, either week* or ref*
                basedark_name = os.path.join( folder.replace(wk, 'all'), 'basedark.fits' )
                REFSTIS_weekdark.make_weekdark( raw_files, weekdark_name, weekbias_name , basedark_name)
    

#-------------------------------------------------------------------------------

def clean_directory( root_path ):
    """ Cleans directory of any fits files that do not end in _raw.fits 

    This WILL remove ANY files that do not match *_raw.fits.  This includes
    any plots, txt files, other fits files, anything.

    Use with caution

    """

    for root, dirs, files in os.walk( root_path ):
        for filename in files:
            if not filename.endswith('_raw.fits'):
                print 'Removing: ', filename
                os.remove( os.path.join( root, filename ) )

#-------------------------------------------------------------------------------

def get_new_obs(file_type, start, end):

    if file_type == 'DARK':
        proposal_list = dark_proposals
        MIN_EXPTIME = 1000
        MAX_EXPTIME = 1200
    elif file_type == 'BIAS':
        proposal_list = bias_proposals
        MIN_EXPTIME = -1
        MAX_EXPTIME = 100
    else:
        print 'file type not recognized: ', file_type

    query = support.SybaseInterface("ZEPPO", "dadsops")

    OR_part = "".join(["science.sci_pep_id = %d OR "%(proposal) for proposal in proposal_list])[:-3]

    #obs_name_query = "SELECT science.sci_data_set_name FROM science WHERE ( " + OR_part + " ) AND  science.sci_targname ='%s' AND science.sci_actual_duration BETWEEN %d AND %d "%(file_type, MIN_EXPTIME, MAX_EXPTIME)
    #start_time_query = "SELECT science.sci_start_time FROM science WHERE ( " + OR_part + " ) AND  science.sci_targname ='%s' AND science.sci_actual_duration BETWEEN %d AND %d "%(file_type, MIN_EXPTIME, MAX_EXPTIME)

    data_query = "SELECT science.sci_start_time,science.sci_data_set_name FROM science WHERE ( " + OR_part + " ) AND  science.sci_targname ='%s' AND science.sci_actual_duration BETWEEN %d AND %d "%(file_type, MIN_EXPTIME, MAX_EXPTIME)
    query.doQuery(query=data_query)
    new_dict = query.resultAsDict()
 
    #query.doQuery(query=obs_name_query)
    #new_dict = query.resultAsDict()
    #obs_names = new_dict[new_dict.keys()[0]][2:]  #remove non-obs entries in dictionary

    #query.doQuery(query=start_time_query)
    #new_dict = query.resultAsDict()
    #start_times = new_dict[new_dict.keys()[0]][2:]  #remove non-obs entries in dictionary
    
    obs_names = np.array( new_dict['sci_data_set_name'] )
    
    start_times_MJD = np.array( map(translate_date_string, new_dict['sci_start_time'] ) )
    
    index = np.where( (start_times_MJD > start) & (start_times_MJD < end) )[0]

    if not len( index ):
        print "WARNING: didn't find any datasets, skipping"
        return []

    assert start_times_MJD[index].min() > start, 'Data has mjd before period start'
    assert start_times_MJD[index].max() < end, 'Data has mjd after period end'

    datasets_to_retrieve = obs_names[index]

    return list(datasets_to_retrieve)

#-----------------------------------------------------------------------

def collect_new(observations_to_get):
    '''
    Function to find and retrieve new datasets for given proposal.
    Uses modules created by B. York: DADSAll.py and SybaseInterface.py.
    '''

    xml = createXmlFile(ftp_dir= retrieve_directory, 
                        set=observations_to_get, file_type='RAW')

    response = submitXmlFile(xml, 'dmsops1.stsci.edu')
    if ('SUCCESS' in response):
        success=True
    else:
        success=False
    success = True

    return success



#-----------------------------------------------------------------------

def separate_obs( base_dir, month_begin, month_end  ):
    all_files = glob.glob( os.path.join( base_dir, '*raw.fits') )
    #all_files = glob.glob( os.path.join( retrieve_directory, '*raw.fits') )

    print 'Separating', base_dir
    print
    print 'Data run from', month_begin, ' to ',  month_end

    mjd_times = np.array( [ pyfits.getval(item, 'EXPSTART', ext=1) for item in all_files ] )
    print 'Data goes from', mjd_times.min(),  ' to ',  mjd_times.max()

    print 'Making Lists'
    bias_111_list = [item for item in all_files if ( pyfits.getval(item,'TARGNAME', ext=0)=='BIAS') & ( pyfits.getval(item, 'CCDGAIN', ext=0) == 1) ]
    bias_411_list = [item for item in all_files if ( pyfits.getval(item, 'TARGNAME', ext=0)=='BIAS') & ( pyfits.getval(item, 'CCDGAIN', ext=0) == 4) ]
    dark_111_list = [item for item in all_files if ( pyfits.getval(item, 'TARGNAME', ext=0)=='DARK') & ( pyfits.getval(item, 'CCDGAIN', ext=0) == 1) ]
    print 'Done'

    for obs_list, file_type, mode in zip( [bias_111_list, dark_111_list, bias_411_list], 
                                        ['BIAS', 'DARK', 'BIAS'], 
                                        ['WK', 'WK', 'BIWK'] ):

        if len(obs_list) == 0:
            print '%s No obs to move.  Skipping'%(mode)
            continue
        gain = list( set( [ pyfits.getval(item, 'CCDGAIN', ext=0) for item in obs_list ] ) )
        print obs_list, 'gain = ', gain
        assert len(gain) == 1, 'ERROR: Not everything has the same gain'
        gain = gain[0]


        if mode == 'WK':
            N_periods = 4
        elif mode == 'BIWK':
            N_periods = 2

        anneal_weeks = REFSTIS_functions.divide_anneal_month(month_begin, month_end,'/grp/hst/stis/calibration/anneals/', N_periods)

        print
        print file_type, mode, 'will be broken up into %d periods as follows:'%(N_periods)
        print '\tWeek start, Week end'
        for a_week in anneal_weeks: print '\t', a_week
        print

        for period in range(N_periods):
            begin, end = anneal_weeks[period]
            week = str(period + 1) ##weeks from 1-4, not 0-3
            while len(week) < 2:
                week = '0'+week            

            output_path = base_dir
            if file_type == 'BIAS':
                output_path = os.path.join(output_path, 'biases/%d-1x1/%s%s/'%(gain, mode.lower(), week) )
            elif file_type == 'DARK':
                output_path = os.path.join(output_path, 'darks/%s%s/'%(mode.lower(), week) )
            else: print 'File Type not recognized'
            
            print output_path
            if not os.path.exists( output_path ): 
                os.makedirs( output_path )

            print 'week goes from: ', begin, end
            obs_to_move = [ item for item in obs_list if 
                            ( (pyfits.getval(item, 'EXPSTART', ext=1) >= begin) and 
                              (pyfits.getval(item, 'EXPSTART', ext=1) < end) ) ]
            print begin, end,  obs_to_move
            if not len(obs_to_move):
                print 'error, empty list to move'

            for item in obs_to_move:
                print 'Moving ', item,  ' to:', output_path
                shutil.move( item,  output_path )
                if 'IMPHTTAB' not in pyfits.getheader(os.path.join(output_path, item.split('/')[-1]), 0).keys():
                    pyfits.setval(os.path.join(output_path, item.split('/')[-1]), 'IMPHTTAB', ext = 0, value = 'oref$x9r1607mo_imp.fits')
                obs_list.remove( item )

#-----------------------------------------------------------------------

def move_obs(new_obs, base_output_dir):
    print 'Files not yet delivered.'
    delivered_set = set( [ os.path.split(item)[-1][:9].upper() for item in glob.glob( os.path.join(retrieve_directory, '*raw*.fits') ) ] )
    new_set = set( new_obs )

    while not new_set.issubset( delivered_set ):
        wait_minutes = 2
        time.sleep(wait_minutes * 60) #sleep for 2 min
        delivered_set = set( [ os.path.split(item)[-1][:9].upper() for item in glob.glob( os.path.join(retrieve_directory, '*raw*.fits') ) ] )

    assert len(new_obs) > 0, 'Empty list of new observations to move'

    if not os.path.exists( base_output_dir):
        os.makedirs( base_output_dir )

    list_to_move = [ os.path.join( retrieve_directory, item.lower()+'_raw.fits') for item in new_obs ]

    for item in list_to_move:
        print 'Moving ', item,  ' to:', base_output_dir
        shutil.move( item, base_output_dir )

    list_to_remove = glob.glob( os.path.join(retrieve_directory, '*.fits') )
    for item in list_to_remove:
        os.remove(item)

#-------------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
Description
------------------------------------
  Automated script for producing 
  dark and bias reference files for 
  STIS CCD data.

Locations
------------------------------------
  Product directory = EMPTY

Procedure
------------------------------------
1: pass

2: pass

------------------------------------
''' ) )

    parser.add_argument("-r",  "--redo_all", 
                      action='store_true', dest="redo_all",  default=False,
                      help="Re-run analysis on all past anneal months.")
    parser.add_argument("-c",  "--no_collect", 
                      action='store_false', dest="collect_new",  default=True, 
                      help="Turn off data collection function.")
    parser.add_argument("-p",  "--plots_only", 
                      action='store_true', dest="plots_only", default=False, 
                      help="Only remake plots and update the website.")
    parser.add_argument("-u", "--user_information", 
                        action='store', dest="user_information", default=None,
                        help="info string needed to request data")
    return parser.parse_args()

#-----------------------------------------------------------------------

def run():
    """ Run the reference file pipeline """

    args = parse_args()

    REFSTIS_pop_db.main()

    all_folders = get_new_periods()

    for folder in all_folders:
        make_ref_files( folder, clean=args.redo_all )

#-----------------------------------------------------------------------

if __name__ == "__main__":
    run()
