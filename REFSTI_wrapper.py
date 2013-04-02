#!/usr/bin/env python

"""Script to fill database with cos anneal month start and end times.  
"""

from __future__ import division

import sqlite3
import glob
import os
import sys
import pyfits
import support
import time
import shutil
import pdb
import numpy as np
from support import SybaseInterface
from support import createXmlFile, submitXmlFile
from REFSTI_functions import figure_days_in_period,figure_number_of_periods,translate_date_string
import REFSTI_pop_db


#### Needs 
#mjd converter
#get_directories
#anneal_dir

#products_directory = '/grp/hst/stis/referencefiles/darks_biases/'
#retrieve_directory = '/grp/hst/stis/referencefiles/darks_biases/requested/'
products_directory = '/user/ely/STIS/refstis_mark2/darks_biases/'
retrieve_directory = '/user/ely/STIS/refstis_mark2/requested/'

#dark_proposals = [7600, 7601, 8408, 8437, 8837, 8864, 8901, 8902, 9605, 9606, 
#                  10017, 10018, 11844, 11845, 12401, 12402, 12741, 12742]
#bias_proposals = [7600, 7601, 8409, 8439, 8838, 8865, 8903, 8904, 9607, 9608, 
#                  10019, 10020, 11846, 11847, 12403, 12404, 12743, 12744]

dark_proposals = [ 12741, 12742]
bias_proposals = [ 12743, 12744]

#-------------------------------------------------------------------------------

def list_retreieved_files( folder ):
    for root,dirs,files in os.walk( folder ):
        file_list = glob.glob('?????????_raw.fits')

#-------------------------------------------------------------------------------

def get_new_periods():
    print 'Reading from database'
    db = sqlite3.connect("/user/ely/STIS/refstis_mark2/my_scripts/anneal_info")
    c = db.cursor()
    table = 'anneals'

    c.execute("""SELECT * FROM %s """%(table))

    all_info = [row for row in c]

    table_id_all = [row[0] for row in all_info]
    proposal_id_all = [row[1] for row in all_info]
    visit_id_all = [ int(row[2]) for row in all_info ]
    anneal_start_all = [row[3] for row in all_info]
    anneal_end_all = [row[4] for row in all_info]

    for i in range( len(table_id_all) )[::-1]:
        if i == len(table_id_all)-1: continue
        ref_begin = anneal_end_all[i]
        ref_end = anneal_start_all[i+1]
        proposal = proposal_id_all[i]
        visit = visit_id_all[i+1]
        year,month,day,dec_year = support.mjd_to_greg(ref_begin)
        
        if visit < 10:
            visit = '0'+str(visit)
        else:
            visit = str(visit)

        print '%d_%d_%s'%(year,proposal,visit)
        products_folder = os.path.join( products_directory,'%d_%d_%s'%(year,proposal,visit) )
        if not os.path.exists( products_folder ): os.mkdir( products_folder )
        already_retrieved = [ os.path.split(item)[1][:9] for item in glob.glob( os.path.join(products_folder,'?????????_raw.fits') ) ]

        new_obs = get_new_obs('DARK',ref_begin,ref_end) + get_new_obs('BIAS',ref_begin,ref_end)  ###Grab both darks and biases before moving on to collect and move.
        obs_to_get = [ obs for obs in new_obs if not obs in already_retrieved ]
        #print new_obs

        collect_new( obs_to_get )
        move_obs( obs_to_get, output_folder) 

        #make_ref_files(products_folder, proposal, visit)

        '''
        if not os.path.exists(products_folder):
            print 'Making files for:'
            print products_folder
            print
            make_ref_files(products_folder, proposal, visit, ref_begin, ref_end)
        else:
            print 'Folders already exists for '
            print products_folder
            print 'I assume ref_files have already been created for this period.'
            print 'Skipping'
            print
         '''
#-------------------------------------------------------------------------------

def make_ref_files(output_folder, proposal, visit):
    pass
    #print new_obs

#-------------------------------------------------------------------------------

def get_new_obs(file_type, start, end):

    if file_type == 'DARK':
        proposal_list = dark_proposals
        MIN_EXPTIME = 1000
        MAX_EXPTIME = 1200
    elif file_type == 'BIAS':
        proposal_list = bias_proposals
        MIN_EXPTIME = 0
        MAX_EXPTIME = 100
    else:
        print 'file type not recognized: ',file_type

    query = support.SybaseInterface("ZEPPO","dadsops")

    OR_part = "".join(["science.sci_pep_id = %d OR "%(proposal) for proposal in proposal_list])[:-3]

    obs_name_query = "SELECT science.sci_data_set_name FROM science WHERE ( " + OR_part + " ) AND  science.sci_targname ='%s' AND science.sci_actual_duration BETWEEN %d AND %d "%(file_type,MIN_EXPTIME,MAX_EXPTIME)
    start_time_query = "SELECT science.sci_start_time FROM science WHERE ( " + OR_part + " ) AND  science.sci_targname ='%s' AND science.sci_actual_duration BETWEEN %d AND %d "%(file_type,MIN_EXPTIME,MAX_EXPTIME)

    print obs_name_query
    print start_time_query

    query.doQuery(query=obs_name_query)
    new_dict = query.resultAsDict()
    obs_names = new_dict[new_dict.keys()[0]][2:]  #remove non-obs entries in dictionary

    query.doQuery(query=start_time_query)
    new_dict = query.resultAsDict()
    start_times = new_dict[new_dict.keys()[0]][2:]  #remove non-obs entries in dictionary

    obs_names = np.array( obs_names )
    start_times_MJD = np.array( [ translate_date_string(item) for item in start_times ] )
    
    index = np.where( (start_times_MJD > start) & (start_times_MJD < end) )[0]

    datasets_to_retrieve = obs_names[index]
    dataset_times = start_times_MJD[index]

    return list(datasets_to_retrieve)

#-----------------------------------------------------------------------

def collect_new(observations_to_get):
    '''
    Function to find and retrieve new datasets for given proposal.
    Uses modules created by B. York: DADSAll.py and SybaseInterface.py.
    '''
    print '#----------------------------#'
    print 'Searching for new observations'
    print '#----------------------------#'

    xml = createXmlFile(ftp_dir= retrieve_directory, 
                      set=observations_to_get, file_type='RAW', 
                      archive_user='ely', 
                      archive_pwd='etacar', email='ely@stsci.edu', 
                      host='science3.stsci.edu', 
                      ftp_user='ely', ftp_pwd='5correctmonkeys')

    #response = submitXmlFile(xml,'dmsops1.stsci.edu')
    #print response
    #if ('SUCCESS' in response):
    #    success=True
    #else:
    #    success=False
    success = True

    return success

#-----------------------------------------------------------------------

def move_obs(new_obs, base_output_dir):
    print 'Files not yet delivered.'
    while len(glob.glob(os.path.join(retrieve_directory,'*raw*.fits')))!=len(new_obs):
        wait_minutes=2
        time.sleep(wait_minutes*60) #sleep for 2 min

    list_to_move = glob.glob(os.path.join(retrieve_directory,'*raw*.fits'))
    assert len(list_to_move)>0, 'Empty list of new observations to move'

    print 'Making Lists'
    bias_111_list = [item for item in list_to_move if ( pyfits.getval(item,'TARGNAME',ext=0)=='BIAS') & ( pyfits.getval(item,'CCDGAIN',ext=0) == 1) ]
    bias_411_list = [item for item in list_to_move if ( pyfits.getval(item,'TARGNAME',ext=0)=='BIAS') & ( pyfits.getval(item,'CCDGAIN',ext=0) == 4) ]
    dark_111_list = [item for item in list_to_move if ( pyfits.getval(item,'TARGNAME',ext=0)=='DARK') & ( pyfits.getval(item,'CCDGAIN',ext=0) == 1) ]
    print 'Done'

    if not os.path.exists(base_output_dir):
        print 'Making output folders'
        os.mkdir(base_output_dir)
        for ending in ['biases','biases/1-1x1',
                       'biases/4-1x1',
                       'delivery','darks']:
            os.mkdir( os.path.join(base_output_dir,ending) )

    for obs_list,file_type,mode in zip( [bias_111_list,dark_111_list,bias_411_list],
                                        ['BIAS','DARK','BIAS'],
                                        ['WK','Wk','BIWK'] ):
        gain = list( set( [ pyfits.getval(item,'CCDGAIN',ext=0) for item in obs_list ] ) )
        assert len(gain) == 1, 'ERROR: Not everything has the same gain'
        gain = gain[0]

        obs_day_list = list( set( [ pyfits.getval(item,'TDATEOBS',ext=0) for item in obs_list ] ) )  #forgive the confusion
        obs_day_list.sort()

        N_days_total = len(obs_day_list)
        N_periods = figure_number_of_periods(N_days_total, mode)

        N_days_per_period = figure_days_in_period(N_periods, N_days_total)
        print file_type,mode,N_days_total,N_periods
        begin = 0
        end = 0
        for period in range(N_periods):
            period += 1 ##weeks from 1-4, not 0-3
            if period < 10:
                week = '0%d'%(period)
            else:
                week = '%d'%(period)

            output_path = base_output_dir
            if file_type == 'BIAS':
                output_path += 'biases/%d-1x1/%s%s/'%(gain,mode.lower(),week)
            elif file_type == 'DARK':
                output_path += 'darks/%s%s/'%(mode.lower(),week)
            else: print 'File Type not recognized'
            if not os.path.exists( output_path ): os.mkdir( output_path )

            end += N_days_per_period[period]
            which_days = obs_day_list[begin:end]
            begin += end

            obs_to_move = [ item for item in obs_list if pyfits.getval(item,'TDATEOBS',ext=0) in which_days ]

            for item in obs_to_move:
                print 'Moving ',item, ' to:',output_path
                shutil.move( item, output_path )

    #list_to_remove=glob.glob(os.path.join(retrieve_directory,'*.fits'))
    #for item in list_to_remove:
    #    os.remove(item)

#-------------------------------------------------------------------------------

def main():
    REFSTI_pop_db.main()

    get_new_periods()

if __name__ == "__main__":
    main()
