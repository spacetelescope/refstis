"""Create STIS Superdarks and Superbiases for the CCD detector.

"""

import sqlite3
import glob
import os
import sys
from astropy.io import fits
import argparse
import textwrap
import time
import shutil
import re
import numpy as np
import yaml
import getpass

from .delivery import check_all
from .functions import figure_number_of_periods, translate_date_string, mjd_to_greg
from .retrieval import submit_xml_request, build_xml_request, everything_retrieved
from . import pop_db
from . import basedark
from . import weekdark
from . import refbias
from . import weekbias
from . import basejoint
from . import functions

#-------------------------------------------------------------------------------

def get_new_periods(products_directory, settings):
    print('#-------------------#')
    print('Reading from database')
    print('#-------------------#\n')
    db = sqlite3.connect( "anneal_info.db" )
    c = db.cursor()
    table = 'anneals'

    c.execute("""SELECT * FROM %s """ % (table))

    all_info = [row for row in c]

    table_id_all = [row[0] for row in all_info]
    proposal_id_all = [row[1] for row in all_info]
    visit_id_all = [int(row[2]) for row in all_info]
    anneal_start_all = [row[3] for row in all_info]
    anneal_end_all = [row[4] for row in all_info]

    dirs_to_process = []

    for i in range(len(table_id_all))[::-1]:
        if i == len(table_id_all) - 1: continue
        ref_begin = anneal_end_all[i]
        ref_end = anneal_start_all[i + 1]
        #-- defined from proposal of the next anneal
        proposal = proposal_id_all[i + 1]

        visit = visit_id_all[i + 1]
        year, month, day, dec_year = mjd_to_greg(ref_begin)
        end_year, end_month, end_day, dec_year = mjd_to_greg(ref_end)

        if visit < 10:
            visit = '0' + str(visit)
        else:
            visit = str(visit)

        print('#--------------------------------#')
        print('Searching for new observations for')
        print('Period: %d_%d_%s'%(year, proposal, visit))
        print('MJD %5.5f %5.5f'%(ref_begin, ref_end))
        print(month, day, year, ' to ', end_month, end_day, end_year)
        print('#--------------------------------#')

        products_folder = os.path.join(products_directory,
                                       '%d_%s' % (proposal, visit))

        #-- maybe put this somewhere else?
        #-- to only process datasets if new data has been retrieved?
        dirs_to_process.append(products_folder)

        if not os.path.exists(products_folder):
            os.makedirs(products_folder)

        already_retrieved = []
        for root, dirs, files in os.walk(products_folder):
            for filename in files:
                if filename.endswith('_raw.fits'):
                    already_retrieved.append(filename[:9])

        new_obs = get_new_obs('DARK', ref_begin, ref_end, settings) + \
                  get_new_obs('BIAS', ref_begin, ref_end, settings)

        obs_to_get = [obs for obs in new_obs if not obs in already_retrieved]

        if not len(obs_to_get):
            print('No new obs to get, skipping this period\n\n')
            continue
        else:
            print('Found new observations for this period')
            print(obs_to_get, '\n\n')

        response = collect_new(obs_to_get, settings)
        move_obs(obs_to_get, products_folder, settings['retrieve_directory'])
        separate_obs(products_folder, ref_begin, ref_end)

    return dirs_to_process

#-------------------------------------------------------------------------------

def split_files(all_files):
    """Split file list into two smaller lists for iraf tasks

    Each list will have a selection from both early and late in time.

    Parameters
    ----------
    all_files : list
        full list of files

    Returns
    -------
    super_list : list
        list of lists containing the split list of all files

    """

    all_info = [(fits.getval(filename,'EXPSTART', 1), filename)
                for filename in all_files]
    all_info.sort()
    all_files = [line[1] for line in all_info]

    super_list = [all_files[0::2],
                  all_files[1::2]]

    return super_list

#-------------------------------------------------------------------------------

def pull_out_subfolders(root_folder):
    """ Walk through input folder and use regular expressions to return lists of
    the folders for each gain and for each week

    Parameters
    ----------
    root_folder
        string, the folder to walk through

    Returns
    -------
    gain_folders
        list, containing folders for each gain
    week_folders
        list, containing folders for each week

    """

    gain_folders = []
    week_folders = []
    for root, dirs, files in os.walk(root_folder):
        tail = os.path.split(root)[-1]

        if 'wk' in tail:
            week_folders.append(root)

        # ex: 1-1x1 or 4-1x1
        if re.search('([0-4]-[0-4]x[0-4])', tail):
            gain_folders.append(root)

    return gain_folders, week_folders

#-------------------------------------------------------------------------------

def grab_between(file_list, mjd_start, mjd_end):
    """Select files between given start/end times

    Parameters
    ----------
    file_list : list
        files on which to filter
    mjd_start : float
        earliest time
    mjd_end : float
        latest time

    Yields
    ------
        filenames with mjd_start < TEXPSTRT < mjd_end
    """

    for filename in file_list:
        with pyfits.open(filename) as hdu:
            data_start = hdu[0].header['TEXPSTRT']
            data_end = hdu[0].header['TEXPEND']

        if mjd_start < data_start < mjd_end:
            yield filename

#-------------------------------------------------------------------------------

def pull_info(foldername):
    """ Pull proposal and week number from folder name

    A valid proposal is a string of 5 numbers from 0-9
    A valid week is a string of 'wk' + 2 numbers ranging from 0-9.
    A valid biweek is the string 'biwk' + 2 numbers from 0-9.

    Parameters
    ----------
    foldername
        string, name of folder to search in

    Returns
    -------
    proposal
        string, proposal number as string
    week
        string, week of the anneal (wk01, etc)

    """

    try:
        proposal, visit = re.findall('([0-9]{5})_([0-9]{2})', foldername)[0]
    except:
        proposal, visit = '', ''

    try:
        week = re.findall('([bi]*wk0[0-9])', foldername)[0]
    except:
        week = ''

    return proposal, week, visit

#-------------------------------------------------------------------------------

def get_anneal_month(proposal_id, anneal_id):

    db = sqlite3.connect( "anneal_info.db" )
    c = db.cursor()
    table = 'anneals'


    c.execute("""SELECT id,start FROM {} WHERE proposid={} AND visit={}""".format(table,
                                                                                  proposal_id,
                                                                                  anneal_id))

    all_info = [row for row in c]
    if len(all_info) > 1:
        raise ValueError("Too many values returned: {}".format(all_info))
    pri_key, anneal_end = all_info[0]



    c.execute("""SELECT end FROM {} WHERE id={}""".format(table, pri_key - 1))

    all_info = [row for row in c]
    if len(all_info) > 1:
        raise ValueError("Too many values returned: {}".format(all_info))
    anneal_start = all_info[0][0]

    return anneal_start, anneal_end

#-------------------------------------------------------------------------------

def make_pipeline_reffiles(root_folder, last_basedark=None, last_basebias=None):
    """Make reference files like the refstis pipeline

    1.  Separate dark and bias datasets into week folders
    2.

    """

    if not 'oref' in os.environ:
        raise ValueError("oref hasn't been defined in the environment")

    bias_threshold = {(1, 1, 1) : 98,
                      (1, 1, 2) : 25,
                      (1, 2, 1) : 25,
                      (1, 2, 2) : 7,
                      (1, 4, 1) : 7,
                      (1, 4, 2) : 4,
                      (4, 1, 1) : 1}

    print('#-----------------------------#')
    print('#  Making all ref files for   #')
    print(root_folder)
    print('#-----------------------------#')

    if not os.path.exists(root_folder):
        raise IOError('Root folder does not exist')

    # separate raw files in folder into periods
    separate_period(root_folder)

    print('###################')
    print(' make the basebias ')
    print('###################')
    raw_files = []

    for root, dirs, files in os.walk(os.path.join(root_folder, 'biases')):
        if not '1-1x1' in root:
            continue

        for filename in files:
            if filename.startswith('o') and filename.endswith('_raw.fits'):
                raw_files.append(os.path.join(root, filename))

    basebias_name = os.path.join(root_folder, 'basebias.fits')
    if os.path.exists(basebias_name):
        print('{} already exists, skipping')
    else:
        basejoint.make_basebias(raw_files, basebias_name)

    print('#######################')
    print(' make the weekly biases')
    print('#######################')
    #-- Find the premade folders if they exist
    gain_folders, week_folders = pull_out_subfolders(root_folder)

    #-- use last basefiles if supplied
    basebias_name = last_basebias or basebias_name

    bias_folders = [item for item in week_folders if '/biases/' in item]
    print(bias_folders)

    for folder in sorted(bias_folders):
        print('Processing {}'.format(folder))

        proposal, wk, visit = pull_info(folder)

        raw_files = glob.glob(os.path.join(folder, '*raw.fits'))
        n_imsets = functions.count_imsets(raw_files)

        gain = functions.get_keyword(raw_files, 'CCDGAIN', 0)
        xbin = functions.get_keyword(raw_files, 'BINAXIS1', 0)
        ybin = functions.get_keyword(raw_files, 'BINAXIS2', 0)

        weekbias_name = os.path.join(folder,
                                     'weekbias_%s_%s_%s_bia.fits'%(proposal, visit, wk))
        if os.path.exists(weekbias_name):
            print('{} already exists, skipping')
            continue

        #make weekbias if too few imsets

        if n_imsets < bias_threshold[(gain, xbin, ybin)]:
            weekbias.make_weekbias(raw_files, weekbias_name, basebias_name)
        else:
            if n_imsets > 120:
                super_list = split_files(raw_files)
                all_subnames = []
                for i, sub_list in enumerate(super_list):
                    subname = weekbias_name.replace('.fits', '_grp0'+str(i+1)+'.fits')
                    print('Making sub-file for datasets')
                    print(sub_list)
                    refbias.make_refbias(sub_list, subname)
                    all_subnames.append(subname)
                functions.refaver(all_subnames, weekbias_name)

            else:
                refbias.make_refbias(raw_files, weekbias_name)


    print('#######################')
    print(' make the weekly darks ')
    print('#######################')
    all_flt_darks = []
    dark_folders = [item for item in week_folders if '/darks/' in item]
    print(dark_folders)
    for folder in sorted(dark_folders):
        print('Prepping {}'.format(folder))

        proposal, wk, visit = pull_info(folder)

        weekbias_name = os.path.join(root_folder,
                                     'biases/1-1x1',
                                     wk,
                                     'weekbias_%s_%s_%s_bia.fits'%(proposal, visit, wk))

        raw_files = glob.glob(os.path.join(folder, '*raw.fits'))
        for item in raw_files:
            print("bias subtracting")
            flt_name = functions.bias_subtract_data(item, weekbias_name)
            all_flt_darks.append(flt_name)


    for folder in sorted(dark_folders):
        print('Processing {}'.format(folder))
        raw_files = glob.glob(os.path.join(folder, '*flt.fits'))
        n_imsets = functions.count_imsets(raw_files)

        gain = functions.get_keyword(raw_files, 'CCDGAIN', 0)
        xbin = functions.get_keyword(raw_files, 'BINAXIS1', 0)
        ybin = functions.get_keyword(raw_files, 'BINAXIS2', 0)

        proposal, wk, visit = pull_info(folder)
        weekdark_name = os.path.join(folder,
                                     'weekdark_%s_%s_%s_drk.fits'%(proposal, visit, wk))
        if os.path.exists(weekdark_name):
            print('{} already exists, skipping')
            continue


        #-- create a base-dark for everyweek
        weekbias_name = os.path.join(root_folder,
                                     'biases/1-1x1',
                                     wk,
                                     'weekbias_%s_%s_%s_bia.fits'%(proposal, visit, wk))
        basedark_name = os.path.join(folder, 'basedark_%s_%s_%s.fits'%(proposal, visit, wk))

        if not os.path.exists(basedark_name):
            basedark.make_basedark(all_flt_darks, basedark_name, weekbias_name)

        basedark_name = last_basedark or basedark_name
        weekdark.make_weekdark(raw_files,
                               weekdark_name,
                               basedark_name,
                               weekbias_name)

#-------------------------------------------------------------------------------

def clean_directory(root_path):
    """ Cleans directory of any fits files that do not end in _raw.fits

    .. warning::

      This WILL remove ANY files that do not match \*_raw.fits.  This includes
      any plots, txt files, other fits files, anything.

    """

    for root, dirs, files in os.walk(root_path):
        for filename in files:
            if not filename.endswith('_raw.fits'):
                print('Removing: ', filename)
                os.remove(os.path.join( root, filename))

#-------------------------------------------------------------------------------

def reset(folder):
    for root, dirs, files in os.walk(folder):
        for filename in files:
            if '_raw.fits' in filename:
                full = os.path.join(root, filename)
                if not os.path.exists(os.path.join(folder, filename)):
                    shutil.move(full, folder)

    os.system('rm -fR {}'.format(os.path.join(folder, 'darks')))
    os.system('rm -fR {}'.format(os.path.join(folder, 'biases')))
    for item in glob.glob('{}/basedark.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/basedark_?????_*.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/basebias.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/*_flt.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/*_crj.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/*_blev.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/weekdark*.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/weekbias*.fits'.format(folder)):
        os.remove(item)

    for item in glob.glob('{}/*.txt'.format(folder)):
        os.remove(item)

#-------------------------------------------------------------------------------

def get_new_obs(file_type, start, end, settings):

    if file_type == 'DARK':
        proposal_list = settings['dark_proposals']
        MIN_EXPTIME = 1000
        MAX_EXPTIME = 1200
    elif file_type == 'BIAS':
        proposal_list = settings['bias_proposals']
        MIN_EXPTIME = -1
        MAX_EXPTIME = 100
    else:
        print('file type not recognized: ', file_type)

    # Gather configuration settings
    mast_server = settings['mast_server']
    mast_database = settings['mast_database']
    mast_account = settings['mast_account']
    mast_password = settings['mast_password']

    # Connect to server
    connection = "tsql -S {0} -D '{1}' -U '{2}' -P '{3}' -t '|'".format(mast_server, mast_database, mast_account, mast_password)
    transmit, receive = os.popen2(connection)

    OR_part = "".join(["science.sci_pep_id = %d OR "%(proposal) for proposal in proposal_list])[:-3]
    data_query = """SELECT science.sci_start_time,science.sci_data_set_name
                           FROM science INNER JOIN stis_ref_data
                                ON science.sci_data_set_name = stis_ref_data.ssr_data_set_name
                            WHERE ( {} ) AND
                                  stis_ref_data.ssr_ccdamp = 'D' AND
                                  science.sci_targname = '{}' AND
                                  science.sci_actual_duration BETWEEN {} AND {}
                                  \ngo\n""".format(OR_part, file_type, MIN_EXPTIME, MAX_EXPTIME)

    # Perform query and capture results
    transmit.write(data_query)
    transmit.close()
    mast_results = receive.readlines()
    receive.close()

    # Prune mast_results of unwanted information
    mast_results = mast_results[7:-2]
    mast_start_times = [item.strip().split('|')[0] for item in mast_results]
    mast_rootnames = [item.lower().strip().split('|')[1] for item in mast_results]

    obs_names = np.array(mast_rootnames)
    start_times_MJD = np.array(list(map(translate_date_string, mast_start_times)))
    index = np.where((start_times_MJD > start) & (start_times_MJD < end))[0]

    if not len(index):
        print("WARNING: didn't find any datasets, skipping")
        return []

    assert start_times_MJD[index].min() > start, 'Data has mjd before period start'
    assert start_times_MJD[index].max() < end, 'Data has mjd after period end'

    datasets_to_retrieve = obs_names[index]

    return list(datasets_to_retrieve)

#-----------------------------------------------------------------------

def collect_new(observations_to_get, settings):
    '''
    Function to find and retrieve new datasets for given proposal.
    '''

    xml = build_xml_request(observations_to_get, settings)
    response = submit_xml_request(xml, settings)

    username = getpass.getuser()
    tracking_id = re.search("("+username+"[0-9]{5})", response).group()

    if not 'SUCCESS' in response:
        return False

    done = False
    killed = False

    while not done:
        print("waiting for files to be delivered")
        time.sleep(60)
        done, killed = everything_retrieved(tracking_id)

        if killed:
            return False

    return True

#-----------------------------------------------------------------------

def separate_period(base_dir):
    """Separate observations in the base dir into needed folders.

    Parameters
    ----------
    base_dir, str
        directory containing darks and biases to be split.

    """


    print('Separating', base_dir)
    all_files = glob.glob(os.path.join(base_dir, 'o*_raw.fits'))
    if not len(all_files):
        print("nothing to move")
        return

    mjd_times = np.array([fits.getval(item, 'EXPSTART', ext=1)
                          for item in all_files])
    month_begin = mjd_times.min()
    month_end = mjd_times.max()
    print('All data goes from', month_begin, ' to ',  month_end)

    select_gain = {'WK' : 1,
                   'BIWK' : 4}

    for file_type, mode in zip(['BIAS', 'DARK', 'BIAS'],
                               ['WK', 'WK', 'BIWK']):

        gain = select_gain[mode]

        obs_list = []
        for item in all_files:
            with fits.open(item) as hdu:
                if (hdu[0].header['TARGNAME'] == file_type) and (hdu[0].header['CCDGAIN'] == gain):
                    obs_list.append(item)

        if not len(obs_list):
            print('{} No obs to move.  Skipping'.format(mode))
            continue
        else:
            print(file_type,  mode, len(obs_list), 'files to move, ', 'gain = ', gain)

        N_days = int(month_end - month_begin)
        N_periods = figure_number_of_periods(N_days, mode)
        week_lengths = functions.figure_days_in_period(N_periods, N_days)

        #--Add remainder to end
        week_lengths[-1] += (month_end - month_begin) - N_days

        #-- Translate to MJD
        anneal_weeks = []
        start = month_begin
        end = start + week_lengths[0]
        anneal_weeks.append((start, end))
        for item in week_lengths[1:]:
            start = end
            end += item
            anneal_weeks.append((start, end))

        print()
        print(file_type, mode, 'will be broken up into %d periods as follows:'%(N_periods))
        print('\tWeek start, Week end')
        for a_week in anneal_weeks:
            print('\t', a_week)
        print()

        for period in range(N_periods):
            begin, end = anneal_weeks[period]
            # weeks from 1-4, not 0-3
            week = str(period + 1)
            while len(week) < 2:
                week = '0' + week

            output_path = base_dir
            if file_type == 'BIAS':
                output_path = os.path.join(output_path,
                                           'biases/%d-1x1/%s%s/'%(gain,
                                                                  mode.lower(),
                                                                  week))
            elif file_type == 'DARK':
                output_path = os.path.join(output_path,
                                           'darks/%s%s/'%(mode.lower(), week))
            else:
                print('File Type not recognized')

            print(output_path)
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            print('week goes from: ', begin, end)
            obs_to_move = [item for item in obs_list if
                           (begin <= fits.getval(item, 'EXPSTART', ext=1) <= end)]

            if not len(obs_to_move):
                raise ValueError('error, empty list to move')

            for item in obs_to_move:
                print('Moving ', item,  ' to:', output_path)
                shutil.move(item,  output_path)
                if not 'IMPHTTAB' in fits.getheader(os.path.join(output_path,
                                                                 item.split('/')[-1]), 0):
                    ###Dynamic at some point
                    fits.setval(os.path.join(output_path, item.split('/')[-1]),
                                'IMPHTTAB',
                                ext=0,
                                value='oref$x9r1607mo_imp.fits')

                obs_list.remove(item)
                all_files.remove(item)

#-----------------------------------------------------------------------


def separate_obs(base_dir, month_begin, month_end):
    print('Separating', base_dir)
    print()
    print('Period runs from', month_begin, ' to ',  month_end)

    all_files = glob.glob(os.path.join(base_dir, '*raw.fits'))

    mjd_times = np.array([fits.getval(item, 'EXPSTART', ext=1)
                          for item in all_files])
    print('All data goes from', mjd_times.min(), ' to ',  mjd_times.max())

    select_gain = {'WK' : 1,
                   'BIWK' : 4}

    for file_type, mode in zip(['BIAS', 'DARK', 'BIAS'],
                               ['WK', 'WK', 'BIWK']):

        gain = select_gain[mode]

        obs_list = []
        for item in all_files:
            with fits.open(item) as hdu:
                if (hdu[0].header['TARGNAME'] == file_type) and (hdu[0].header['CCDGAIN'] == gain):
                    obs_list.append(item)

        if not len(obs_list):
            print('%s No obs to move.  Skipping'%(mode))
            continue

        print(file_type,  mode, len(obs_list), 'files to move, ', 'gain = ', gain)

        N_days = int(round(month_end - month_begin))
        N_periods = figure_number_of_periods(N_days, mode)
        week_lengths = functions.figure_days_in_period(N_periods, N_days)

        anneal_weeks = [(month_begin + item - week_lengths[0], month_begin + item) for item in np.cumsum(week_lengths)]

        print()
        print(file_type, mode, 'will be broken up into %d periods as follows:'%(N_periods))
        print('\tWeek start, Week end')
        for a_week in anneal_weeks:
            print('\t', a_week)
        print()

        for period in range(N_periods):
            begin, end = anneal_weeks[period]
            # weeks from 1-4, not 0-3
            week = str(period + 1)
            while len(week) < 2:
                week = '0'+week

            output_path = base_dir
            if file_type == 'BIAS':
                output_path = os.path.join(output_path,
                                           'biases/%d-1x1/%s%s/'%(gain,
                                                                  mode.lower(),
                                                                  week))
            elif file_type == 'DARK':
                output_path = os.path.join(output_path,
                                           'darks/%s%s/'%(mode.lower(), week))
            else:
                print('File Type not recognized')

            print(output_path)
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            print('week goes from: ', begin, end)
            obs_to_move = [item for item in obs_list if
                            ((fits.getval(item, 'EXPSTART', ext=1) >= begin) and
                             (fits.getval(item, 'EXPSTART', ext=1) < end))]

            if not len(obs_to_move):
                print('error, empty list to move')

            for item in obs_to_move:
                print('Moving ', item,  ' to:', output_path)
                shutil.move(item,  output_path)
                if not 'IMPHTTAB' in fits.getheader(os.path.join(output_path,
                                                                 item.split('/')[-1]), 0):
                    ###Dynamic at some point
                    fits.setval(os.path.join(output_path, item.split('/')[-1]),
                                'IMPHTTAB',
                                ext=0,
                                value='oref$x9r1607mo_imp.fits')
                obs_list.remove(item)
                all_files.remove(item)

#-----------------------------------------------------------------------

def move_obs(new_obs, base_output_dir, retrieve_directory):
    assert len(new_obs) > 0, 'Empty list of new observations to move.'

    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    list_to_move = [os.path.join(retrieve_directory, item.lower()+'_raw.fits') for item in new_obs]

    for item in list_to_move:
        print('Moving ', item,  ' to:', base_output_dir)
        shutil.move( item, base_output_dir )

    list_to_remove = glob.glob(os.path.join(retrieve_directory, '*.fits'))
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

    parser.add_argument("-r",
                        "--redo_all",
                        action='store_true',
                        dest="redo_all",
                        default=False,
                        help="Re-run analysis on all past anneal months.")
    parser.add_argument("-c",
                        "--no_collect",
                        action='store_false',
                        dest="collect_new",
                        default=True,
                        help="Turn off data collection function.")
    parser.add_argument("-p",
                        "--plots_only",
                        action='store_true',
                        dest="plots_only",
                        default=False,
                        help="Only remake plots and update the website.")
    parser.add_argument("-u",
                        "--user_information",
                        action='store',
                        dest="user_information",
                        default=None,
                        help="info string needed to request data")
    parser.add_argument("-m",
                        "--reprocess_month",
                        action='store',
                        nargs=2,
                        type=str,
                        help="Which dates during which you'd like to process",
                        default=None,
                        dest = "reprocess_month")

    return parser.parse_args()

#-----------------------------------------------------------------------

def run(config_file='config.yaml'):
    """Run the reference file pipeline """

    args = parse_args()

    print(args)

    #-- start pipeline configuration
    print(os.path.abspath(config_file))
    if not os.path.exists(config_file):
        raise IOError("Can't open configure file: {}".format(config_file))

    with open(config_file, 'r') as f:
        data = yaml.load(f)
        products_directory = data['products_directory'] #AER 23 August 2016

    for location in [data['products_directory'], data['retrieve_directory'], data['delivery_directory']]:
        if not os.path.isdir(location):
            os.makedirs(location)

    if not 'oref' in os.environ:
        if not 'oref' in data:
            raise KeyError("oref environment must be set to run the pipeline.")
        else:
            os.environ['oref'] = data['oref']

    pop_db.main()

    all_folders = get_new_periods(data['products_directory'], data)


    if args.redo_all: # AER 11 Aug 2016: in an attempt to use the commandline args.
        print("----------------------------------")
        print("Processing all past anneal months")
        print("----------------------------------")
        print('all_folders: {}'.format(all_folders))
        for folder in all_folders:
            make_pipeline_reffiles(folder)
            tail = folder.rstrip(os.sep).split(os.sep)[-1]
            destination = os.path.join(data['delivery_directory'], tail)
            check_all(folder, destination)

    # AER 11 Aug 2016
    if not args.redo_all and not args.reprocess_month:
        print("-----------------------------------")
        print("Processing most recent anneal month")
        print("-----------------------------------")
        make_pipeline_reffiles(all_folders[0])
        tail = all_folders[0].rstrip(os.sep).split(os.sep)[-1]
        destination = os.path.join(data['delivery_directory'], tail)
        check_all(all_folders[0], destination)

    # AER 12 Aug 2016
    if args.reprocess_month:
        print("------------------------------------------------")
        print("Processing files between {} and {}".format(args.reprocess_month[0], args.reprocess_month[1]))
        print("------------------------------------------------")


        filestoprocess = []
        print('products_directory:  {}'.format(products_directory))
        for all_anneals in glob.glob(''.join([products_directory, '?????_??/darks/'])):
            for root, directories, files_all in os.walk(all_anneals):
                print(directories)
                if not directories:
                    fltfiles = glob.glob(''.join([root, '/*_flt.fits']))
                    if len(fltfiles) != 0:
                        onefile = np.sort(fltfiles)[0]
                        obsdate = fits.getval(onefile, 'TDATEOBS', ext = 0)
                        if (obsdate <= args.reprocess_month[1] and obsdate >= args.reprocess_month[0]):
                            filestoprocess.append(root)
        print('filestoprocess:  {}'.format(filestoprocess))
        folders1 = []
        for f in filestoprocess:
            print('/'.join(f.split('/')[:7]))
            folders1.append('/'.join(f.split('/')[:7]))
        all_folders = set(folders1)

        for folder in all_folders:
            make_pipeline_reffiles(folder)
            tail = folder.rstrip(os.sep).split(os.sep)[-1]
            destination = os.path.join(data['delivery_directory'], tail)
            check_all(folder, destination)
#-----------------------------------------------------------------------
