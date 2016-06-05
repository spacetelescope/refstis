"""
Script to populate database with STIS anneal month start and end times.

This is used by the STIS darks and bias reference file pipeline to determine
the start and end of each anneal period, and to retrieve the appropriate dark
and bias files.

"""

import sqlite3
from astropy.io import fits as pyfits
import os
import glob

anneal_dir = '/grp/hst/stis/calibration/monitors/ccd/anneals/'

#-------------------------------------------------------------------------------

def get_directories():
    """ Loop over the anneal monitor directory and return list of directories
    containing anneal observations

    returns
    -------
    directories
        list, containing directory paths

    """

    directories = []
    for year in range(2011, 2020):   #change back to 1996 to create older darks
        for month in ('01', '02', '03', '04', '05', '06',
                      '07', '08', '09', '10', '11', '12'):
            for last in ('/', 'a/', 'b/'):
                path = anneal_dir + str(year) + '_' + month + last
                if os.path.exists(path) and path != anneal_dir+'2010_01/':
                    crj_list = glob.glob(path + '?????????_crj.fits')
                    crj_list.sort()
                    if len(crj_list) == 2:
                        print(path)
                        directories.append( path )
    return directories

#-------------------------------------------------------------------------------

def grab_anneal_mjds():
    """  Loop over anneal directories and return a list of tuples containing
    useful mjd info to be populated into anneal database

    returns
    -------
    anneal_info
        list, full of tupes containing information for database

    """

    print('Getting anneal information')
    anneal_info = []
    for directory in get_directories():
        anneal_obs = glob.glob( os.path.join(directory, '?????????_crj.fits') )
        anneal_obs.sort()
        if len(anneal_obs) != 2:
            print('Error in ', directory)
            continue

        proposid = pyfits.getval( anneal_obs[0], 'PROPOSID', ext=0 )
        visit_number = int(pyfits.getval(anneal_obs[1], 'OBSET_ID')) - 1
        anneal_start = pyfits.getval(anneal_obs[0], 'TEXPSTRT', ext=0)
        anneal_end = pyfits.getval(anneal_obs[1], 'TEXPSTRT', ext=0)

        anneal_info.append( (proposid, visit_number, anneal_start, anneal_end) )

    return anneal_info

#-------------------------------------------------------------------------------

def pop_database(anneal_info):
    """  Populate anneal database with information gained pulled from anneal
    monitor observations

    parameters
    ----------
    anneal_info
        list, list of tuples to be populated into db

    returns
    -------
    None

    """

    print('Populating Anneal database')
    db = sqlite3.connect("anneal_info.db")
    c = db.cursor()
    table = 'anneals'
    try:
        c.execute("""CREATE TABLE %s (id integer PRIMARY KEY, proposid integer, visit real, start real, end real)""" % (table))
    except sqlite3.OperationalError:
        c.execute("""DROP TABLE %s""" % (table))
        c.execute("""CREATE TABLE %s (id integer PRIMARY KEY, proposid integer, visit real, start real, end real)""" % (table))

    for i, line in enumerate(anneal_info):
        proposid = line[0]
        visit_number = line[1]
        start = line[2]
        end = line[3]
        c.execute("""INSERT INTO %s VALUES (?,?,?,?,?)""" %
                  (table), (i+1, proposid, visit_number, start, end))
    db.commit()
    c.execute("""SELECT * FROM %s """ % (table))
    print('Database Populated')
    print('------------------')
    for row in c:
        print(row)

#-------------------------------------------------------------------------------

def main():
    """ Main function to retrieve anneal info and populate database """
    anneal_stats = grab_anneal_mjds()
    pop_database(anneal_stats)

#-------------------------------------------------------------------------------
