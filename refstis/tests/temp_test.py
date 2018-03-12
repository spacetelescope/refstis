#!/usr/bin/python

from refstis import functions, basedark, weekdark
from astropy.io import fits
import os
import yaml
import numpy as np
from subprocess import call

SERVER = 'plhstins2.stsci.edu'
DATA_DIR = '/ifs/archive/ops/hst/public/'


def getdata(ipppssoot, raw=False):
    """Grab darks from MAST"""
    print('Retrieving {}...'.format(ipppssoot.lower()))
    if len(ipppssoot) != 9:
        raise IOError('IPPPSSOOT not proper length:  {}'.format(ipppssoot))

    ippp = ipppssoot[0:4].lower()

    if raw:
        for rawtype in ['raw', 'flt']:
            os.system('scp {}:"{}{}/{}/*_{}.fits" data/'.format(
                SERVER, DATA_DIR, ippp, ipppssoot.lower(), rawtype))
    else:
        for filetype in ['flt']:
            os.system('scp {}:"{}{}/{}/*_{}.fits" data/'.format(SERVER, DATA_DIR, ippp, ipppssoot.lower(), filetype))
    print('')


if __name__ == '__main__':
    source = "ocqq9tq6q_flt.fits"
    # Select Source Dark
    dark = fits.open(os.path.join("data/sources/", source))

    # Grab Keyword Values
    print(dark[1].header['OCCDHTAV'])
    print(dark[0].header['TDATEOBS'])
    reffile = dark[0].header['DARKFILE'].split('$')[-1]

    # Grab Basedark Files (anneal period)
    with open('data/anneal_boundaries_v2.yml') as f:
        boundaries = yaml.load(f)

    dataset = source.split('_')[0].upper()
    anneal_mask = [np.any(np.array([dataset in d.values() for d in boundaries[i]['darks']]))
                   for i in range(len(boundaries))]
    basedark_files = [d.get('exposure').lower() for d in np.array(boundaries)[anneal_mask][0]['darks']]
    for ipppssoot in basedark_files:
        try:
            getdata(ipppssoot, raw=False)
        except IOError as e:
            print(e)
            print('')
            continue
    import pdb; pdb.set_trace()
    call('rm data/*.fits', shell=True)  # Delete all basedark inputs
    import pdb; pdb.set_trace()
    # Grab Weekdark Files
    reffile = '/grp/hst/cdbs/oref/'+reffile
    weekdark_files = str(fits.open(reffile)[0].header['HISTORY']).split("The following input files were used:")[-1].split("_flt.fits")[:-1]

    # Grab darks used in generation of the source darks weekdark ref file (dark dark dark dark...)
    for ipppssoot in weekdark_files:
        try:
            getdata(ipppssoot, raw=False)
        except IOError as e:
            print(e)
            print('')
            continue



