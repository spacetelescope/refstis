import os
try:
    from urllib.request import urlopen
    from urllib.parse import urlparse, urlencode
    from http.client import HTTPSConnection
except ImportError:
    from urlparse import urlparse
    from urllib import urlopen, urlencode
    from httplib import HTTPSConnection
import string

from .SignStsciRequest import SignStsciRequest

REQUEST_TEMPLATE = string.Template('\
<?xml version=\"1.0\"?> \n \
<!DOCTYPE distributionRequest SYSTEM \"http://dmswww.stsci.edu/dtd/sso/distribution.dtd\"> \n \
  <distributionRequest> \n \
    <head> \n \
      <requester userId = \"$archive_user\" email = \"$email\" source = "starview"/> \n \
      <delivery> \n \
        <ftp hostName = \"$ftp_host\" loginName = \"$ftp_user\" directory = \"$ftp_dir\" secure = \"true\" /> \n \
      </delivery> \n \
      <process compression = \"none\"/> \n \
    </head> \n \
    <body> \n \
      <include> \n \
        <select> \n \
          $suffix\n \
        </select> \n \
        $datasets \n \
      </include> \n \
    </body> \n \
  </distributionRequest> \n' )

# ------------------------------------------------------------------------------

def everything_retrieved(tracking_id):
    '''Check request status page to see if retreival is done.

    Parameters
    ----------
    tracking_id : string
        A submission ID string.

    Returns
    -------
    done : bool
        Boolean specifying is data is retrieved or not.
    killed : bool
        Boolean specifying is request was killed.
    '''

    done = False
    killed = False
    status_url = "http://archive.stsci.edu/cgi-bin/reqstat?reqnum=={0}".format(tracking_id)
    for line in urlopen(status_url).readlines():

        if isinstance(line, bytes):
            line = line.decode()

        if "State" in line:
            if "COMPLETE" in line:
                done = True
            elif "KILLED" in line:
                killed = True
    return done, killed

# ------------------------------------------------------------------------------

def build_xml_request(datasets, settings):
    """Build the XML request for the given datasets

    Parameters
    ----------
    datasets : list
        A list of rootnames to download from MAST.

    Returns
    -------
    xml_request : string
        The XML request string.
    """

    # The list of datasets must be one string
    datasets = ''.join(['<rootname>{0}</rootname>\n'.format(rootname) for rootname in datasets])

    request_string = REQUEST_TEMPLATE.safe_substitute(
                                            archive_user=settings['archive_user'],
                                            email=settings['email'],
                                            ftp_host=settings['host'],
                                            ftp_dir=settings['retrieve_directory'],
                                            ftp_user=settings['ftp_user'],
                                            datasets=datasets,
                                            suffix="<suffix name=\"RAW\" />")

    xml_request = string.Template(request_string)
    xml_request = xml_request.template

    return xml_request

#-------------------------------------------------------------------------------

def submit_xml_request(xml_request, settings):
    """Submit the XML request to the MAST archive.

    Parameters
    ----------
    xml_request : string
        The request XML string.

    Returns
    -------
    submission_results : httplib object
        The XML request submission results.
    """

    user = os.environ.get("USER")
    home = os.environ.get("HOME")

    signer = SignStsciRequest()
    request_xml_str = signer.signRequest('{0}/.ssh/privkey.pem'.format(home), xml_request)
    params = urlencode({'dadshost': settings['dads_host'],
                        'dadsport': 4703,
                        'mission':'HST',
                        'request': request_xml_str})
    headers = {"Accept": "text/html", "User-Agent":"{0}PythonScript".format(user)}
    req = HTTPSConnection(settings['archive'])
    req.request("POST", "/cgi-bin/dads.cgi", params, headers)
    response = req.getresponse().read()
    req.close()

    return response

#-------------------------------------------------------------------------------
