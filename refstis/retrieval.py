import os
try:
    from urllib.request import urlopen
    from urllib.parse import urlparse
except ImportError:
    from urlparse import urlparse
    from urllib import urlopen
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
                    archiveUserEmail=settings['email'],
                    ftpHost=settings['host'],
                    ftpDir=settings['retrieve_directory'],
                    ftpUser=settings['ftp_user'],
                    datasets=datasets,
                    suffix="<suffix name=\"*\" />")

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
    params = urlparse.urlencode({'dadshost': settins['dads_host'],
                                'dadsport': 4703,
                                'mission':'HST',
                                'request': request_xml_str})
    headers = {"Accept": "text/html", "User-Agent":"{0}PythonScript".format(user)}
    req = httplib.HTTPSConnection(settings['archive'])
    req.request("POST", "/cgi-bin/dads.cgi", params, headers)
    response = req.getresponse().read()
    req.close()

    return response

#-------------------------------------------------------------------------------
