''' Copied from Phil Hodge's itools, 2017 Oct 23.
    Designed to support the PyRAF msarith function.
'''

import glob
import os
import re
import numpy as np
import astropy.io.fits as fits

# This file contains utility functions for Python modules to replace IRAF
# tasks.

def splitOnComma(images):
    """Split an input string on commas, ignoring commas within brackets.

    This function also expands environment variables and wildcards in the
    image names.

    Parameters
    ----------
    images: str
        One or more comma-separated image names, possibly including an
        extension number (e.g. [1]) and/or extension name and version
        (e.g. [sci,3]).

    Returns
    -------
    list of strings
        Each string in the list is an image name.
    """

    MAX_COUNT = 100

    if isinstance(images, str):
        temp = []
        i = 0
        in_brackets = False
        for j in range(len(images)):
            if not in_brackets and images[j] == '[':
                in_brackets = True
            elif in_brackets and images[j] == ']':
                in_brackets = False
            if images[j] == ',' and not in_brackets:
                temp.append(images[i:j].strip())
                i = j + 1
        if not images.endswith(","):
            temp.append(images[i:].strip())
    elif isinstance(images, list) or isinstance(images, tuple):
        temp = images
    else:
        temp = [images]

    # Expand environment variables and wildcards.
    temp2 = []
    for image in temp:
        if isinstance(image, str):
            for i in range(MAX_COUNT):
                exp_image = os.path.expandvars(image)
                if exp_image == image:
                    break
                image = exp_image
            # Look for expressions in brackets at the end of the file name.
            len_image = len(image)
            k = len_image
            done = False
            while not done:
                if image[:k].endswith("]"):
                    n = image[0:k].rfind("[")
                    if n > 0:
                        k = n
                    else:
                        done = True
                else:
                    done = True
            if k < len_image:
                save = image[k:]
            else:
                save = ""
            if "[" in image[0:k] or "*" in image[0:k] or "?" in image[0:k]:
                # Expand wildcards.
                all = glob.glob(image[0:k])
                all.sort()
                for im in all:
                    temp2.append(im + save)
            else:
                temp2.append(image)             # no wildcard character
        else:
            temp2.append(image)

    return temp2

def getNameEtc(image):
    """Split an image name into a file name, extension, etc.

    Parameters
    ----------
    image: str
        An image name, i.e. a file name, possibly with one or more
        expressions in brackets giving an extension number (e.g. [1]) or
        extension name and/or version (e.g. [sci,3]).  There may also be
        brackets containing qualifiers such as append=yes or inplace+.

    Returns
    -------
    (filename, extension, section, qualifier): tuple
        `filename` is the name of the file.
        `extension` may be an integer (the extension number) or a tuple;
        if a tuple, the first element will be a string (the EXTNAME value),
        and if there is a second element it will be an integer (the
        EXTVER value).  `extension` will be None if it wasn't specified.
        `section` is a string containing the image section in numpy
        notation, or None if not specified.
        `qualifier` will be a dictionary containing keywords specified
        by the user and True or False values, or None if not specified.
        The keywords may be "append" or "inplace".
    """

    # If specified, the extension, image section, and qualifiers are
    # enclosed in brackets.
    temp_parts = image.split("[")
    parts = []
    for part in temp_parts:
        part = part.replace("]", "")
        parts.append(part)
    if not parts:
        raise RuntimeError("don't understand image name '%s'" % image)

    nparts = len(parts)

    extension = None
    python_section = None
    qualifier = None
    got_extension = False
    got_section = False
    got_qualifier = False
    garbled = False
    filename = parts[0]

    # If an extension was specified, assume it will be the first
    # bracketed expression.
    k = 1
    if nparts > k:
        (extension, got_extension) = checkExtension(parts[k])
        if got_extension:
            k += 1

    if nparts > k:
        (python_section, got_section, garbled) = checkSection(parts[k])
        if got_section:
            k += 1
        if garbled:
            raise RuntimeError("don't understand image name '%s'" % image)

    done = False
    qualifier = {}
    while not done:
        if nparts > k:
            (q, got_qualifier, garbled) = checkQualifier(parts[k])
            if got_qualifier:
                qualifier.update(q)
                k += 1
            else:
                done = True
            if garbled:
                raise RuntimeError("don't understand image name '%s'" % image)
        else:
            done = True

    # Check that all bracketed expressions have been examined.
    if nparts > k:
        if nparts - k > 1:
            msg = "can't interpret these parts of image name: "
        else:
            msg = "can't interpret this part of image name: "
        for i in range(k, nparts):
            msg = msg + " [%s]" % parts[i]
        raise RuntimeError(msg)

    return (filename, extension, python_section, qualifier)

def checkInt(charval):
    """Test whether the input string is an integer value.

    Parameters
    ----------
    charval: str
        A string that might contain an integer value.  The value will also
        be accepted as an integer if it is floating point but with no
        fractional part.

    Returns
    -------
        If the value in `charval` is an integer, the integer value will
        be returned; otherwise, None will be returned.
    """

    try:
        intval = int(charval)
    except ValueError:
        intval = None

    return intval

def checkExtension(part):
    """Extract an extension specification from the input string.

    Parameters
    ----------
    part: str
        A string that may or may not contain a specification of the
        extension of a FITS file.  For example:
            "3"
            "extension=3"
            "sci,1"
            "'sci',1"
            "extname=sci,extver=1"
            "extname='sci',extver=1"

    Returns
    -------
    (extension, got_extension): tuple
        If the current part of an input image name does actually contain
        an extension specification in a recognized format, `extension` will
        be either an integer or a tuple, and `got_extension` will be True.
        Otherwise, `extension` will be None and `got_extension` will be
        False.
        If `extension` is an integer, it will be the extension number (0
        is the primary HDU).  If it is a tuple, it will contain one or two
        elements; the first element will be a string, the extension name
        (extname), and the second element (if present) will be an integer,
        the extension version number (extver).
    """

    extension = None
    got_extension = False

    part = part.lower()

    test_parts = part.split(",")
    nwords = len(test_parts)
    if nwords == 1:
        extn = checkExtn(test_parts[0], "extension")
        if extn is not None:
            extension = extn
            got_extension = True
        else:
            extname = checkExtn(test_parts[0], "extname")
            if extname is not None:
                extension = (extname,1)
                got_extension = True
            else:
                extension = None
    elif nwords == 2:
        extname = checkExtn(test_parts[0], "extname")
        extver = checkExtn(test_parts[1], "extver")
        if extname is not None and extver is not None:
            extension = (extname,extver)
            got_extension = True

    return (extension, got_extension)

def checkExtn(word, keyword):
    """Extract an extension specification from the input string.

    Parameters
    ----------
    word: str
        A string that may or may not contain a specification of the
        extension of a FITS file.  For example:
            "3"
            "extension=3"
            "sci,1"
            "'sci',1"
            "extname=sci,extver=1"
            "extname='sci',extver=1"
    keyword: str
        A name to look for in `word`.  `keyword` should be "extension",
        "extname", or "extver".

    Returns
    -------
    value: int or str
        An extension specification.  If it is an integer, it will be
        either the extension number or extension version (extver) number.
        If it is a string, it will be the extension name (extname).
    """

    value = None
    try:
        value = int(word)
    except ValueError:
        locn_equals = word.rfind("=")
        if keyword == "extname":
            if locn_equals >= 0:
                if locn_equals <= len(keyword) and \
                   word[0:locn_equals] == keyword[0:locn_equals]:
                    value = word[locn_equals+1:]
                else:
                    value = None
            else:
                value = word
            if value is not None:
                value = value.replace("'", "")
                value = value.replace('"', '')
        else:
            if locn_equals >= 0 and locn_equals <= len(keyword) and \
               word[0:locn_equals] == keyword[0:locn_equals]:
                value = int(word[locn_equals+1:])

    return value

def checkSection(part):
    """Test whether the input string is an image section.

    Parameters
    ----------
    part: str
        A string that may or may not contain an image section.

    Returns
    -------
    (python_section, got_section, garbled): tuple
        `python_section` is the image section converted to Python syntax.
        `got_section` will be True if the input does contain an image
        section, False otherwise.  `garbled` will be True if the input
        looks something like an image section (e.g. contains ":" or "*")
        but doesn't really match what was expected, and it will be False
        if the input does actually look like an image section.
    """

    python_section = ""
    got_section = False
    garbled = False

    # quick check:  an '=' sign would not be found in an image section
    if part.rfind("=") >= 0:
        return (python_section, got_section, garbled)

    part = part.lower()

    is_section = False

    # Split into one string per dimension, and check whether it looks
    # like an image section.
    sections = part.split(",")
    ndim = len(sections)
    for section in sections:
        test_int = checkInt(section)
        if test_int is not None:
            is_section = True
        else:
            if section.find(":") >= 0 or section.find("*") >= 0:
                is_section = True
                test_words = section.split(":")
                for word in test_words:
                    test_int = checkInt(word)
                    if word != "*" and test_int is None:
                        # looks something like an image section, but not quite
                        garbled = True
                        is_section = False
            else:
                is_section = False
        if not is_section:
            break

    if is_section:
        # It does look like an image section; convert to Python syntax.
        python_section = "["
        connector = ""                  # will be changed at end of loop
        for i in range(ndim-1, -1, -1):
            if sections[i] == "*":
                temp_section = ":"
            elif sections[i].find(":") >= 0:
                words = sections[i].split(":")
                nwords = len(words)
                if nwords < 2 or nwords > 3:
                    garbled = True
                    break
                if nwords == 2:
                    temp_section = "%d:%d" % (int(words[0]) - 1,
                                              int(words[1]))
                else:
                    temp_section = "%d:%d:%d" % (int(words[0]) - 1,
                                                 int(words[1]),
                                                 int(words[2]))
            else:
                test_int = checkInt(sections[i])
                if test_int is None:
                    garbled = True
                    break
                else:
                    temp_section = str(test_int - 1)
            python_section = python_section + connector + temp_section
            connector = ","
        if garbled:
            python_section = None
        else:
            python_section = python_section + "]"
            got_section = True
    else:
        python_section = None

    return (python_section, got_section, garbled)

def checkQualifier(part):
    """Test whether the input string is a qualifier for an image.

    Parameters
    ----------
    part: str
        A string that may or may not contain a qualifier, such as
        "append=True" or "inplace=True".

    Returns
    -------
    (qualifier, got_qualifier, garbled): tuple
        `qualifier` is a dictionary; if the input was not in a format
        recognized as a qualifier for an image, `qualifier` will be empty
        and `got_qualifier` will be False.  If the input does look like a
        qualifier (contains "=", "+", "-", "True", "False", "yes", "no",
        "append" or "inplace"), then the key for `qualifier` will be the
        string (e.g. "append" or "inplace"), and the value will be either
        True or False.
        If the input looks something like a qualifier but doesn't quite
        match the expected pattern, `garbled` will be set to True.
    """

    qualifier = {}
    got_qualifier = False
    garbled = False

    part = part.lower()

    qualifiers = part.split(",")
    for qual in qualifiers:
        is_qualifier = False

        # First find the value (which should be equivalent to True or False).
        locn_equals = qual.rfind("=")
        locn_plus = qual.rfind("+")
        locn_minus = qual.rfind("-")
        if locn_equals >= 0:
            locn_break = locn_equals
            if qual.find("yes") >= 0 or qual.find("true") >= 0:
                value = True
                is_qualifier = True
            elif qual.find("no") >= 0 or qual.find("false") >= 0:
                value = False
                is_qualifier = True
            else:
                garbled = True
        elif qual.endswith("+"):
            locn_break = locn_plus
            value = True
            is_qualifier = True
        elif qual.endswith("-"):
            locn_break = locn_minus
            value = False
            is_qualifier = True
        else:
            is_qualifier = False
            break

        # Now check whether we recognize the keyword.
        key = None
        if is_qualifier:
            checklist = ["append", "inplace"]
            for keyword in checklist:
                if locn_break <= len(keyword) and \
                   qual[0:locn_break] == keyword[0:locn_break]:
                    key = keyword
                    break

        if key is not None:
            qualifier[key] = value
            got_qualifier = True

    return (qualifier, got_qualifier, garbled)

def getLtmLtv(hdr):
    """Get all the LTMi_j and LTVj keywords from the header.

    Parameters
    ----------
    hdr: pyfits Header object
        Header from an input image.  Keywords of the form LTMi_j and LTVj
        will be read, and their values will be returned in arrays ltm and
        ltv respectively.

    Returns
    -------
    (ltm, ltv): tuple
        ltm is the 2-D array of LTMi_j keyword values, and ltv is the 1-D
        array of LTVj keyword values.  For LTMi_j, the values are stored
        as follows (for example):
            [[ LTM1_1    LTM1_2]
             [ LTM2_1    LTM2_2]]
            >>> for j in range(2):
            ...     for i in range(2):
            ...         print(j, i, ltm[j,i])
            0 0 11.0    # LTM1_1
            0 1 12.0    # LTM1_2
            1 0 21.0    # LTM2_1
            1 1 22.0    # LTM2_2
    """

    min_i = 9999
    max_i = -1
    min_j = 9999
    max_j = -1
    found_it = False

    # key is an LTM or LTV header keyword, value is a tuple giving
    # the indices (minus one) from the header keyword and the header value.
    ltm_dict = {}
    ltv_dict = {}
    r1 = re.compile(r"""(?x)
                    LTM(?P<index1>\d+)_(?P<index2>\d+)
                    """)
    r2 = re.compile(r"""(?x)
                    LTV(?P<index1>\d+)
                    """)
    for (keyword, value) in hdr.items():
        m1 = r1.match(keyword)
        m2 = r2.match(keyword)
        if m1:
            i = int(m1.group("index1")) - 1
            j = int(m1.group("index2")) - 1
            ltm_dict[keyword] = (i, j, value)
            min_i = min(min_i, i)
            max_i = max(max_i, i)
            min_j = min(min_j, j)
            max_j = max(max_j, j)
            found_it = True
        if m2:
            j = int(m2.group("index1")) - 1
            ltv_dict[keyword] = (j, value)
            min_j = min(min_j, j)
            max_j = max(max_j, j)
            found_it = True

    if not found_it:
        return (None, None)
    if min_i < 0 or min_j < 0:
        print("Warning:  Some LTM or LTV keywords have invalid indices")
        return (None, None)

    # The LTM matrix should be square.
    max_j = max(max_j, max_i)

    ltm = np.identity(max_j+1, dtype=np.float64)
    ltv = np.zeros(max_j+1, dtype=np.float64)
    if ltm_dict:
        for (keyword, keyval) in ltm_dict.items():
            (i, j, value) = keyval
            ltm[i,j] = value
    if ltv_dict:
        for (keyword, keyval) in ltv_dict.items():
            (j, value) = keyval
            ltv[j] = value

    return (ltm, ltv)

def writeLtmLtv(hdr, ltm, ltv):
    """Write all the LTMi_j and LTVj keywords to the header.

    Parameters
    ----------
    hdr: pyfits Header object
        Header from an input image.  Keywords of the form LTMi_j and LTVj
        will be read, and their values will be returned in arrays ltm and
        ltv respectively.

    ltm: array_like, or None
        The array (matrix) of LTMi_j keyword values, or None.

    ltv: array_like, or None
        The array (vector) of LTVj keyword values, or None.
    """

    if ltm is not None:
        for j in range(ltm.shape[1]):
            for i in range(ltm.shape[0]):
                keyword = "LTM%d_%d" % (i+1, j+1)
                hdr[keyword] = ltm[i,j]

    if ltv is not None:
        for j in range(ltv.shape[0]):
            keyword = "LTV%d" % (j+1,)
            hdr[keyword] = ltv[j]

def imageSectionLtmLtv(python_section):
    """Get the ltm and ltv corresponding to an image section.

    Parameters
    ----------
    python_section: str
        A string containing an image section in numpy slice notation, e.g.
        "[500:520,:]", or None.

    Returns
    -------
    (ltm, ltv): tuple
        ltm is the 2-D array of LTMi_j keyword values, and ltv is the 1-D
        array of LTVj keyword values.  See getLtmLtv for further details.
        The returned tuple will be (None, None) if python_section is None.
    """

    if python_section is None:
        return (None, None)

    np_section = python_section.replace("[", "")
    np_section = np_section.replace("]", "")
    sections = np_section.split(",")
    nelem = len(sections)
    start = []
    incr = []
    ltm = np.identity(nelem, dtype=np.float64)
    ltv = np.zeros(nelem, dtype=np.float64)
    for section in sections:
        if section == ":":
            start.append(1.)
            incr.append(1.)
        else:
            words = section.split(":")
            nwords = len(words)
            if nwords == 0:
                start.append(1.)
                incr.append(1.)
            else:
                if words[0] == "":
                    start.append(1.)
                else:
                    test_int = checkInt(words[0])
                    if test_int is None:
                        start.append(1.)
                    else:
                        start.append(float(test_int) + 1.)
                if nwords >= 3:
                    if words[2] == "":
                        incr.append(1.)
                    else:
                        test_int = checkInt(words[2])
                        if test_int is None:
                            incr.append(1.)
                        else:
                            incr.append(float(test_int))
                else:
                    incr.append(1.)

    j = nelem - 1
    for i in range(nelem):
        ltm[i,i] = 1. / incr[j]
        ltv[i] = 1. - start[j] / incr[j]
        j -= 1

    return (ltm, ltv)

def binningLtmLtv(bin):
    """Get the ltm and ltv corresponding to a list of binning factors.

    Parameters
    ----------
    bin: list of int
        A list of binning factors.  A factor of n means that n pixels will
        be added (or averaged) together to generate one output pixel.

    Returns
    -------
    (ltm, ltv): tuple
        ltm is the 2-D array of LTMi_j keyword values, and ltv is the 1-D
        array of LTVj keyword values.  See getLtmLtv for further details.
        The returned tuple will be the identity array and ltv will be an
        array of zeros if all the bin factors are 1.
    """

    nelem = len(bin)
    ltm = np.identity(nelem, dtype=np.float64)
    ltv = np.zeros(nelem, dtype=np.float64)

    j = nelem - 1
    for i in range(nelem):
        ltm[i,i] = 1. / bin[j]
        ltv[i] = (bin[j] - 1.) / (2. * bin[j])
        j -= 1

    return (ltm, ltv)

def combineLtm(ltm1, ltv1, ltm2, ltv2):
    """Combine (ltm1, ltv1) and (ltm2, ltv2).

    Parameters
    ----------
    ltm1: 2-D numpy array
        First 2-D array of LTMi_j values.

    ltv1: 1-D numpy array
        First 1-D array of LTVj values.

    ltm2: 2-D numpy array
        Second 2-D array of LTMi_j values.

    ltv2: 1-D numpy array
        Second 1-D array of LTVj values.

    Returns
    -------
    (ltm, ltv): tuple
        The result of applying (ltm2, ltv2) to (ltm1, ltv1).  In matrix
        notation:

            im1 = ltm1 * ref + ltv1
            im2 = ltm2 * im1 + ltv2
        or:
            im2 = ltm2 * (ltm1 * ref + ltv1) + ltv2
                = (ltm2 * ltm1) * ref + (ltm2 * ltv1 + ltv2)
        so:
            ltm = ltm2 * ltm1
            ltv = ltm2 * ltv1 + ltv2

        where ref is a vector of pixel coordinates in the reference image
        (before taking an image section [slice] or binning), im1 is a
        vector of pixel coordinates in the first image (may include a
        slice), and im2 is a vector of pixel coordinates in the second
        image or after binning.

        ltm is the 2-D array of LTMi_j keyword values, and ltv is the 1-D
        array of LTVj keyword values.  See getLtmLtv for further details.
        The returned tuple will be (None, None) if all the input variables
        are None.
    """

    nelem = None
    if ltm1 is not None:
        nelem = len(ltm1)
    elif ltv1 is not None:
        nelem = len(ltv1)
    elif ltm2 is not None:
        nelem = len(ltm2)
    elif ltv2 is not None:
        nelem = len(ltv2)

    if nelem is None:
        return (None, None)

    if ltm1 is None:
        ltm1 = np.identity(nelem, dtype=np.float64)
    if ltv1 is None:
        ltv1 = np.zeros(nelem, dtype=np.float64)
    if ltm2 is None:
        ltm2 = np.identity(nelem, dtype=np.float64)
    if ltv2 is None:
        ltv2 = np.zeros(nelem, dtype=np.float64)

    ltm1 = np.matrix(ltm1)
    ltv1 = np.matrix(ltv1)
    ltm2 = np.matrix(ltm2)
    ltv2 = np.matrix(ltv2)
    ltm = ltm2 * ltm1
    ltv = ltm2 * ltv1.T + ltv2.T

    ltm = ltm.A
    ltv = ltv.T
    ltv = ltv.A[0]

    return (ltm, ltv)

def extnameExtver(hdr, extension):
    """Update or assign extname and extver.

    Parameters
    ----------
    input: hdr
        A pyfits Header object, which may be update in-place.

    extension: int or tuple
        If `extension` is an integer, it is interpreted as the extension
        number, and no change will be made to `hdr`.
        If `extension` is a tuple, it is interpreted as containing the
        extension name and/or version number, and it should have only one
        or two elements.  If an element is a string, keyword EXTNAME in
        the header will be assigned that value.  If an element is an
        integer, keyword EXTVER in the header will be assigned that value.
    """

    if extension is not None:
        if isinstance(extension, int):
            return
        for i in range(len(extension)):
            if isinstance(extension[i], str):
                hdr["extname"] = extension[i]
            elif isinstance(extension[i], int):
                hdr["extver"] = extension[i]
            else:
                print("Warning:  don't know what to do with '%s'" % \
                      repr(extension[i]))

def writeOutput(out_filename, out_extension, out_section,
                ltm, ltv, input_is_image, in_phdr, in_hdr,
                data, pixtype=None):
    """Write the output data to a FITS file.

    Parameters
    ----------
    out_filename: str
        The name of the output file.

    out_extension: int or str or tuple
        A specification of the FITS extension for the output image.

    out_section: str, or None
        An image section for the output image.  This may be given if
        writing to an existing image; only this portion will be modified.

    ltm: array_like, or None
        The array (matrix) of LTMi_j keyword values, or None.

    ltv: array_like, or None
        The array (vector) of LTVj keyword values, or None.

    in_phdr: pyfits Header object
        Primary header from an input file to use as a template when (if)
        creating a new output file.

    in_hdr: pyfits Header object
        Header from the input image to use as a template when creating a
        new output header/data unit.

    data: numpy array
        Image data to be written to the output file.

    pixtype: numpy data type, or None
        Data type to use when writing the data to the output file.  If
        `pixtype` is None, the data type of `data` will not be modified.
    """

    exists = os.access(out_filename, os.W_OK)
    if exists:
        ofd = fits.open(out_filename, mode="update")
        try:
            out_hdu = ofd[out_extension]
        except KeyError:
            out_hdu = None
    else:
        out_hdu = None

    if pixtype is not None:
        data = data.astype(pixtype)

    if exists:
        if out_hdu is None:                     # not found in output file
            if input_is_image:
                out_hdu = fits.ImageHDU(header=in_hdr, data=data)
                writeLtmLtv(out_hdu.header, ltm, ltv)
            else:
                out_hdu = fits.BintableHDU(header=in_hdr, data=data)
            extnameExtver(out_hdu.header, out_extension)
            ofd.append(out_hdu)
        else:
            if out_section is None:
                if input_is_image:
                    out_hdu.header = fits.ImageHDU(header=in_hdr).header
                    writeLtmLtv(out_hdu.header, ltm, ltv)
                else:
                    out_hdu.header = fits.BintableHDU(header=in_hdr).header
                out_hdu.data = data
                extnameExtver(out_hdu.header, out_extension)
            else:
                # Since an image section was specified for the output
                # image, only modify the data, not the header.
                out_data = eval("out_hdu.data" + out_section)
                out_data[:] = data.copy()
        ofd.close()
    else:
        if out_extension == 0 or out_extension is None:
            # We can put an image in the primary HDU, but if we're
            # writing a table it must go in an extension.
            if input_is_image:
                primary_hdu = fits.PrimaryHDU(header=in_phdr, data=data)
                ofd = fits.HDUList(primary_hdu)
                writeLtmLtv(ofd[0].header, ltm, ltv)
            else:
                primary_hdu = fits.PrimaryHDU(header=in_phdr)
                ofd = fits.HDUList(primary_hdu)
                ofd.append(in_hdu)
        else:
            # An output extension was specified.
            primary_hdu = fits.PrimaryHDU(header=in_phdr)
            ofd = fits.HDUList(primary_hdu)
            if input_is_image:
                out_hdu = fits.ImageHDU(header=in_hdr, data=data)
                writeLtmLtv(out_hdu.header, ltm, ltv)
            else:
                out_hdu = fits.BintableHDU(header=in_hdr, data=data)
            extnameExtver(out_hdu.header, out_extension)
            ofd.append(out_hdu)
        ofd.writeto(out_filename)
