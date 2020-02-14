''' Copied from Phil Hodge's itools, 2017 Oct 23.
    Designed to emulate the PyRAF msarith function.
'''

import math
import os
import numpy as np
import astropy.io.fits as fits
from . import iutil


def msarith(operand1, op, operand2, result, divzero=0., verbose=True):
    """Combine input files to one output file.

    This is a partial emulation of the STSDAS msarith task.
    Not supported:
        image sections
        SAMP and TIME extensions
        explicit lists of which image sets should be operated on
        different behavior for count rate vs counts

    Parameters
    ----------
    operand1: str or list of str
        One or more comma-separated input file names, or a list of input
        file names, or a numerical value.  A numerical value may be given
        either as a single number or as a two-element tuple, the value and
        its error estimate.

    op: str
        The operation to be performed, "+", "-", "*", or "/".

    operand2: str or list of str
        One or more comma-separated input file names, or a list of input
        file names, or a numerical value.  See `operand1`.

    result: str or list of str
        One or more comma-separated output file names, or a list of output
        file names.  If both `operand1` and `operand2` are numerical
        values, `result` will be ignored and the computed result will be
        printed (if verbose) and returned.

    divzero: float
        Replacement value for division by zero.

    verbose: boolean
        If True, print info.

    Returns
    -------
    value or None
        If both `operand1` and `operand2` are numerical values, the
        computed value and its error estimate will be returned as a
        tuple.  If either operand is a file name, None will be returned.
    """

    op = op.strip()

    # Check for numerical constant.
    (operand1, op1_const) = isNumConst(operand1)
    (operand2, op2_const) = isNumConst(operand2)
    out_const = (op1_const and op2_const)

    if op1_const:
        len_op1 = 1
    else:
        operand1 = iutil.splitOnComma(operand1)
        len_op1 = len(operand1)
    if op2_const:
        len_op2 = 1
    else:
        operand2 = iutil.splitOnComma(operand2)
        len_op2 = len(operand2)
    if out_const:
        len_out = 1
    else:
        result = iutil.splitOnComma(result)
        len_out = len(result)

    bad_length = False
    if len_out == 1:
        if len_op1 > 1 or len_op2 > 1:
            bad_length = True
    else:
        if len_op1 > 1 and len_op1 != len_out:
            bad_length = True
        if len_op2 > 1 and len_op2 != len_out:
            bad_length = True
    if bad_length:
        raise RuntimeError("Lengths of input and output lists "
                           "must be the same.")

    for i in range(len_out):
        if op1_const:
            input1 = operand1
        elif len_op1 == 1:
            input1 = operand1[0]
        else:
            input1 = operand1[i]
        if op2_const:
            input2 = operand2
        elif len_op2 == 1:
            input2 = operand2[0]
        else:
            input2 = operand2[i]
        if out_const:
            output = None
        else:
            output = result[i]

        result_dict = oneFileArith(input1, op1_const, op, input2, op2_const,
                                   output, divzero)
        value = result_dict["value"]
        if verbose:
            if value is None:
                if result_dict["status"] == 0:
                    print("%s, %s --> %s" % (str(input1), str(input2), output))
            else:
                print("%.15g, %.15g)" % (value[0], value[1]))

    return value

def isNumConst(operand):
    """Check for a numerical constant.

    Parameters
    ----------
    operand: str or list of str, or numerical value or a two-element tuple
        `operand` may be a file name or a list of file names, or it may be
        a numerical constant (float), or it may be a two-element tuple
        with a numerical value and its error estimate.

    Returns
    -------
    tuple
        (value, flag)
        `flag` is a boolean flag, True if `operand` is a numerical constant
        (or a constant with an error estimate), False if `operand` is a
        file name or a list of file names.
        If `flag` is True, `value` will be a two-element tuple of the
        numerical constant and its error estimate.  If `flag` is False,
        `value` will be a string or list of strings, the name(s) of the
        input files.
    """

    if isinstance(operand, float):
        input = (operand, 0.)
        flag = True
    elif isinstance(operand, int):
        input = (float(operand), 0.)
        flag = True
    elif isinstance(operand, str):
        try:
            # a string containing a number?
            value = float(operand)
            flag = True
            input = (value, 0.)
        except ValueError:
            input = operand
            flag = False
    else:
        # The remaining possibilities are that `operand` could be a two-
        # element tuple of numerical values, or it could be a list of
        # file names.
        try:
            input = [float(operand[0]), float(operand[1])]
            flag = True
        except ValueError:
            input = operand
            flag = False

    return (input, flag)

def oneFileArith(input1, op1_const, op, input2, op2_const, output, divzero):
    """Do image arithmetic on one set of input/output files.

    Parameters
    ----------
    input1: str or two-element tuple
        Either the name of a file, or a tuple giving a numerical constant
        and its error estimate.

    op1_const: bool
        True if `input1` is a numerical constant.

    op: str
        The operation to be performed, "+", "-", "*", or "/".  For example,
        if op = "-", the operation would be to subtract input2 from input1.

    input2: str or two-element tuple
        Either the name of a file, or a tuple giving a numerical constant
        and its error estimate.

    op2_const: bool
        True if `input2` is a numerical constant.

    output: str
        The name of the output file.  If both `input1` and `input2` are
        numerical constants, `output` will be ignored and the result
        will be returned in the function value.

    divzero: float
        If op = "/", `divzero` will be assigned to the output for zero
        pixels in `input2`.

    Returns
    -------
    result_dict: dict
        result_dict["value"] will be None if at least one of the input
        operands was a file name; if both operands were numerical
        constants, the value will be a two-element tuple of the computed
        value and its error estimate.
        result_dict["status"] will be 1 if the files could not be
        processed; 0 is OK.
    """

    status = 0                  # initialize to "OK"

    samp1 = None                # not supported yet
    time1 = None
    samp2 = None
    time2 = None

    if op1_const:
        value1 = input1[0]
        err1 = input1[1]
        dq1 = None
        fd1 = None
        phdu = None
    else:
        fd1 = fits.open(input1)
        phdu = fd1[0]

    if op2_const:
        value2 = input2[0]
        err2 = input2[1]
        dq2 = None
        fd2 = None
    else:
        fd2 = fits.open(input2)
        if phdu is None:
            phdu = fd2[0]

    # If both operands are numeric, compute the value and return.
    if op1_const and op2_const:
        value = constArith(input1, op, input2, divzero)
        return {"value": value, "status": 0}

    # It is assumed that if there are two input files, they are both for
    # the same instrument/detector, or at least that they have the same
    # number of extensions per image set.
    ext_per_imset = getImsetType(phdu.header)

    bad = False
    if op1_const:
        nimsets1 = 1
    else:
        nextend1 = len(fd1) - 1
        nimsets1 = nextend1 // ext_per_imset
        if nextend1 == 0:
            print("%s has no extensions, must be a multiple of %d ..." %
                  (input1, ext_per_imset))
            bad = True
        elif nimsets1 * ext_per_imset != nextend1:
            print("%s has %d extensions, must be divisible by %d ..." %
                  (input1, nextend1, ext_per_imset))
            bad = True
    if op2_const:
        nimsets2 = 1
    else:
        nextend2 = len(fd2) - 1
        nimsets2 = nextend2 // ext_per_imset
        if nextend2 == 0:
            print("%s has no extensions, must be a multiple of %d ..." %
                  (input2, ext_per_imset))
            bad = True
        elif nimsets2 * ext_per_imset != nextend2:
            print("%s has %d extensions, must be divisible by %d ..." %
                  (input2, nextend2, ext_per_imset))
            bad = True
    if bad:
        print(" ... skipping")
        status = 1
        return {"value": None, "status": status}

    phdu.header["filename"] = os.path.basename(output)
    out_fd = fits.HDUList(phdu)

    # Loop over image sets
    for i in range(max(nimsets1, nimsets2)):
        extver = i + 1
        if not op1_const:
            value1 = fd1[("SCI",extver)].data
            err1 = getFloatData(fd1, ("ERR", extver))
            dq1 = getShortData(fd1, ("DQ", extver))
        if not op2_const:
            value2 = fd2[("SCI",extver)].data
            err2 = getFloatData(fd2, ("ERR", extver))
            dq2 = getShortData(fd2, ("DQ", extver))
        (value, err, dq) = imageArith(value1, err1, dq1, op,
                                      value2, err2, dq2, divzero)
        if not op1_const:
            del value1, err1, dq1
        if not op2_const:
            del value2, err2, dq2

        # Construct header/data units for the output file using the above
        # value, err, dq.  Use headers from the first operand if that is
        # a file; otherwise, use headers from the second operand.
        hdu_appended = False
        if not op1_const:
            hdr = fd1[("SCI",extver)].header
            out_fd.append(fits.ImageHDU(data=value, header=hdr, name="SCI"))
            hdr = fd1[("ERR",extver)].header
            out_fd.append(fits.ImageHDU(data=err, header=hdr, name="ERR"))
            hdr = fd1[("DQ",extver)].header
            out_fd.append(fits.ImageHDU(data=dq, header=hdr, name="DQ"))
            if ext_per_imset > 3:
                # just copy without change; this needs to be fixed xxx
                out_fd.append(fd1[("SAMP",extver)])
                out_fd.append(fd1[("TIME",extver)])
            hdu_appended = True

        if not op2_const and not hdu_appended:
            hdr = fd2[("SCI",extver)].header
            out_fd.append(fits.ImageHDU(data=value, header=hdr, name="SCI"))
            hdr = fd2[("ERR",extver)].header
            out_fd.append(fits.ImageHDU(data=err, header=hdr, name="ERR"))
            hdr = fd2[("DQ",extver)].header
            out_fd.append(fits.ImageHDU(data=dq, header=hdr, name="DQ"))
            if ext_per_imset > 3:
                # just copy without change; this needs to be fixed xxx
                out_fd.append(fd2[("SAMP",extver)])
                out_fd.append(fd2[("TIME",extver)])
            hdu_appended = True

    out_fd.writeto(output, overwrite=True)
    out_fd.close()
    if fd1 is not None:
        fd1.close()
    if fd2 is not None:
        fd2.close()

    return {"value": None, "status": status}

def getImsetType(phdr):
    """Find out how many extensions there are per image set.

    The number of extensions in each image set is assumed to be 3 for all
    cases except the following, for which the number is 5:
        INSTRUME = "NICMOS"
        INSTRUME = "WFC3" and DETECTOR = "IR"

    Parameters
    ----------
    phdr: pyfits Header object
        The primary header from one of the input images.

    Returns
    -------
    ext_per_imset: int
        The number of extensions (3 or 5) per image set.
    """

    instrume = phdr.get("instrume", "missing")
    detector = phdr.get("detector", "missing")
    if instrume == "NICMOS":
        ext_per_imset = 5
    elif instrume == "WFC3" and detector == "IR":
        ext_per_imset = 5
    else:
        ext_per_imset = 3

    return ext_per_imset

def getFloatData(fd, extn):
    """Get the data, or return the PIXVALUE.

    Parameters
    ----------
    fd: pyfits HDUList object
        The file handle for one of the input files.

    extn: tuple
        Two-element tuple containing the EXTNAME and EXTVER, for example
        ("ERR", 1).

    Returns
    -------
    numpy array
        If the data for the specified extension exists, then the data array
        will be returned.  Otherwise, a one-element numpy array of type
        numpy.float32 with the value of keyword PIXVALUE will be returned.
    """

    if fd[extn].data is None:
        npix1 = fd[extn].header.get("npix1", 1)
        npix2 = fd[extn].header.get("npix2", 1)
        value = np.zeros((npix2,npix1), dtype=np.float32)
        value[:,:] = fd[extn].header.get("pixvalue", 0.)
    else:
        value = fd[extn].data

    return value

def getShortData(fd, extn):
    """Get the data, or return the PIXVALUE.

    Parameters
    ----------
    fd: pyfits HDUList object
        The file handle for one of the input files.

    extn: tuple
        Two-element tuple containing the EXTNAME and EXTVER, for example
        ("DQ", 1).

    Returns
    -------
    numpy array
        If the data for the specified extension exists, then the data array
        will be returned.  Otherwise, a one-element numpy array of type
        numpy.int16 with the value of keyword PIXVALUE will be returned.
    """

    if fd[extn].data is None:
        npix1 = fd[extn].header.get("npix1", 1)
        npix2 = fd[extn].header.get("npix2", 1)
        value = np.zeros((npix2,npix1), dtype=np.int16)
        value[:,:] = fd[extn].header.get("pixvalue", 0)
    else:
        value = fd[extn].data

    return value

def constArith(input1, op, input2, divzero):
    """Arithmetic for two constants."""

    val1 = input1[0]            # value of the first constant
    err1 = input1[1]            # error estimate for the first constant
    val2 = input2[0]
    err2 = input2[1]

    if op == "+":
        value = val1 + val2
        err = math.sqrt(err1**2 + err2**2)

    elif op == "-":
        value = val1 - val2
        err = math.sqrt(err1**2 + err2**2)

    elif op == "*":
        value = val1 * val2
        err = math.sqrt((val1 * err2)**2 + (val2 * err1)**2)

    elif op == "/":
        if val2 == 0.:
            value = divzero
            err = 0.
        else:
            value = val1 / val2
            a_expr = err1 / val2
            b_expr = err2 * val1 / val2**2
            err = math.sqrt(a_expr**2 + b_expr**2)

    return (value, err)

def imageArith(value1, err1, dq1, op, value2, err2, dq2, divzero):
    """Arithmetic between two images, or an image and a constant."""

    if op == "+":
        value = value1 + value2
        err = np.sqrt(err1**2 + err2**2)

    elif op == "-":
        value = value1 - value2
        err = np.sqrt(err1**2 + err2**2)

    elif op == "*":
        value = value1 * value2
        err = np.sqrt((value1 * err2)**2 + (value2 * err1)**2)

    elif op == "/":
        if isinstance (value2, np.ndarray):
            zero_locn = np.where(value2 == 0.)
            value2[zero_locn] = 1.      # protect against dividing by zero
            value = value1 / value2
            a_expr = err1 / value2
            b_expr = err2 * value1 / value2**2
            err = np.sqrt(a_expr**2 + b_expr**2)
            value[zero_locn] = divzero
            err[zero_locn] = 0.
        elif value2 == 0.:
            value = divzero
            err = 0.
        else:
            value = value1 / value2
            a_expr = err1 / value2
            b_expr = err2 * value1 / value2**2
            err = np.sqrt(a_expr**2 + b_expr**2)

    if dq1 is None:
        dq = dq2
    elif dq2 is None:
        dq = dq1
    else:
        dq = np.bitwise_or(dq1, dq2)

    return (value, err, dq)
