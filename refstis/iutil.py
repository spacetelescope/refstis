''' Copied from Phil Hodge's itools, 2017 Oct 23.
    Designed to support the PyRAF msarith function.
'''

import glob
import os

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
