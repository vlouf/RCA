import os
import datetime

import numpy as np


def get_files(inpath, date=None):
    '''
    Find the list of files with the supported extension in the given
    path. Will recursively search in subdirectories too. If provided a date
    (string or datetime object) it will only returns the files whose
    filename matches.

    Parameters
    ==========
        inpath: str
            General path for the data.
        date: str or datetime
            Look for files with a specific date. If date is None, it will look
            for all files with the supported extension.

    Returns
    =======
        flist: list
            List of files found.
    '''

    supported_extension = ['.nc', '.NC']
    flist = []

    # Check date type
    if isinstance(date, datetime.datetime):
        date = date.strftime("%Y%m%d")

    for dirpath, dirnames, filenames in os.walk(inpath):
        for filenames_slice in filenames:

            # If no date provided, nothing new under the sun
            if date is None:
                pass  # pretends there was no if statement
            elif date in filenames_slice:
                pass  # pretends there was no if statement
            else:
                continue  # A date was given and we didn't found it.

            file_extension = os.path.splitext(str(filenames_slice))[1]
            # Get extension

            if np.any(np.in1d(supported_extension, file_extension)):
                # Check if file extension is in the list of supported ones
                the_path = os.path.join(dirpath, filenames_slice)
            else:  # If not test next file.
                continue

            # File does have the supported extension, we keep it for returning
            # list
            flist.append(the_path)

    to_return = flist

    return sorted(to_return)  # Type: List[str, ...]
