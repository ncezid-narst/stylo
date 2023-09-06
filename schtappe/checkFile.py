# Joseph S. Wirth
# September 2023

import glob, os


def __getExpectedFiles(fn:str) -> list[str]:
    """gets a list of expected files from an input file containing the filenaems

    Args:
        fn (str): the input file

    Returns:
        list[str]: a list of the filenames contained within it
    """
    # constants
    SEP = "\t"
    FILE_IDX = 1
    
    # initialize output
    outL = list()
    
    # go through each line of the file; skip headers
    header = True
    with open(fn, 'r') as fh:
        for line in fh:
            if header:
                header = False
            else:
                # save the filename
                outL.append(line.split(SEP)[FILE_IDX])
    
    return outL


def __findMissingFiles(existing:set[str], expected:list[str]) -> list[str]:
    """finds any files that should exist that do not

    Args:
        existing (set[str]): a set of all the files that currently exist
        expected (list[str]): a set of the files that should exist

    Returns:
        list[str]: a list of files that should exist but do not
    """
    # initialize output
    missing = list()
    
    # find any expected files that do not exist
    for file in expected:
        if file not in existing:
            missing.append(file)

    return missing


def main(fn:str, filePattern:str) -> None:
    """main runner function:
         * checks that all files in the filename exist
         * raises an error if any expected files do not exist

    Args:
        fn (str): the filename containing the expected files
        filePattern (str): a glob pattern for all the existing files

    Raises:
        FileNotFoundError: one or more expected files do not exist
    """
    # constants
    ERR_MSG_A = 'the following files were not found in the '
    ERR_MSG_B = "directory:\n\n"
    
    # get the existing files and the expected files
    existingFiles = {os.path.basename(x) for x in glob.glob(filePattern)}
    expectedFiles = __getExpectedFiles(fn)
    
    missingFiles = __findMissingFiles(existingFiles, expectedFiles)
    
    if missingFiles != []:
        raise FileNotFoundError(ERR_MSG_A + os.path.dirname(filePattern) + \
                                ERR_MSG_B + "\n".join(missingFiles))
