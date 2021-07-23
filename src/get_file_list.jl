"""
    get_file_list(sca, test; [base="/att/projrepo/wfirst/H4RG/HyC/"])

You can probably ignore this, unless you are working with
Nancy Grace Roman Space Telescope test data arcived on the 
NASA ADAPT system. It assumes a specific file hierarchy
and Goddard Detector Characterization Laboratory (DCL) file
naming conventions.

    Parameters: test, String
                  The test description part of a DCL filename;
                  e.g. "_90k_1p1m0p1_noise_"
                SCA, String
                  A Roman SCA serial number; e.g. "20663"
                base, String (optional)
                  Base directory containing the data. The
                  default is for ADAPT. This must be terminated
                  with a "/" character.
    Returns:
      A list of full path filenames. If no files are found, then []
      is returned.
"""
function get_file_list(sca, test; base="/att/projrepo/wfirst/H4RG/HyC/")
    files = []
    dirs = glob("*", base)
    for dir in dirs
        # Check to be sure it is a directory and that we have permission
        # to search
        if isdir(dir) && (gperm(dir) & 0x04 != 0)
            files = glob("*" * test * sca * "_*.fits", dir*"/")
            if length(files) > 0
                return(files)
            end
        end
    end
    return(files)
end