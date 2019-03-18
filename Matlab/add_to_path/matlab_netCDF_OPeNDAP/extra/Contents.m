% MATLAB interface to netCDF files and DODS/OPEnDAP datasets
%
% The following m-files were written and developed by Jim Mansbridge (C)
% CSIRO 1992-2005.  They use the mex-files maintained by John Evans (see
% http://mexcdf.sourceforge.net/). For the DODS/OPEnDAP operations they use
% the Matlab Structs Tool (see http://www.opendap.org/download/ml-structs).
%
%==========================================================================
% The following commands are available (use help for further
% information) and will be supported in the future.
%==========================================================================
%
% attnc         - imports attributes of a netCDF file or DODS/OPEnDAP dataset
% ddsnc         -  gets information about a netcdf file or DODS/OPEnDAP dataset
% getnc         - imports variables from a netCDF file or DODS/OPEnDAP dataset
% inqnc         - interactive inquiry of netCDF file or DODS/OPEnDAP dataset
% timenc        - finds the time vector and the base date for a netcdf file or
%                 DODS/OPEnDAP dataset that meets COARDS standards 
% whatnc        - lists netCDF files in current directory
%
%==========================================================================
% The following commands are often used by the routines above (use help for
% further information) and will be supported in the future.
%==========================================================================
% check_nc          - checks whether a netcdf or opendap file is accessible
% check_st          - checks if a string is stored in a cell or character array
% choose_mexnc_opendap - decides if a file should be read using mexnc or loaddap
% fill_att          - returns the names of the attributes of a variable
% fill_var          - fills a row of an array with the name of a variable
% get_calendar_date - converts Julian day numbers to calendar dates
% get_day_of_week   - gets the day of the week given the Julian day number
% get_dods_dds      - finds the Dataset Descriptor Structure (DDS) for a
%                       DODS/OPEnDAP data set
% get_julian_day    - converts calendar dates to Julian day numbers
% getnc_s           - command line version of getnc
% get_serial_time   - find the absolute time value given a relative time
% gregorian         - converts Julian day numbers to Gregorian calendar date
% hms2h             - converts hours, minutes, and seconds to hours
% julian            - converts Gregorian calendar dates to Julian day numbers
% loaddap           - direct calls to get DODS/OPEnDAP data and metadata
% menu_old          - generate a menu of choices for user input
% mexnc             - direct calls to netcdf functions
% parsetnc          - parses the COARDS string that specifies time units
% pos_cds           - returns the path to the common data set directory
% return_v          - prints out a string to prompt the user for keyboard input
% s2hms             - converts seconds to integer hour,minute,seconds
% tmexnc            - test of mexnc Mex-file interface to NetCDF
% uncmp_nc          - offers to uncompress a netCDF file
% y_rescal          - returns the scalars rescale_var and rescale_att
%
%==========================================================================

%     $Id: Contents.m Mon, 03 Jul 2006 17:16:40 $
