function [mex_name, file_full, desc_das, file_status, exe_name] = choose_mexnc_opendap(file)
% CHOOSE_MEXNC_OPENDAP decides if a file should be read using mexnc or loaddap
%
% [mex_name, file_full, desc_das, file_status, exe_name] = ...
%      choose_mexnc_opendap(file)
%
% INPUT:
% file is the may be a url to an opendap file or it may be the name of a
%  netCDF file, with or without the .cdf or .nc extent. The file may also be
%  in a compressed form, in which the user is offered the choice of having
%  the file uncompressed; this is included for backwards compatibility and
%  its use is deprecated.
%
% OUTPUT:
% mex_name:
%   the name of the mexfile relevant to the given file and it depends on what
%   is available. It may be 'mexnc', 'loaddap', 'loaddods' or 'none'. Of
%   course 'none' means that we can't deal with the file.
% file_full:
%   the same as file but .nc or .cdf may be added to it if that was left off
%   in the original version. If the file is found in the common data set then
%   the path will be prepended.
% desc_das:
%   if we have an opendap dataset then this contains the dds of the file as
%   returned by a call to loaddap or loaddods. Otherwise this will be empty. 
% file_status: a status flag; 
%   = 0 if the netCDF file is in the current directory.
%   = 1 if the file cannot be found anywhere.
%   = 2 if the file is in the directory specified by a call to the
%       m-function pos_cds.
%   = 3 if a compressed version of the file is in the current directory.
%   = 4 if this is a url.
% exe_name:
%   the name of the executable file associated with the mexfile. It contains
%   the full path name so that it can be called from other functions (notably
%   get_dods_dds). It may be writedap, writeval, writedap.exe or
%   writeval.exe. If the mexfile is mexnc then there is no associated
%   executable and so exe_name is 'none'.
%
% NOTES:
% 1) If file is a url then the function will attempt to get the dds first by
%    a call to loaddap. If this works then mex_name is set to 'loaddap'. If
%    this fails it tries loaddods next. If this fails also however then it
%    will be assumed that mexnc is capable of reading the dataset. Although
%    the existence of mexnc will be checked only an elementary test will be
%    done and so the whole thing may fail later on.

% $Id: choose_mexnc_opendap.m Mon, 03 Jul 2006 17:16:40 $
% Copyright J. V. Mansbridge, CSIRO, Fri Oct 28 17:37:44 EST 2005

persistent FUNC_FOR_LOCAL_FILE FUNC_FOR_OPENDAP EXE_NAME_FOR_LOAD

% Get the full name of the file (i.e., including the path name) and find out
% what type of file it is.

[file_full, file_status] = CSIRO_get_more_file_info(file);

% Find out mex_name, desc_das & exe_name as appropriate for the type of
% file. The best way is to use a user supplied routine CSIRO_get_globals but
% if this is not possible then get_mex_name is called. This has a number of
% ways that it can try to find the information.
% This must be the first time choose_mexnc_opendap has been called and so we
% must figure out how the file can be read. Try to find a user supplied
% function that will write out the global variables.

if exist('CSIRO_get_globals')
  if isempty(FUNC_FOR_LOCAL_FILE)
    % This must be the first time choose_mexnc_opendap has been called and so we
    % must figure out how the file can be read. The required information is
    % put into the global variables so that CSIRO_get_globals will not need
    % to be called again.

    [func_for_local_file, func_for_opendap, exe_name_for_load] = ...
	CSIRO_get_globals;
    FUNC_FOR_LOCAL_FILE = func_for_local_file;
    FUNC_FOR_OPENDAP = func_for_opendap;
    EXE_NAME_FOR_LOAD = exe_name_for_load;
  end
  
  % Now use file_status to determine mex_name, desc_das and exe_name. This is
  % particularly messy if loaddap or loaddods is being used.
  
  switch file_status
   case 4
    mex_name = FUNC_FOR_OPENDAP;
    desc_das = [];
    exe_name = 'none';
    switch mex_name
     case 'loaddap'
      try
	desc_das = loaddap('-A +v', file_full);
      catch
	error(['loaddap could not find desc_das for ' file_full])
      end
      exe_name = EXE_NAME_FOR_LOAD;
     case 'loaddods'
      try
	desc_das = loaddods('-A +v', file_full);
      catch
	error(['loaddods could not find desc_das for ' file_full])
      end
      exe_name = EXE_NAME_FOR_LOAD;
     case 'java'
      desc_das = urlread([file '.das']);
      exe_name = 'none';
    end
   otherwise
    mex_name = FUNC_FOR_LOCAL_FILE;
    desc_das = [];
    exe_name = 'none';
  end
else
  opt_choice_method = 2;
  [mex_name, desc_das, exe_name] = get_mex_name(opt_choice_method, ...
						  file_full, file_status);
end

function [mex_name, desc_das, exe_name] = get_mex_name(opt_choice_method, ...
						  file_full, file_status)

% Set opt_choice_method to choose how the method by which the data will be
% retrieved. This can be complicated because not all methods may be
% available. The methods of data retrieval (returned in mex_name) are:
%
% 'java': this uses the matlab java interface to access toolsUI.jar and can
%   read local netcdf files as well as opendap files.
% 'loaddap': uses the latest version of the The Matlab Structs Tool and
%   provides a way to read any OPeNDAP-accessible data into Matlab. This is
%   the second major update to the older 'Command Line Tool' 
% 'loaddods': From the older 'Command Line Tool' - it can read
%   OPeNDAP-accessible data into Matlab.
% 'mexnc': a mex file that can read netcdf files and, in some versions
%   only, can read OPeNDAP-accessible data.
%
%       opt_choice_method - meaning of different settings:
% 1: Always use 'java'
% 2: Use 


% Set mex_name which determines whether we use mexnc, loaddap or loaddods to
% read the dataset.

switch opt_choice_method
 case 1
  mex_name = 'java';
  desc_das = [];
  exe_name = 'none';
 case 2
  switch file_status
   case 4
    % Try to use loaddap to get a description of the opendap dataset. On
    % failure try loaddods. On failure again assume that you can use mexnc.

    mex_name = 'none';
    try
      desc_das = loaddap('-A +v', file_full);
      if exist('desc_das')
	mex_name = 'loaddap';
      end
    catch
      try
	desc_das = loaddods('-A +v', file_full);
	if exist('desc_das')
	  mex_name = 'loaddods';
	end
      catch
	try
	  version_num = mexnc('get_mexnc_info');
	  mex_name = 'mexnc';
	catch
	  mex_name = 'none';
	end
      end
    end
   otherwise
    
    % Do an elementary check for mexnc.
    
    if exist('mexnc')
      mex_name = 'mexnc';
    else
      mex_name = 'none';
    end
    exe_name = 'none';
  end

  % Now find the name and path for the executable called by loaddap or
  % loaddods.

  switch mex_name
   case {'loaddap', 'loaddods'}
    
    % First find the name of the executable.
    
    if isunix
      switch mex_name
       case 'loaddap'
	exe_name_short = 'writedap';
       case 'loaddods'
	exe_name_short = 'writeval';
      end
    else
      switch mex_name
       case 'loaddap'
	exe_name_short = 'writedap.exe';
       case 'loaddods'
	exe_name_short = 'writeval.exe';
      end
    end

    % Now look for the location of the executable. Look first in the same
    % directory as the mex file. If that fails try to find a parallel 'bin'
    % directory as happens in the windows version of the command line tools.

    dir_mex = which(mex_name);
    ff = findstr(dir_mex, mex_name);
    dir_mex = dir_mex(1:(ff(end) - 1));
    exe_name = [dir_mex exe_name_short];
    if exist(exe_name)
      return
    else
      cd(dir_mex)
      cd('..')
      if exist('bin') == 7
	cd('bin')
	if exist(exe_name)
	  if isunix
	    exe_name = [pwd '/' exe_name_short];
	  else
	    exe_name = [pwd '\' exe_name_short];
	  end
	else
	  error(['Couldn''t find ' exe_name_short])
	end
      else
	error(['Couldn''t find ' exe_name_short])
      end
    end
   otherwise
    exe_name = 'none';
  end
end
