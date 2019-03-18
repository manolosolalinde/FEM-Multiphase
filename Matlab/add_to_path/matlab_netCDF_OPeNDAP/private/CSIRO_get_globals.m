function [func_for_local_file, func_for_opendap, exe_name_for_load] = ...
    CSIRO_get_globals
  
  % CSIRO_get_globals returns 3 strings which specify the drivers to be used
  % for reading local and opendap files.
  %
  %    func_for_local_file:
  % The name of the driver for reading local netcdf files. It may be
  % 1) 'mexnc': the mexnc file
  % 2) 'java': the toolsUi.jar/java system
  %
  %    func_for_opendap:
  % The name of the driver for reading the opendap files. It may be
  % 1) 'java': the toolsUi.jar/java system
  % 2) 'loaddap': the latest version of the Matlab Command Line Tool
  % 3) 'loaddods': the old version of the Matlab Command Line Tool
  % 4) 'mexnc': the mexnc file can read oendap files for some machines only
  %             and is rarely used for this.
  %
  %    exe_name_for_load:
  % This is the name of the  executable required if the Matlab Command Line
  % Tool is used to read opendap files. It may be
  % 1) 'none': not using the Matlab Command Line Tool
  % 2) 'writedap': if using 'loaddap' on a linux box
  % 3) 'writedap.exe': if using 'loaddap' on a windows box
  % 4) 'writeval': if using 'loaddods' on a linux box
  % 5) 'writeval.exe': if using 'loaddods' on a windows box
  
  % EDIT THE LINE BELOW TO CHOOSE AN OPTION.
  opt = 2;
  
  % Set the globals that specify how netcdf and opendap files are to be read.
  
  switch opt
   case 1
    func_for_local_file = 'java';
    func_for_opendap = 'java';
    exe_name_for_load = 'none';
   case 2
    func_for_local_file = 'mexnc';
    func_for_opendap = 'java';
    exe_name_for_load = 'none';
   case 3
    func_for_local_file = 'mexnc';
    func_for_opendap = 'loaddap';
    
    % Find the directory containing the files used for opendap access. It is
    % assumed that they are in the same directory as getnc.
    directory = which('getnc');
    path_name = directory(1:(end-7));
    if isunix
      exe_name_for_load = [path_name 'writedap'];
    else
      exe_name_for_load = [path_name 'writedap.exe'];
    end
  end
  