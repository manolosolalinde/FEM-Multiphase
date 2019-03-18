function [file_full, file_status] = CSIRO_get_more_file_info(file)

% Check that the file is accessible.  If it is then its full name will
% be stored in the variable full_nam.  The file may have the extent .cdf or
% .nc and be in the current directory or the common data set (whose
% path can be found by a call to pos_cds.m.  If a compressed form
% of the file is in the current directory then the user is prompted to
% uncompress it.  If, after all this, the netcdf file is not accessible
% then the m file is exited with an error message.
  
file_full_list = { '', '.nc', '.cdf'};
ilim = length(file_full_list);
for ii = 1:ilim 
  file_full = [file file_full_list{ii}];
  file_status = check_nc(file_full);

  switch file_status
   case {0, 4}
    break;
   case 1
    if ii == ilim
      error([ file ' could not be found' ])
    end
   case 2
    path_name = pos_cds;
    file_full = [ path_name file_full ];
    break;
   case 3
    err1 = uncmp_nc(file_full);
    if err1 == 0
      break;
    elseif err1 == 1
      disp([ 'exiting because you chose not to uncompress ' file_full ])
      return;
    elseif err1 == 2
      error([ 'exiting because ' file_full ' could not be uncompressed' ])
    end
  end
end
