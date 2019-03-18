function values = getnc(varargin);
% function values = getnc(file, varid, bl_corner, tr_corner, stride, order, ...
%      change_miss, new_miss, squeeze_it, rescale_opts, err_opt)
%
%  GETNC retrieves data from a NetCDF file or a DODS/OPEnDAP dataset.
%
%  function values = getnc(file, varid, bl_corner, tr_corner, stride, order, ...
%        change_miss, new_miss, squeeze_it, rescale_opts, err_opt)
%
% DESCRIPTION:
%  getnc retrieves data either from a local NetCDF file or a DODS/OPEnDAP
%  dataset. The way getnc behaves depends on how many of the input arguments are
%  passed to it. If no arguments are passed then it returns this help
%  message. If one argument (the name of a netcdf file) is passed then the user
%  is asked questions to determine information necessary for the data
%  retrieval. (DODS/OPEnDAP data cannot be retrieved this way.) Other usage
%  (described below) is for the non-interactive behaviour. If more than one
%  argument is passed then getnc returns the data without needing to ask any
%  questions. The input arguments are listed below.
%
% USAGE:
% getnc retrieves data in two ways. It can be used used interactively to
% retrieve data from a netCDF file by simply typing:
%
% >> values = getnc(file);
%
% getnc is more commonly used as a function call - it can then retrieve data
% from both netCDF and OPeNDAP files. Because many options are available
% getnc can take up to 11 input arguments (although most have default
% values). To make things easier for the user there are various ways of
% specifying these arguments.  Specifying up to 11 arguments to getnc can be
% complicated and confusing. To make the process easier getnc will accept a
% variety of types of input. These are given as follows:
%
%    1) Specify all 11 arguments. Thus we could make a call like:
%
% >> values = getnc(file, varid, bl_corner, tr_corner, stride, order, ...
%             change_miss, new_miss, squeeze_it, rescale_opts, err_opt);
%
%    2) Use default arguments. Only the first 2 arguments are strictly
%    necessary as the other arguments all have defaults. The following call
%    would retrieve the entire contents of the named variable:
%
% >> values = getnc(file, varid);
%
% If you want non-default behaviour for one or more of the later arguments
% then you can do something like:
%
% >> values = getnc(file, varid, -1, -1, -1, -1, change_miss, new_miss);
%
% In this case there are 4 arguments specified and 7 with default values used.
%
%    3) Use a structure as an argument. From version 3.3 onwards it is
%    possible to pass a structure to getnc. This is illustrated below: 
% 
% >> x.file = 'fred.nc';
% >> x.varid = 'foo';
% >> x.change_miss = 1;
% >> values = getnc(x);
% 
% This specifies 3 arguments and causes defaults to be used for the other 8.
% Note that it is possible to mix the usual arguments with the passing of a
% structure - it is only necessary that the structure be the last argument
% passed. We could achieve the same effect as above by doing:
% 
% >> x.change_miss = 1;
% >> values = getnc('fred.nc', 'foo', x);
%
% INPUT ARGUMENTS:
%  1. file: This is a string containing the name of the netCDF file or the
%   URL to the OpenDAP dataset. It does not have a default. If describing a
%   netCDF file it is permissible to drop the ".nc" prefix but this is not
%   recommended.
%
%  2. varid:  This may be a string or an integer. If it is a string then it
%   should be the name of the variable in the netCDF file or OPEnDAP
%   dataset. The use of an integer is a deprecated way of accessing netCDF
%   file data; if used the integer is taken from the order of the variables
%   returned by a call to ddsnc or inqnc (starting from 1).
%
%  3. bl_corner: This is a vector of length n specifying the hyperslab
%   corner with the lowest index values (the bottom left-hand corner in a
%   2-space).  The corners refer to the dimensions in the same order that
%   these dimensions are listed in the inqnc description of the variable. For
%   a netCDF file this is the same order that they are returned in a call to
%   "ncdump". With an OPEnDAP dataset it is the same order as in the
%   DDS. Note also that the indexing starts with 1 - as in matlab and
%   fortran, NOT 0 as in C. A negative element means that all values in that
%   direction will be returned.  If a negative scalar (or an empty array) is
%   used this means that all of the elements in the array will be
%   returned. This is the default, i.e., all of the elements of varid will be
%   returned.
%
%  4. tr_corner: This is a vector of length n specifying the hyperslab
%   corner with the highest index values (the top right-hand corner in a
%   2-space). A negative element means that the returned hyperslab should run
%   to the highest possible index (this is the default). Note, however, that
%   the value of an element in the end_point vector will be ignored if the
%   corresponding element in the corner vector is negative.
%
%  5. stride: This is a vector of length n specifying the interval between
%   accessed values of the hyperslab (sub-sampling) in each of the n
%   dimensions.  A value of 1 accesses adjacent values in the given
%   dimension; a value of 2 accesses every other value; and so on. If no
%   sub-sampling is required in any direction then it is allowable to just
%   pass the scalar 1 (or -1 to be consistent with the corner and end_point
%   notation). Note, however, that the value of an element in the stride
%   vector will be ignored if the corresponding element in the corner vector
%   is negative.
%
%  6. order: 
%     * order == -1 then the n dimensions of the array will be returned in
%     the same order as described by a call to inqnc(file) or "ncdump". It
%     therefore corresponds to the order in which the indices are specified
%     in corner, end_point and stride. This is the default.
%     * order == -2 will reverse the above order. Because matlab's array
%     storage is row-dominant this is actually a little quicker but the
%     difference is rarely significant.
%
%  7. change_miss: Missing data are indicated by the attributes _FillValue,
%   missing_value, valid_range, valid_min and valid_max. The action to be
%   taken with these data are determined by change_miss.
%     * change_miss == 1 causes missing values to be returned unchanged.
%     * change_miss == 2 causes missing values to be changed to NaN (the
%     default).
%     * change_miss == 3 causes missing values to be changed to new_miss
%     (after rescaling if that is necessary).
%     * change_miss < 0 produces the default (missing values to be changed to
%     NaN).
%
%  8. new_miss: This is the value given to missing data if change_miss == 3.
%
%  9. squeeze_it: This specifies whether the matlab function "squeeze"
%   should be applied to the returned array. This will eliminate any
%   singleton array dimensions and possibly cause the returned array to have
%   less dimensions than the full array.
%     * squeeze_it ~= 0 causes the squeeze function to be applied.  This is
%     the default. Note also that a 1-d array is returned as a column
%     vector.
%     * squeeze_it == 0 so that squeeze will not be applied.
%
% 10. rescale_opts: This is a 2 element vector specifying whether or not
%  rescaling is carried out on retrieved variables and certain
%  attributes. The relevant attributes are _FillValue', 'missing_value',
%  'valid_range', 'valid_min' and 'valid_max'; they are used to find missing
%  values of the relevant variable. The option was put in to deal with files
%  that do not follow the netCDF conventions (usually because the creator of
%  the file has misunderstood the convention). For further discussion of the
%  problem see here. Only use this option if you are sure that you know what
%  you are doing.
%     * rescale_opts(1) == 1 causes a variable read in by getnc.m to be
%     rescaled by 'scale_factor' and  'add_offset' if these are attributes of
%     the variable; this is the default.
%     * rescale_opts(1) == 0 disables rescaling of the retrieved variable.
%     * rescale_opts(2) == 1 causes the attributes '_FillValue', etc to be
%     rescaled by 'scale_factor' and 'add_offset'; this is the default.
%     * rescale_opts(2) == 0 disables the rescaling of the attributes
%     '_FillValue', etc.
%
% 11. err_opt: This is an integer that controls the error handling in a call
%  to getnc.
%     * err_opt == 0 on error this prints an error message and aborts.
%     * err_opt == 1 prints a warning message and then returns an empty
%     array. This is the default.
%     * err_opt == 2 returns an empty array. This is a dangerous option and
%     should only be used with caution. It might be used when getnc is called
%     in a loop and you want to do your own error handling without being
%     bothered by warning messages.
%
% OUTPUT:
%  values is a scalar, vector or array of values that is read in
%     from the NetCDF file or DODS/OPEnDAP dataset
%
% NOTES:
%   1) In order for getnc to work non-interactively it is only strictly
% necessary to pass the first 2 input arguments to getnc - sensible
% defaults are available for the rest.
% The defaults are:
% bl_corner, tr_corner == [-1 ... -1], => all elements retrieved
% stride == 1, => all elements retrieved
% order == -1
% change_miss == 2, => missing values replaced by NaNs
% new_miss == 0;
% squeeze_it == 1; => singleton dimensions will be removed
% rescale_opts == [1 1]; => the obvious rescaling
% error_opt == 1 prints a warning message and then returns an empty array.
%
%   2) It is not acceptable to pass only 3 input arguments since there is
% no default in the case of the corner points being specified but the
% end points not.
%
%   3) By default the order of the dimensions of a returned array will be the
% same as they appear in the relevant call to 'inqnc' (from matlab) or
% 'ncdump -h' (from the command line).  (This is the opposite to what
% happened in an earlier version of getnc.)  For a netcdf file this actually
% involves getnc re-arranging the returned array because the netCDF utilities
% follow the C convention for data storage and matlab follows the fortran
% convention. For a DODS/OPEnDAP dataset it is even weirder.
%
%   4) If the values are returned in a one-dimensional array then this will
% always be a column vector.
%
%   5) A strange 'feature' of matlab 5 is that it will not tolerate a singleton
% dimension as the final dimension of a multidimensional array.  Thus, if
% you chose to have only one element in the final dimension this dimension
% will be 'squeezed' whether you want it to be or not - this seems to be
% unavoidable.
%
%   6) Some earlier versions of this function the argument "order" to be an
% array. This option has been removed because it was so confusing - the
% matlab function "permute" can be used to do the same thing.
%
% EXAMPLES:
% 1) Get all the elements of the variable, note the order of the dimensions:
% >> airtemp = getnc('oberhuber.nc', 'airtemp');
% >> size(airtemp)
% ans =
%     12    90   180
%
% 2) Get a subsample of the variable, note the stride:
% >> airtemp = getnc('oberhuber.nc', 'airtemp', [-1 1 3], [-1 46 6], [1 5 1]);
% >> size(airtemp)
% ans =
%     12    10     4
%
% 3) Get all the elements of the variable, but with missing values
%    replaced with 1000.  Note that the bl_corner, tr_corner, stride and
%    order vectors have been replaced by -1 to choose the defaults:
% >> airtemp = getnc('oberhuber.nc', 'airtemp', -1, -1, -1, -1, 3, 1000); 
% >> size(airtemp)
% ans =
%     12    90   180
%
% 4) Get a subsample of the variable, a singleton dimension is squeezed:
% >> airtemp = getnc('oberhuber.nc', 'airtemp', [-1 7 -1], [-1 7 -1]);   
% >> size(airtemp)                                                         
% ans =
%     12   180
% 
% 5) Get a subsample of the variable, a singleton dimension is not squeezed:
% >> airtemp = getnc('oberhuber.nc','airtemp',[-1 7 -1],[-1 7 -1],-1,-1,-1,-1,0);
% >> size(airtemp)                                                            
% ans =
%     12     1   180
%
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

% This function calls: check_nc.m, check_st.m, fill_att.m, fill_var.m,
%                      getnc_s.m, inqnc.m, menu_old.m, mexnc, pos_cds.m,
%                      return_v.m, uncmp_nc.m, y_rescal.m

%     Copyright (C), J.V. Mansbridge, 1992
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: getnc.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

% In November 1998 some code was added to deal better with byte type data in a
% netcdf file. Note that any values greater than 127 will have 256 subtracted
% from them. This is because on some machines (an SGI running irix6.2 is an
% example) values are returned in the range 0 to 255. Note that in the fix the
% values less than 128 are unaltered and so we do not have to know whether the
% particular problem has occurred or not; for machines where there is no problem
% no values will be altered. This is applied to byte type attributes (like
% _FillValue) as well as the variable values.

% Check the number of arguments.  If there are no arguments then return
% the help message.  If there is more than one argument then call
% getnc_s which reads the netcdf file in a non-interactive way.
% If there is only one argument then drop through and find the values
% interactively.

num_args = length(varargin);
if num_args == 0
  % No input argument
  help getnc
  disp('You must pass an input argument to getnc')
  return
elseif num_args == 1
  if isstruct(varargin{1})
    % Send on to getnc_s
    values = getnc_s(varargin);
    return
  else
    % Carry on with the interactive version of getnc
    file = varargin{1};
  end
else
  % Send on to getnc_s
  values = getnc_s(varargin);
  return
end

% Find the full path name of the file.

[cdf, file_status] = CSIRO_get_more_file_info(file);

% Find out whether values should be automatically rescaled or not.

[rescale_var, rescale_att] = y_rescal;

% Get the dds for the file.

desc = ddsnc(cdf);

% varstring = fill_var(cdfid, nvars);

% Prompt the user for the name of the variable containing the hyperslab.

nvars = length(desc.variable);
list_variable_name = cell(nvars, 1);
for ii = 1:nvars
  list_variable_name{ii} = desc.variable(ii).name;
end
i_var = menu('Choose a variable', list_variable_name);

% Get attribute and dimension information for the chosen variable.

var_name = list_variable_name{i_var};
[att_val, att_name_list] = attnc(cdf, var_name);
dim_idents = desc.variable(i_var).dim_idents;
num_dims_of_var = length(dim_idents);
list_dim_name = cell(num_dims_of_var, 1);
list_dim_length = zeros(num_dims_of_var, 1);
is_dim_a_variable = zeros(num_dims_of_var, 1);
for ii = 1:num_dims_of_var
  list_dim_name{ii} = desc.dimension(dim_idents(ii)).name;
  list_dim_length(ii) = desc.dimension(dim_idents(ii)).length;
  for jj = 1:nvars
    if strcmp(list_dim_name{ii}, list_variable_name{jj})
      is_dim_a_variable(ii) = 1;
      break
    end
  end
end

% Call simple input dialog so user can select the features they need. Note
% that this is in a loop in case the user messes up the answer.

dlg_title = [var_name ': input details'];
num_lines = 1;
for ii = 1:num_dims_of_var
  prompt{ii} = [list_dim_name{ii} ...
		': [starting index, finishing index, stride]'];
  defAns{ii} = ['[1, ' num2str(list_dim_length(ii)) ', 1]'];
end
ndvp = num_dims_of_var + 1;
prompt{ndvp} = 'replace missing values with NaN? (y or n)';
defAns{ndvp} = 'y';

Resize = 'on';
repeat_loop = 1;
while repeat_loop
  
  % Get user respnse
  
  answer = inputdlg(prompt, dlg_title, num_lines, defAns, Resize);

  % Now parse the reply to see if it is o.k.
  
  repeat_loop = 0;
  for ii = 1:num_dims_of_var
    try
      xx = str2num(answer{ii});
      bl_corner(ii) = xx(1);
      tr_corner(ii) = xx(2);
      stride(ii) = xx(3);
    catch
      repeat_loop = repeat_loop + 1;
      err_msg{repeat_loop} = ['Dimension ' list_dim_name{ii} ...
		    ': specify [starting index, finishing index, stride]'];
    end
    switch lower(answer{ndvp})
     case 'y'
      change_miss = 2;
     case 'n'
      change_miss = 1;
     otherwise
      repeat_loop = repeat_loop + 1;
      err_msg{repeat_loop} = 'Replacement of missing values specified wrongly';
    end
  end
  
  % Deal with any errors that the user made in the dialog box.
  
  if repeat_loop == 1
    errordlg(err_msg{1}, 'Error specifying array')
  elseif repeat_loop > 1
    clear err_list
    err_list{1} = ['There were ' num2str(repeat_loop) ' errors in the input' ...
		   ' dialog box'];
    for ii = 1:repeat_loop
      err_list{ii + 1} = err_msg{ii};
    end
    errordlg(err_list, 'Error specifying array')
  end
end
values = getnc(cdf, var_name, bl_corner, tr_corner, stride, -1, change_miss, ...
		    0, 1);
