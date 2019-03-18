function inqnc(file)
% INQNC interactively returns info about a netcdf file or DODS/OPEnDAP dataset
%--------------------------------------------------------------------
%     Copyright (C) J. V. Mansbridge, CSIRO, january 24 1992
%     Revision $Revision: 1.20 $
%
%  function inqnc(file)
%
% DESCRIPTION:
%  inqnc('file') is an interactive function that returns information
%  about a netcdf file or DODS dataset. 
% INPUT:
%  file may be the name of a netCDF file with or it may be the URL of a
%  DODS/OPEnDAP dataset. 
%
% OUTPUT:
%  information is written to the user's terminal.
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

%     Copyright (C), J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: inqnc.m Mon, 03 Jul 2006 17:16:40 $
% 
% Note that the netcdf functions are accessed by reference to the mex
% function mexnc. The DODS/OPEnDAP use the Matlab Structs tool.
%--------------------------------------------------------------------

% This function calls: check_nc.m, get_dods_dds.m, loaddap or loaddods,
%                      mexnc, pos_cds.m, uncmp_nc.m, 
% This function is called by:

% Check the number of arguments.

if nargin < 1
  help inqnc
  return
end

% Do some initialisation.

blank = abs(' ');
CSIRO_add_jar_file_maybe;
[mex_name, full_name, desc_das, file_status, exe_name] = ...
    choose_mexnc_opendap(file);

% Find the full path name of the file.

[cdf, file_status] = CSIRO_get_more_file_info(file);

% Get the dds for the file.

try
  desc = ddsnc(cdf);
catch
  error(['Could not open ' file])
end

% Print out iinformation about global attributes.

[att_val_global, att_name_list_global] = attnc(cdf, 'global');
num_global_atts = length(att_val_global);
if num_global_atts > 0
  disp('                ---  Global attributes  ---')
  for ii = 1:num_global_atts
    att = att_val_global{ii};
    if isnumeric(att)
      att = num2str(att);
    end
    disp([att_name_list_global{ii} ': ' att])
  end
else
  disp('   ---  There are no Global attributes  ---')
end

% Print out information about dimensions.

num_dims = length(desc.dimension);
disp(' ')
if num_dims > 0
  str = ['The ' num2str(num_dims) ' dimensions are '];
  for ii = 1:num_dims
    str = [str ' '  num2str(ii) ') ' desc.dimension(ii).name ' = ' ...
	   num2str(desc.dimension(ii).length)];
  end
  disp(str)
else
  disp('   ---  There are no dimensions  ---')
end

% Initialise for variable calculations.

nvars = length(desc.variable);
list_variable_name = cell(nvars, 1);
for ii = 1:nvars
  list_variable_name{ii} = desc.variable(ii).name;
end

% Let user choose the info they get.

i_choice = menu('get information about', 'all variables', 'some variables');

% Loop to prompt the user for the name of the variable containing the
% hyperslab.

if nvars > 0
  repeat_loop = 1;
else
  repeat_loop = 0;
end
while repeat_loop
  switch i_choice
   case 1
    i_var = repeat_loop;
   case 2
    i_var = menu('Choose a variable', list_variable_name);
  end

  % Get attribute and dimension information for the chosen variable.

  var_name = list_variable_name{i_var};
  [att_val, att_name_list] = attnc(cdf, var_name);
  num_atts = length(att_val);
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
  
  % Print out the info about the variable.
  
  s = [ '   ---  Information about ' var_name '(' ];	
  index_list = desc.variable(i_var).dim_idents;
  nvdims = length(index_list);
  for j = 1:nvdims
    s = [ s  desc.dimension(index_list(j)).name ' ' ];
  end
  s = [ s ')  ---' ];
  disp(' ')
  disp(s)
  if num_atts > 0
    for ii = 1:num_atts
      att = att_val{ii};
      if isnumeric(att)
	att = num2str(att);
      end
      disp([att_name_list{ii} ': ' att])
    end
  else
    disp([var_name ' has no attributes'])
  end
  
  % Decide about printing info about the next variable.

  switch i_choice
   case 1
    if repeat_loop == nvars
      repeat_loop = 0;
    else
      repeat_loop = repeat_loop + 1;
    end
   case 2
    ii = menu('Print information about another variable?', 'yes', 'no');
    switch ii
     case 1
      repeat_loop = 1;
     case 2
      repeat_loop = 0;
    end
  end
end

