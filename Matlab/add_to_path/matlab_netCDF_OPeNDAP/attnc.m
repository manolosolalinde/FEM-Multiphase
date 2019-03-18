function [att_val, att_name_list] = attnc(file, var_name, att_name, verbose)
%attnc returns selected attributes of a netcdf file or DODS/OPEnDAP dataset.
%--------------------------------------------------------------------
% DESCRIPTION:
%
% [att_val, att_name_list] = attnc(file, var_name, att_name, verbose)
%       INPUT
% file:  may be the name of a netCDF file with or without the .cdf or .nc
%   extent. file may also be the URL of a DODS/OPEnDAP dataset.
% var_name: a string containing the name of the variable whose
%   attribute is required.  A global attribute may be specified by
%   passing the string 'global'.
% att_name: a string containing the name of the attribute that is
%   required.  If att_name is not specified then it is assumed that the
%   user wants all of the attributes.
% verbose: if == 0 (the default) then no messages about the attributes are
%   displayed. Otherwise some simple messages will be displayed.
%       OUTPUT
% att_val: If att_name is specified then its value is returned in att_val.
%   If att_name is not specified then the values of all of the attributes
%   are returned in the cell att_val.
% att_name_list: If att_name is specified then the same name is returned in
%   att_name_list provided that the attribute is found (otherwise it is
%   empty).  If att_name is not specified then the names of all of the
%   attributes are returned in the cell att_name_list.
%
% Note 1) If only 2 arguments are passed to attnc then all of the
%   attributes and their names (for the specified variable) are returned.
% Note 2) If only 1 argument is passed to attnc then all of the global
%   attributes and their names are returned.
% Note 3) If "file" is a compressed netcdf file then the user will be
%   given an option for automatic uncompression. Of course, the
%   uncompressiion can only occur if the user has appropriate write
%   permission.
%
% This function calls: check_nc.m, loaddap or loaddods, mexnc, pos_cds.m,
%                      uncmp_nc.m
% This function is called by: NONE
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

%     Copyright (C), 1992, J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: attnc.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

% Check the number of arguments and put in the defaults

if ( (nargin < 1) || (nargin > 4) )
   help attnc
   return
end
if nargin == 0
   help attnc
   return  
elseif nargin == 1
  var_name = 'global';
  att_name = [];
  verbose = 0;
elseif nargin == 2
  att_name = [];
  verbose = 0;
elseif nargin == 3
  verbose = 0;
end
if isempty(var_name)
  var_name = 'global';
end

% Do some initialisation.

err_opt = 1;
CSIRO_add_jar_file_maybe;
[mex_name, full_name, desc_das] = choose_mexnc_opendap(file);

switch mex_name
 case 'mexnc'

  % Open the netcdf file.
  
  [cdfid, rcode] = mexnc('ncopen', full_name, 'nowrite');

  % don't print out netcdf warning messages

  mexnc('setopts',0);

  if rcode < 0
    error(['mexnc: ncopen: rcode = ' int2str(rcode)])
  end

  % Find the values of the attributes.

  att_val = [];
  att_name_list = [];
  if strcmp(var_name, 'global')
    %Collect information about the cdf file.

    [num_dims, nvars, ngatts, recdim, rcode] =  mexnc('ncinquire', cdfid);
    if rcode < 0
      error([ 'mexnc: ncinquire: rcode = ' int2str(rcode) ])
    end
    if isempty(att_name)
      if ngatts > 0
	for i = 0:ngatts-1
	  [attnam, rcode] = mexnc('attname', cdfid, 'global', i);
	  [attype, attlen, rcode] = mexnc('ncattinq', cdfid, 'global', attnam);
	  [values, rcode] = mexnc('ncattget', cdfid, 'global', attnam);
	  att_val{i+1} = values;
	  att_name_list{i+1} = attnam;
	end
      else
	if verbose
	  disp('   ---  There are no Global attributes  ---')
	end
      end
    else
      if ngatts > 0
	found_it = 0;
	for i = 0:ngatts-1
	  [attnam, rcode] = mexnc('attname', cdfid, 'global', i);
	  if strcmp(attnam, att_name)
	    found_it = 1;
	    [attype, attlen, rcode] = mexnc('ncattinq', cdfid, 'global', attnam);
	    [values, rcode] = mexnc('ncattget', cdfid, 'global', attnam);
	    att_val = values;
	    att_name_list = attnam;
	    break
	  end
	end
	if found_it == 0
	  warning(['the attribute ' att_name ' was not found'])
	end
      else
	warning(['the attribute ' att_name ' was not found'])
      end
    end
  else
    varid = mexnc('VARID', cdfid, var_name);
    if varid == -1
      error([ 'mexnc: varid: ' var_name ' is not a variable' ])
    end

    [varname, datatype, num_dims, dim, natts, rcode] = ...
	mexnc('VARINQ', cdfid, varid);

    % Now find the attributes and store them. Note that there was previous
    % code to turn a numeric row vector into a column vector. This code was
    % taken out for simplicity and backwards compatibility.
    
    if isempty(att_name)
      if natts > 0
	for i = 0:natts-1
	  [attnam, rcode] = mexnc('attname', cdfid, varid, i);
	  [attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, attnam);
	  [values, rcode] = mexnc('ncattget', cdfid, varid, attnam);
	  att_val{i+1} = values;
	  att_name_list{i+1} = attnam;
	end
      else
	if verbose
	  disp(['   ---   ' varid ' has no attributes  ---'])
	end
      end
    else
      if natts > 0
	found_it = 0;
	for i = 0:natts-1
	  [attnam, rcode] = mexnc('attname', cdfid, varid, i);
	  if strcmp(attnam, att_name)
	    found_it = 1;
	    [attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, attnam);
	    [values, rcode] = mexnc('ncattget', cdfid, varid, attnam);
	    att_val = values;
	    att_name_list = attnam;
	    break
	  end
	end
	if found_it == 0
	  warning(['the attribute ' att_name ' was not found'])
	end
      else
	warning(['the attribute ' att_name ' was not found'])
      end
    end
  end

  % Close the netcdf file.

  [rcode] = mexnc('ncclose', cdfid);
  if rcode < 0
    error(['** ERROR ** ncclose: rcode = ' num2str(rcode)])
  end
  
 case {'loaddap', 'loaddods'}
  
  % Dealing with a dods file

  % Special case. This fiddle is presumably necessary because matlab names
  % cannot start with _ .

  if ~isempty(att_name)
    if strcmp(att_name(1), '_')
      att_name = ['ml_' att_name];
    end
  end

  att_val = [];
  att_name_list = [];

  % Find and store all of the global attributes of the dods file or only one
  % specific value.

  if strcmp(var_name, 'global')
    
    % Find and print out the global attributes of the dods file. We use the
    % information about the DAS that was previously obtained by a call to
    % loaddap or loaddods.

    try
      %desc_das_down_1 = desc_das.('Global_Attributes');
      desc_das_down_1 = getfield(desc_das, 'Global_Attributes');
    catch
      error('desc_das does not have a field Global_Attributes')
    end
    
    fnames_1 = fieldnames(desc_das_down_1);
    found_global = 0;
    for ii = 1:length(fnames_1)
      if length(fnames_1{ii}) > 7
	if strcmpi(fnames_1{ii}(end-6:end), '_global')
	  ii_1 = ii;
	  found_global = 1;
	  break
	end
      end
    end
    
    if found_global
      
      % Find values for all attributes or a specified one.
      
      if isempty(att_name)
	
	% Find values for all attributes.
	
	ngatts = 0;
	%desc_das_down_2 = desc_das_down_1.(fnames_1{ii_1});  
	desc_das_down_2 = getfield(desc_das_down_1, fnames_1{ii_1});  
	if isstruct(desc_das_down_2)
	  list_gatts = fieldnames(desc_das_down_2);
	  for ii = 1:length(list_gatts)
	    if ((strcmp(list_gatts(ii), 'DODS_ML_Real_Name') == 0) && ...
		(strcmp(list_gatts(ii), 'DODS_ML_Type') == 0))
	      ngatts = ngatts + 1;
	      att_name_list{ngatts} = list_gatts{ii};
	      %att_val{ngatts} = desc_das_down_2.(att_name_list{ngatts});
	      att_val{ngatts} = getfield(desc_das_down_2,att_name_list{ngatts});
	    end
	  end
	end
      else
	
	% Find values for a specified attribute.
	
	%desc_das_down_2 = desc_das_down_1.(fnames_1{ii_1});
	desc_das_down_2 = getfield(desc_das_down_1, fnames_1{ii_1});
	if isfield(desc_das_down_2, att_name)
	  %att_val = desc_das_down_2.(att_name);
	  att_val = getfield(desc_das_down_2, att_name);
	  att_name_list = att_name;
	else
	  if verbose
	    warning([ 'the attribute ' att_name ' was not found'])
	  end
	end
      end 
    else
      if verbose
	disp('   ---  There are no Global attributes  ---')
      end
      return
    end
  else
    
    % Check that there are some attributes.
    
    if isfield(desc_das, var_name)
      %fi = desc_das.(var_name);
      fi = getfield(desc_das, var_name);
    else
      error(['There is no variable ' var_name ' in ' file])
    end
    
    % Find values for all attributes or a specified one. Note that a numeric
    % column vector is transformed into a numeric row vector. This is to be
    % compatible with the output from the equivalent action when a netcdf
    % file is read directly from the local disk.
    
    if isempty(att_name)
      att_list = fieldnames(fi);
      natts = 1;
      while 1
	if isempty(findstr(att_list{natts}, 'DODS_ML_'))
	  att_name_list{natts} = att_list{natts};
	  %values = fi.(att_list{natts});
	  values = getfield(fi, att_list{natts});
	  
	  % A numeric column vector is transformed into a numeric row vector.
	  
	  if isnumeric(values) && (ndims(values) == 2)
	    [m_temp, n_temp] = size(values);
	    if n_temp == 1
	      values = values';
	    end
	  end
	  att_val{natts} = values;
	  natts = natts + 1;
	else
	  natts = natts - 1;
	  break
	end
      end
      if natts == 0
	if verbose
	  disp(['   ---  ' var_name ' has no attributes  ---'])
	end
      end  
    else
      if isfield(fi, att_name)
	%att_val = fi.(att_name);
	att_val = getfield(fi, att_name);
	
	% A numeric column vector is transformed into a numeric row vector.
	
	if isnumeric(att_val) && (ndims(att_val) == 2)
	  [m_temp, n_temp] = size(att_val);
	  if n_temp == 1
	    att_val = att_val';
	  end
	end
	att_name_list = att_name;
      else
	if verbose
	  warning([var_name ' does not have an attribute ' att_name])
	  return
	end
      end
    end
  end
  
  % Call clean_up_string to replace any control characters with a # to avoid
  % messing up the display - null characters make a major mess. clean_up_string
  % also strip off extraneous quote marks. We also check att_name_list to see if
  % any attributes start with 'ml__'. If they do then we cut off the 'ml_' under
  % the assumption that it has been added to deal with the problem of matlab
  % variables starting with an underscore.

  if iscell(att_val)
    for ii = 1:length(att_val)
      att_val{ii} = clean_up_string(att_val{ii});
      if length(att_name_list{ii}) > 4
	if strcmp(att_name_list{ii}(1:4), 'ml__')
	  att_name_list{ii} = att_name_list{ii}(4:end);
	end
      end
    end
  else
    att_val = clean_up_string(att_val);
    if length(att_name_list) > 4
      if strcmp(att_name_list(1:4), 'ml__')
	att_name_list = att_name_list(4:end);
      end
    end
  end
 case 'java'
  % Open dataset so that all attributes will be accessible.

  try
    enhance = false;
    ncdJ = ucar.nc2.dataset.NetcdfDataset.openDataset(file, enhance, []);
  catch
    ss = lasterror;
    mess_str = ss.message;
    rcode = -1000000;
    values = error_handle('java', mess_str, rcode, err_opt);
  end
  
  % Initialise some flags.
  
  if strfind(var_name, 'global')
    var_is_global = 1;
  else
    var_is_global = 0;
  end
  if isempty(att_name)
    get_all_atts = 1;
  else
    get_all_atts = 0;
  end
  
  % Get the variable object if this is necessary.
  
  if ~var_is_global
    try
      varJ = ncdJ.findVariable(var_name);
    catch
      ss = lasterror;
      mess_str = ss.message;
      rcode = -1000000;
      values = error_handle(ncdJ, mess_str, rcode, err_opt);
    end
  end
  
  % Actually get the attribute values.
  
  if var_is_global
    if get_all_atts
      att_list = ncdJ.getGlobalAttributes();
      num_atts = att_list.size();
      att_name_list = cell(num_atts, 1);
      att_val = cell(1, num_atts);
      for ii = 1:num_atts
	att = att_list.get(ii - 1);
	[att_name_list{ii}, att_val{ii}] = get_att_val(att, err_opt);
      end
    else
      try
	att = ncdJ.findGlobalAttribute(att_name);
      catch
	ss = lasterror;
	mess_str = ss.message;
	rcode = -1000000;
	values = error_handle(ncdJ, mess_str, rcode, err_opt);
      end
      [att_name_list, att_val] = get_att_val(att, err_opt);      
    end
  else
    if get_all_atts
      att_list = varJ.getAttributes();
      num_atts = att_list.size();
      att_name_list = cell(num_atts, 1);
      att_val = cell(1, num_atts);
      for ii = 1:num_atts
	att = att_list.get(ii - 1);
	[att_name_list{ii}, att_val{ii}] = get_att_val(att, err_opt);
      end
    else
      try
	att = varJ.findAttribute(att_name);
      catch
	ss = lasterror;
	mess_str = ss.message;
	rcode = -1000000;
	values = error_handle(varJ, mess_str, rcode, err_opt);
      end
      [att_name_list, att_val] = get_att_val(att, err_opt);
    end
  end
  
 case 'none'
  error(['Couldn''t find a suitable mex-file for reading ' file])
end

function new_val = clean_up_string(val)
% If val is a string then we replace any control characters with a # to
% avoid messing up the display - null characters make a major mess. We
% also strip off extraneous quote marks.
  if ischar(val)
    s = abs(val);
    if (s(1) == s(end)) && ((s(1) == 34) || (s(1) == 39))
      s = s(2:(end - 1));
    end
     s(find(s < 32)) = 35;
    new_val = char(s);
  else
    new_val = val;
  end

function [att_name_list, att_val] = get_att_val(att, err_opt)
  try
    att_name_list =  char(att.getName());
  catch
    ss = lasterror;
    mess_str = ss.message;
    rcode = -1000000;
    values = error_handle(att, mess_str, rcode, err_opt);
  end
  if att.isString()
    att_val = char(att.getStringValue());
  else
    num = double(att.getLength);
    vec = zeros(1, num);
    for mm = 1:num
      vec(mm) = double(att.getNumericValue(mm - 1));
    end
    att_val = vec;
  end
function values = error_handle(cdfid, mess_str, rcode, err_opt)

% error_handle is called after a mexnc call has failed. It ensures
% that an open netcdf file is closed. The value of err_opt determines what
% else is done.
%    err_opt == 1 prints an error message and then aborts
%            == 2 prints a warning message and then returns an empty
%                 array. This is the default.
%            == 3 returns an empty array. This is a very dangerous option and
%                 should only be used with caution. It might be used when
%                 getnc_s is called in a loop and you want to do your own
%                 error handling without being bothered by warning messages.

% Decide what part of the code made the call.
  
  if isempty(cdfid)
    called_by = 'loadd';
  else
    if isnumeric(cdfid)
      called_by = 'mexnc';
    else
      called_by = 'java';
    end
  end

  % Close an open netcdf of java file.
  
  switch called_by
   case 'mexnc'
    if cdfid >= 0
      [rcode_sub] = mexnc('ncclose', cdfid);
    end
   case 'java'
    if isjava(cdfid)
      ncdJ.close();
    end
  end
  
  % Handle the errors according to the value of err_opt. If rcode is empty
  % then this is probably because loaddap or loaddods was called.
  
  values = [];
  switch err_opt
   case 1
    if isempty(rcode)
      str = ['ERROR: ' mess_str];
    else
      str = ['ERROR: ' mess_str ' : rcode = ' num2str(rcode)];
    end
    error(str)
   case 2
    if isempty(rcode)
      str = ['WARNING: ' mess_str];
    else
      str = ['WARNING: ' mess_str ' : rcode = ' num2str(rcode)];
    end
    disp(str)
   case 3
    return
   otherwise
    error(['error_handle was called with err_opt = ' num2str(err_opt)])
  end

  