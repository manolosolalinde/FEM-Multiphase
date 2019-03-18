function values = getnc_s(cell_args)
% getnc_s is called by getnc in the script or non-interactive case. The
% arguments passed to getnc have been bundled into the cell named cell_args
% which is passed to getnc_s. For documentation do a help on getnc.
  
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------
%
%     Copyright (C), J.V. Mansbridge, 1992
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: getnc_s.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

% The function getnc_s calls the following external functions:
%      check_nc.m, check_st.m, fill_att.m, get_dods_dds.m, loaddap or
%      loaddods, mexnc.m, pos_cds.m, uncmp_nc.m, y_rescal.m.
% getnc_s also calls the internal functions :
%      remove_quotes, error_handle, parse_args, substitute_new.
%
% This function is called by: getnc.m

% Written by Jim Mansbridge december 10 1991

% In November 1998 some code was added to deal better with byte type data. Note
% that any values greater than 127 will have 256 subtracted from them. This is
% because on some machines (an SGI running irix6.2 is an example) values are
% returned in the range 0 to 255. Note that in the fix the values less than 128
% are unaltered and so we do not have to know whether the particular problem has
% occurred or not; for machines where there is no problem no values will be
% altered. This is applied to byte type attributes (like _FillValue) as well as
% the variable values.

% Check that there are the correct number of arguments. If the 6th, 7th
% or 8th arguments are missing then they are set here.  If the 3rd, 4th
% or 5th arguments are missing then their defaults will be set later
% when we find out the dimensions of the variable.  It produces an error
% if the bl_corner is set but the tr_corner is not.

[file, varid, bl_corner, tr_corner, stride, order, change_miss, new_miss, ...
 squeeze_it, rescale_opts, err_opt] = parse_args(cell_args);

% Set whether values will be automatically rescaled or not.

if isnumeric(rescale_opts)
  if ndims(rescale_opts) ~= 2
    error('rescale_opts must be [], -1 or a 2 element vector');
  end
  len_rescale_opts = length(rescale_opts);
  if len_rescale_opts == 2
    [rescale_var, rescale_att] = y_rescal(rescale_opts);
  elseif len_rescale_opts < 2
    [rescale_var, rescale_att] = y_rescal;
  else
    error('rescale_opts must be [], -1 or a 2 element vector')
  end
else 
  error('rescale_opts must be numeric');
end

% Set some constants.

blank = abs(' ');

% Do some initialisation.

CSIRO_add_jar_file_maybe;
try
  [mex_name, full_name, desc_das, file_status, exe_name] = ...
      choose_mexnc_opendap(file);
catch
  cdfid = -1000;
  mess_str = lasterr;
  rcode = -1001;
  values = error_handle(cdfid, mess_str, rcode, err_opt);
  return
end
  

switch mex_name
 case 'mexnc'
  
  % I make mexnc calls to find the integers that specify the attribute
  % types

  nc_byte = mexnc('parameter', 'nc_byte'); %1
  nc_char = mexnc('parameter', 'nc_char'); %2
  nc_short = mexnc('parameter', 'nc_short'); %3
  nc_long = mexnc('parameter', 'nc_long'); %4
  nc_float = mexnc('parameter', 'nc_float'); %5
  nc_double = mexnc('parameter', 'nc_double'); %6

  % Set the value of imap.  Note that this is used simply as a
  % placeholder in calls to vargetg - its value is never used.

  imap = 0;

  % Open the netcdf file.

  [cdfid, rcode] = mexnc('ncopen', full_name, 'NC_NOWRITE');
  if rcode ~= 0
    mess_str = ['ncopen ' full_name];
    values = error_handle(cdfid, mess_str, rcode, err_opt);
    return
  end

  % Suppress all error messages from netCDF 

  [rcode] = mexnc('setopts', 0);

  % Collect information about the cdf file.

  [ndimens, nvars, ngatts, recdim, rcode] =  mexnc('ncinquire', cdfid);
  if rcode < 0
    mess_str = ['ncinquire: ' full_name];
    values = error_handle(cdfid, mess_str, rcode, err_opt);
    return
  end

  % Determine the netcdf id number for the required variable.  If varid is a
  % string then an appropriate call to mexnc is used to convert it to the
  % relevant integer. If varid is a number then it is decremented.  This is done
  % because inqnc & getnc count the variables from 1 to nvars whereas the calls
  % to the mexnc routines use c-type counting, from 0 to nvars - 1.

  varid_old = varid;
  if ischar(varid)
    [varid, rcode] = mexnc('ncvarid', cdfid, varid_old);
    if rcode < 0
      mess_str = ['ncvarid: ' full_name ': variable named ' varid_old];
      values = error_handle(cdfid, mess_str, rcode, err_opt);
      return
    end
  else
    varid = varid - 1;
  end

  % Check the value of varid

  if varid < 0 | varid >= nvars
    error([ 'getnc_s was passed varid = ' int2str(varid) ])
  end

  % Find out info about the variable, in particular find nvdims, the number of
  % dimensions that the variable has.

  [varnam, vartypv, nvdims, vdims, nvatts, rcode] = ...
      mexnc('ncvarinq', cdfid, varid);
  if rcode < 0
    mess_str = ['ncvarinq: ' full_name ': variable named ' varid_old];
    values = error_handle(cdfid, mess_str, rcode, err_opt);
    return
  end

  % Turn off the rescaling of the byte type data because mexnc does not do this
  % for variables anyway. The rescaling of the VALUES array will be done
  % explicitly.

  if vartypv == nc_byte
    rescale_var = 0;
    rescale_att = 0;
  end
  
  % Do checks on bl_corner, tr_corner and stride. The cases where
  % bl_corner, tr_corner and stride are -1 or -1*ones(nvdims, 1) are checked
  % for and handled here.  Note that stride may also be a scalar 1 when a
  % vector may seem to be required.

  if nvdims == 0
    bl_corner = 1;
    tr_corner = 1;
    stride = 1;
  else
    if length(bl_corner) == 1
      if bl_corner < 0
	bl_corner = ones(nvdims, 1);
	tr_corner = -1*ones(nvdims, 1);
      end
    elseif length(bl_corner) ~= nvdims
      error('The bl_corner vector is the wrong length')
    end
    
    if length(stride) == 1
      if stride < 0 | stride == 1
	stride = ones(nvdims, 1);
      end
    elseif length(stride) == nvdims
      ff = find(stride < 0);
      if ~isempty(ff)
	stride(ff) = 1;
      end
    else
      error('The stride vector is the wrong length')
    end
  end
  
  bl_corner_min = min(size(bl_corner));
  bl_corner_max = max(size(bl_corner));
  tr_corner_min = min(size(tr_corner));
  tr_corner_max = max(size(tr_corner));
  stride_min = min(size(stride));
  stride_max = max(size(stride));
  if bl_corner_min ~= tr_corner_min | bl_corner_min ~= stride_min | ...
	bl_corner_max ~= tr_corner_max | bl_corner_max ~= stride_max
    error('The sizes of bl_corner, tr_corner and stride do not agree')
  end

  % Set take_stride which specifies whether strides need to be taken.

  if max(stride) > 1
    take_stride = 1;
  else
    take_stride = 0;
  end

  % Make bl_corner, tr_corner, stride and order into column vectors.

  bl_corner = bl_corner(:);
  tr_corner = tr_corner(:);
  stride = stride(:);
  order = order(:);

  % Check order

  if length(order) == 1
    if order == 1 % Special case where the netcdf variable is a vector
      order = -1;
    elseif order ~= -1 & order ~= -2
      error('ERROR: if order is a scalar it must be -1 or -2')
    end
  else
    if length(order) ~= nvdims
      error('The order vector is the wrong length')
    elseif sum(abs(sort(order) - (1:nvdims)')) ~= 0
      error(['The order vector must be a rearrangement of the numbers 1 to ' ...
	     num2str(nvdims)])
    end
  end

  % Find out whether to return a scalar, vector or matrix.  It is here
  % that bl_corner is decremented and edge is calculated so that the c-style
  % conventions in mexnc will be followed. We also check whether any of the
  % dimensions have zero length (this is certainly possible for the record
  % dimension).

  if nvdims == 0
    edge = 1;
  else
    edge = ones(nvdims, 1);
    for i = 1:nvdims
      dimid = vdims(i);
      [name, sizem, rcode] = mexnc('ncdiminq', cdfid, dimid);
      if rcode < 0
	mess_str = ['ncdiminq: ' full_name ': dimid = ' dimid];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
      
      % Check for a zero length dimension (this is certainly possible for the
      % record dimension). If so then just return an empty array.
      
      if sizem == 0
	mess_str = [varnam ' has dimension ' name ' which has zero length'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
      
      % Apply the deaults if required.
      
      if bl_corner(i) < 0
	bl_corner(i) = 1;
	tr_corner(i) = sizem;
	stride(i) = 1;
      end
      if tr_corner(i) < 0
	tr_corner(i) = sizem;
      end
      
      % Check that bl_corner & tr_corner are in the correct range.  If they
      % are then calculate edge.  Note that because I am using the
      % matlab & fortran conventions for counting indices I must
      % subtract 1 from the bl_corner and end point values.
	
      bl_corner(i) = bl_corner(i) - 1;
      tr_corner(i) = tr_corner(i) - 1;
      if bl_corner(i) >= sizem | tr_corner(i) < 0 | tr_corner(i) >= sizem
	s = [ 'getnc_s was passed bl_corner = ' int2str(bl_corner(i)+1) ...
	      ' & tr_corner = ' int2str(tr_corner(i)+1) ...
	      ' for dimension ' name ];
	error(s)
      end
      if stride(i) > 1
	edge(i) = fix( ( tr_corner(i) - bl_corner(i) )/stride(i) ) + 1;
      else
	edge(i) = tr_corner(i) - bl_corner(i) + 1;
      end
    end
  end
  num_edge = length( find(edge ~= 1) );

  if num_edge == 0

    % Get the scalar.

    [values, rcode] = mexnc('ncvarget1', cdfid, varid, bl_corner, rescale_var);
    if rcode < 0
      mess_str = ['ncvarget1: ' full_name ': varid =  ' varid];
      values = error_handle(cdfid, mess_str, rcode, err_opt);
      return
    end
    
    % Do possible byte correction.
    
    if vartypv == nc_byte
      ff = find(values > 127);
      if ~isempty(ff)
	values(ff) = values(ff) - 256;
      end
    end

  else

    % Get the full hyperslab and return it as an array of the appropriate
    % dimensions.  Note that we must allow for the C-type notation where
    % the fastest changing index is the last mentioned.

    if take_stride
      [values, rcode] = mexnc('ncvargetg', cdfid, varid, bl_corner, edge, ...
			      stride, imap, rescale_var);
      if rcode < 0
	mess_str = ['ncvargetg: ' full_name ': varid =  ' varid];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
    else
      [values, rcode] = mexnc('ncvarget', cdfid, varid, bl_corner, edge, ...
			      rescale_var);
      if rcode < 0
	mess_str = ['ncvarget: ' full_name ': varid =  ' varid];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
    end
    
    % Do possible byte correction.
    
    if vartypv == nc_byte
      ff = find(values > 127);
      if ~isempty(ff)
	values(ff) = values(ff) - 256;
      end
    end

    % Permute the array as required.  Note that the default behaviour is to
    % reverse the order of the indices to map between the matlab and C
    % conventions for ordering indices. If order is a vector then this
    % reversing must be done before the permutation using the vector.
    
    if length(order) == 1
      if order == -1
	values = permute(values, (ndims(values):-1:1));
      end
    else
      values = permute(values, (ndims(values):-1:1));
      values = permute(values, order);
    end
    
    % Squeeze the array if required.
    
    if squeeze_it
      values = squeeze(values);
    end

    % After squeezing a vector may be a row or column vector and so
    % turn any row vector into a column vector for consistency.
    
    if ndims(values) == 2
      [m_temp, n_temp] = size(values);
      if m_temp == 1
	values = values(:);
      end
    end
  end

  % If the missing values are to be replaced then do it here.

  scalef = [];
  addoff = [];
  if change_miss ~= 1

    % Find any scale factors or offsets.

    attstring = fill_att(cdfid, varid, nvatts);
    if rescale_att == 1 | vartypv == nc_byte
      pos = check_st('scale_factor', attstring, nvatts);
      if pos > 0
	[scalef, rcode] = mexnc('ncattget', cdfid, varid, 'scale_factor');
	if rcode < 0
	  mess_str = ['ncattget: ' full_name ': varid =  ' varid ...
		      ': scale_factor'];
	  values = error_handle(cdfid, mess_str, rcode, err_opt);
	  return
	end
      end
      pos = check_st('add_offset', attstring, nvatts);
      if pos > 0
	[addoff, rcode] = mexnc('ncattget', cdfid, varid, 'add_offset');
	if rcode < 0
	  mess_str = ['ncattget: ' full_name ': varid =  ' varid ...
		      ': add_offset'];
	  values = error_handle(cdfid, mess_str, rcode, err_opt);
	  return
	end
      end
    end

    % check for missing values.  Note that a
    % missing value is taken to be one less than valid_min, greater than
    % valid_max or 'close to' _FillValue or missing_value.
    % Note 1: valid_min and valid_max may be specified by the attribute
    % valid_range and if valid_range exists than the existence of
    % valid_min and valid_max is not checked.
    % Note 2: a missing value must be OUTSIDE the valid range to be
    % recognised.
    % Note 3: a range does not make sense for character arrays.
    % Note 4: By 'close to' _FillValue I mean that an integer or character
    % must equal _FillValue and a real must be in the range
    % 0.99999*_FillValue tp 1.00001*_FillValue.  This allows real*8 
    % rounding errors in moving the data from the netcdf file to matlab;
    % these errors do occur although I don't know why given that matlab
    % works in double precision.
    % Note 5: An earlier version of this software checked for an attribute
    % named missing_value.  This check was taken out because,
    % although in common use, missing_value was not given in the netCDF
    % manual list of attribute conventions.  Since it has now appeared in
    % the netCDF manual I have put the check back in.
    
    % The indices of the data points containing missing value indicators
    % will be stored separately in index_miss_low, index_miss_up, 
    % index_missing_value and index__FillValue.
    
    index_miss_low = [];
    index_miss_up = [];
    index__FillValue = [];
    index_missing_value = [];
    miss_low_orig = [];
    miss_up_orig = [];
    fill_value_orig = [];
    
    % First find the indices of the data points that are outside the valid
    % range.
    
    pos_vr = check_st('valid_range', attstring, nvatts);
    if pos_vr > 0
      [attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, 'valid_range');
      if rcode < 0
	mess_str = ['ncattinq: ' full_name ': varid =  ' varid ': valid_range'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
      [miss, rcode] = mexnc('ncattget', cdfid, varid, 'valid_range');
      if rcode < 0
	mess_str = ['ncattget: ' full_name ': varid =  ' varid ': valid_range'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
      
      % Check that valid_range is a 2 element vector.
      
      if length(miss) ~= 2
	error(['The valid_range attribute must be a vector'])
      end
      
      % Correct for possible faulty handling of byte type
      
      if attype == nc_byte
	if miss(1) > 127; miss(1) = miss(1) - 256; end
	if miss(2) > 127; miss(2) = miss(2) - 256; end
      end

      miss_low = miss(1);
      miss_up = miss(2);
      miss_low_orig = miss_low;
      miss_up_orig = miss_up;
      
      % Rescale & add offsets if required.
      
      if rescale_att == 1
	if isempty(scalef) == 0
	  miss_low = miss_low*scalef;
	  miss_up = miss_up*scalef;
	end
	if isempty(addoff) == 0
	  miss_low = miss_low + addoff;
	  miss_up = miss_up + addoff;
	end
      end
      
      index_miss_low = find ( values < miss_low );
      index_miss_up = find ( values > miss_up );
    else
      pos_min = check_st('valid_min', attstring, nvatts);
      if pos_min > 0
	[attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, 'valid_min');
	if rcode < 0
	  mess_str = ['ncattinq: ' full_name ': varid =  ' varid ': vvalid_min'];
	  values = error_handle(cdfid, mess_str, rcode, err_opt);
	  return
	end
	[miss_low, rcode] = mexnc('ncattget', cdfid, varid, 'valid_min');
	if rcode < 0
	  mess_str = ['ncattget: ' full_name ': varid =  ' varid ': valid_min'];
	  values = error_handle(cdfid, mess_str, rcode, err_opt);
	  return
	end

	% Check that valid_min is a scalar
	
	if length(miss_low) ~= 1
	  error(['The valid_min attribute must be a scalar'])
	end
	
	% Correct for possible faulty handling of byte type
	
	if attype == nc_byte
	  if miss_low > 127; miss_low = miss_low - 256; end
	end
	miss_low_orig = miss_low;

	% Rescale & add offsets if required.
	
	if rescale_att == 1
	  if isempty(scalef) == 0
	    miss_low = miss_low*scalef;
	  end
	  if isempty(addoff) == 0
	    miss_low = miss_low + addoff;
	  end
	end
	
	index_miss_low = find ( values < miss_low );
      end
      
      pos_max = check_st('valid_max', attstring, nvatts);
      if pos_max > 0
	[attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, 'valid_max');
	if rcode < 0
	  mess_str = ['ncattinq: ' full_name ': varid =  ' varid ': vvalid_max'];
	  values = error_handle(cdfid, mess_str, rcode, err_opt);
	  return
	end
	[miss_up, rcode] = mexnc('ncattget', cdfid, varid, 'valid_max');
	if rcode < 0
	  mess_str = ['ncattget: ' full_name ': varid =  ' varid ': valid_max'];
	  values = error_handle(cdfid, mess_str, rcode, err_opt);
	  return
	end

	% Check that valid_max is a scalar
	
	if length(miss_up) ~= 1
	  error(['The valid_max attribute must be a scalar'])
	end
	
	% Correct for possible faulty handling of byte type
	
	if attype == nc_byte
	  if miss_up > 127;
	    miss_up = miss_up - 256;
	  end
	end
	miss_up_orig = miss_up;

	% Rescale & add offsets if required.
	
	if rescale_att == 1
	  if isempty(scalef) == 0
	    miss_up = miss_up*scalef;
	  end
	  if isempty(addoff) == 0
	    miss_up = miss_up + addoff;
	  end
	end
	
	index_miss_up = find ( values > miss_up );
      end
    end
    
    % Now find the indices of the data points that are 'close to'
    % _FillValue.  Note that 'close to' is different according to the
    % data type.
    
    pos_missv = check_st('_FillValue', attstring, nvatts);
    if pos_missv > 0
      [attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, '_FillValue');
      if rcode < 0
	mess_str = ['ncattinq: ' full_name ': varid =  ' varid ': _FillValue'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
      [miss_val, rcode] = mexnc('ncattget', cdfid, varid, '_FillValue');
      if rcode < 0
	mess_str = ['ncattget: ' full_name ': varid =  ' varid ': _FillValue'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end

      % Check that _FillValue is a scalar
      
      if length(miss_val) ~= 1
	error(['The _FillValue attribute must be a scalar'])
      end
      
      % Correct for possible faulty handling of byte type
      
      if attype == nc_byte
	if miss_val > 127;
	  miss_val = miss_val - 256;
	end
      end
      fill_value_orig = miss_val;
      
      % Check whether _FillValue is outside the valid range to decide
      % whether to keep going.
      
      keep_going = 1;
      if ~isempty(miss_low_orig)
	if (miss_val < miss_low_orig )
	  keep_going = 0;
	end
      end
      if ~isempty(miss_up_orig)
	if (miss_val > miss_up_orig )
	  keep_going = 0;
	end
      end
      
      if keep_going == 1
	
	% Rescale & add offsets if required.
	
	if rescale_att == 1
	  if isempty(scalef) == 0
	    miss_val = miss_val*scalef;
	  end
	  if isempty(addoff) == 0
	    miss_val = miss_val + addoff;
	  end
	end
	
	if attype == nc_byte | attype == nc_char
	  index__FillValue = find ( values == miss_val );
	elseif attype == nc_short | attype == nc_long
	  need_index_m = 1;
	  if pos_vr > 0 | pos_min > 0
	    if miss_val < miss_low
	      need_index_m = 0;
	    end
	  end
	  if pos_vr > 0 | pos_max > 0
	    if miss_val > miss_up
	      need_index_m = 0;
	    end
	  end
	  if need_index_m
	    index__FillValue = find ( values == miss_val );
	  end
	elseif attype == nc_float | attype == nc_double
	  need_index_m = 1;
	  if miss_val < 0
	    miss_val_low = 1.00001*miss_val;
	    miss_val_up = 0.99999*miss_val;
	  else
	    miss_val_low = 0.99999*miss_val;
	    miss_val_up = 1.00001*miss_val;
	  end
	  
	  if pos_vr > 0 | pos_min > 0
	    if miss_val_up < miss_low
	      need_index_m = 0;
	    end
	  end
	  if pos_vr > 0 | pos_max > 0
	    if miss_val_low > miss_up
	      need_index_m = 0;
	    end
	  end
	  if need_index_m
	    index__FillValue = find ( miss_val_low <= values & ...
				      values <= miss_val_up );
	  end
	end
      end
    end
    
    % Now find the indices of the data points that are 'close to'
    % missing_value.  Note that 'close to' is different according to the
    % data type.  This is only done if the missing_value exists and is
    % different to the _FillValue
    
    pos_missv = check_st('missing_value', attstring, nvatts);
    if pos_missv > 0
      [attype, attlen, rcode] = mexnc('ncattinq', cdfid, varid, 'missing_value');
      if rcode < 0
	mess_str = ['ncattinq: ' full_name ': varid =  ' varid ': missing_value'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end
      [miss_val, rcode] = mexnc('ncattget', cdfid, varid, 'missing_value');
      if rcode < 0
	mess_str = ['ncattget: ' full_name ': varid =  ' varid ': missing_value'];
	values = error_handle(cdfid, mess_str, rcode, err_opt);
	return
      end

      % Check that missing_value is a scalar
      
      if length(miss_val) ~= 1
	error(['The missing_value attribute must be a scalar'])
      end
      
      % Correct for possible faulty handling of byte type
      
      if attype == nc_byte
	if miss_val > 127; miss_val = miss_val - 256; end
      end
      
      % Check whether missing_value is outside the valid range to decide
      % whether to keep going.  Also check whether it equals the original
      % _FillValue.
      
      keep_going = 1;
      if ~isempty(miss_low_orig)
	if (miss_val < miss_low_orig)
	  keep_going = 0;
	end
      end
      if ~isempty(miss_up_orig)
	if (miss_val > miss_up_orig)
	  keep_going = 0;
	end
      end
      if ~isempty(fill_value_orig)
	if (miss_val == fill_value_orig)
	  keep_going = 0;
	end
      end
      
      if keep_going == 1
	
	% Rescale & add offsets if required.
	
	if rescale_att == 1
	  if isempty(scalef) == 0
	    miss_val = miss_val*scalef;
	  end
	  if isempty(addoff) == 0
	    miss_val = miss_val + addoff;
	  end
	end
	
	if attype == nc_byte | attype == nc_char
	  index_missing_value = find ( values == miss_val );
	elseif attype == nc_short | attype == nc_long
	  need_index_m = 1;
	  if pos_vr > 0 | pos_min > 0
	    if miss_val < miss_low
	      need_index_m = 0;
	    end
	  end
	  if pos_vr > 0 | pos_max > 0
	    if miss_val > miss_up
	      need_index_m = 0;
	    end
	  end
	  if need_index_m
	    index_missing_value = find ( values == miss_val );
	  end
	elseif attype == nc_float | attype == nc_double
	  need_index_m = 1;
	  if miss_val < 0
	    miss_val_low = 1.00001*miss_val;
	    miss_val_up = 0.99999*miss_val;
	  else
	    miss_val_low = 0.99999*miss_val;
	    miss_val_up = 1.00001*miss_val;
	  end
	  
	  if pos_vr > 0 | pos_min > 0
	    if miss_val_up < miss_low
	      need_index_m = 0;
	    end
	  end
	  if pos_vr > 0 | pos_max > 0
	    if miss_val_low > miss_up
	      need_index_m = 0;
	    end
	  end
	  if need_index_m
	    index_missing_value = find ( miss_val_low <= values & ...
					 values <= miss_val_up );
	  end
	end
      end
    end
    
    % Combine the arrays of missing value indices into one unordered array.
    % Note that for real numbers the range of the _FillValue and
    % missing_value may intersect both the valid and invalid range and so
    % some indices may appear twice; this does not cause any inaccuracy,
    % although it will result in some inefficiency.  In particular,
    % rescaling is done on the set of indices NOT in index_miss and so is
    % not affected.
    
    index_miss = [ index_miss_low(:); index__FillValue(:); ...
		   index_missing_value(:); index_miss_up(:) ];
    %index_miss = sort(index_miss);
    len_index_miss = length(index_miss);
    
    % If there are any missing values then change them to a
    % more convenient value.
    
    if len_index_miss > 0
      if change_miss == 2
	if vartypv == nc_char
	  values(index_miss) = char(0);
	else
	  values(index_miss) = NaN;
	end
      elseif change_miss == 3
	if vartypv == nc_char
	  if isnumeric(new_miss)
	    values(index_miss) = char(new_miss);
	  else
	    values(index_miss) = new_miss;
	  end
	else
	  values(index_miss) = new_miss;
	end
      else
	s = [ 'getnc_s was passed change_miss = ' int2str(change_miss) ];
	error(s)
      end
    end
  end

  % Rescale the byte type data which was not done automatically. If the option
  % to not rescale has been selected then scalef and addoff will be empty and
  % there will be no rescaling.

  if vartypv == nc_byte
    if isempty(scalef) == 0
      values = values*scalef;
    end
    if isempty(addoff) == 0
      values = values + addoff;
    end
  end

  % Close the netcdf file.

  [rcode] = mexnc('ncclose', cdfid);
  if rcode < 0
    mess_str = ['ncclose: ' full_name];
    values = error_handle(cdfid, mess_str, rcode, err_opt);
    return
  end

 case {'loaddap', 'loaddods'}
  
  % Dealing with a dods file.
  
  % Check the value of varid and get its structure.
  
  if ischar(varid)
    if ~isfield(desc_das, varid)
      mess_str = [varid ' is not a variable in ' file];
      values = error_handle([], mess_str, [], err_opt);
      return
    end
  else
    error('For dods data varid must be a string')
  end
  
  % Get the structure of varid and then use this to extract required
  % information about attributes. Also get the dds of the file. Note that we
  % use the function remove_quotes to eliminate extraneous quotes returned as
  % part of the DAS. 
  
  %desc_varid = desc_das.(varid);
  desc_varid = getfield(desc_das, varid);
  if isfield(desc_varid, 'scale_factor')
    %scalef = desc_varid.('scale_factor');
    scalef = getfield(desc_varid, 'scale_factor');
    scalef = remove_quotes(scalef);
    if length(scalef) ~= 1
      error(['The scalef attribute must be a scalar'])
    end
  else
    scalef = [];
  end
  if isfield(desc_varid, 'add_offset')
    %addoff = desc_varid.('add_offset');
    addoff = getfield(desc_varid, 'add_offset');
    addoff =  remove_quotes(addoff);
    if length(addoff) ~= 1
      error(['The addoff attribute must be a scalar'])
    end
 else
    addoff = [];
  end
  if isfield(desc_varid, 'ml__FillValue')
    %fillvalue = desc_varid.('ml__FillValue');
    fillvalue = getfield(desc_varid, 'ml__FillValue');
    fillvalue =  remove_quotes(fillvalue);
    if rescale_att
      if ~isempty(scalef)
	fillvalue = fillvalue*scalef;
      end
      if ~isempty(addoff)
	fillvalue = fillvalue + addoff;
      end
    end
  else
    fillvalue = [];
  end
  if isfield(desc_varid, 'missing_value')
    %missing_value = desc_varid.('missing_value');
    missing_value = getfield(desc_varid, 'missing_value');
    missing_value =  remove_quotes(missing_value);
    if rescale_att
      if ~isempty(scalef)
	missing_value = missing_value*scalef;
      end
      if ~isempty(addoff)
	missing_value = missing_value + addoff;
      end
    end
  else
    missing_value = [];
  end
  if isfield(desc_varid, 'valid_min')
    %valid_min = desc_varid.('valid_min');
    valid_min = getfield(desc_varid, 'valid_min');
    valid_min =  remove_quotes(valid_min);
    if rescale_att
      if ~isempty(scalef)
	valid_min = valid_min*scalef;
      end
      if ~isempty(addoff)
	valid_min = valid_min + addoff;
      end
    end
  else
    valid_min = [];
  end
  if isfield(desc_varid, 'valid_max')
    %valid_max = desc_varid.('valid_max');
    valid_max = getfield(desc_varid, 'valid_max');
    valid_max =  remove_quotes(valid_max);
    if rescale_att
      if ~isempty(scalef)
	valid_max = valid_max*scalef;
      end
      if ~isempty(addoff)
	valid_max = valid_max + addoff;
      end
    end
  else
    valid_max = [];
  end
  if isfield(desc_varid, 'valid_range')
    %valid_range = desc_varid.('valid_range');
    valid_range = getfield(desc_varid, 'valid_range');
    valid_range =  remove_quotes(valid_range);
    if rescale_att
      if ~isempty(scalef)
	valid_range = valid_range*scalef;
      end
      if ~isempty(addoff)
	valid_range = valid_range + addoff;
      end
    end
  else
    valid_range = [];
  end
    
  % Get the dds of the file and use it to get some required information.
  
  [dds_text, desc_dds] = get_dods_dds(file, exe_name);
  
  for ii = 1:length(desc_dds.variable)
    if strcmp(desc_dds.variable(ii).name, varid)
      break
    end
  end
  if ~strcmp(desc_dds.variable(ii).name, varid)
    error('This error should not ever occur')
  end
  var_type = desc_dds.variable(ii).type;
  var_dim_statement = desc_dds.variable(ii).dim_statement;
  var_dim_idents = desc_dds.variable(ii).dim_idents;
  
  nvdims = length(var_dim_idents);
  dim_length = zeros(nvdims, 1);
  for ii = 1:nvdims
    dim_length(ii) = desc_dds.dimension(var_dim_idents(ii)).length;
    dim_name{ii} = desc_dds.dimension(var_dim_idents(ii)).name;
  end
  
  % Do checks on bl_corner, tr_corner and stride. The cases where
  % bl_corner, tr_corner and stride are -1 or -1*ones(nvdims, 1) are checked
  % for and handled here.  Note that stride may also be a scalar 1 when a
  % vector may seem to be required.

  if nvdims == 0
    bl_corner = 1;
    tr_corner = 1;
    stride = 1;
  else
    if length(bl_corner) == 1
      if bl_corner < 0
	bl_corner = ones(nvdims, 1);
	tr_corner = -1*ones(nvdims, 1);
      end
    elseif length(bl_corner) ~= nvdims
      error('The bl_corner vector is the wrong length')
    end
    
    if length(stride) == 1
      if stride < 0 | stride == 1
	stride = ones(nvdims, 1);
      end
    elseif length(stride) == nvdims
      ff = find(stride < 0);
      if ~isempty(ff)
	stride(ff) = 1;
      end
    else
      error('The stride vector is the wrong length')
    end
  end

  bl_corner_min = min(size(bl_corner));
  bl_corner_max = max(size(bl_corner));
  tr_corner_min = min(size(tr_corner));
  tr_corner_max = max(size(tr_corner));
  stride_min = min(size(stride));
  stride_max = max(size(stride));
  if bl_corner_min ~= tr_corner_min | bl_corner_min ~= stride_min | ...
	bl_corner_max ~= tr_corner_max | bl_corner_max ~= stride_max
    error('The sizes of bl_corner, tr_corner and stride do not agree')
  end

  % Set take_stride.

  if max(stride) > 1
    take_stride = 1;
  else
    take_stride = 0;
  end

  % Make bl_corner, tr_corner, stride and order into column vectors.

  bl_corner = bl_corner(:);
  tr_corner = tr_corner(:);
  stride = stride(:);
  order = order(:);

  % Check order

  if length(order) == 1
    if order == 1 % Special case where the netcdf variable is a vector
      order = -1;
    elseif order ~= -1 & order ~= -2
      error('ERROR: the argument order must be -1 or -2')
    end
  else
    error(['ERROR: the argument order must be a scalar - the vector ' ...
	   ' option is not supported for opendap access'])
  end

  % Generate a string that specifies the hyperslab to be returned. In the
  % process check that the hyperslab is within bounds.
  
  str_h = [];
  for ii = 1:nvdims
    if bl_corner(ii) < 0
      bl_corner(ii) = 1;
      tr_corner(ii) = dim_length(ii);
      stride(ii) = 1;
    elseif bl_corner(ii) > dim_length(ii)
      error(['bl_corner so big that it is outside the hyperslab'])
    end
    if tr_corner(ii) < 0
      tr_corner(ii) = dim_length(ii);
    elseif tr_corner(ii) > dim_length(ii)
      error(['tr_corner so big that it is outside the hyperslab'])
    elseif tr_corner(ii) < bl_corner(ii)
      error(['tr_corner is less than bl_corner'])
    end
    
    if take_stride & (stride(ii) > 1)
      str_h = [str_h '[' num2str(bl_corner(ii) - 1) ':' num2str(stride(ii)) ...
	       ':' num2str(tr_corner(ii) - 1) ']'];
    else
      str_h = [str_h '[' num2str(bl_corner(ii) - 1) ':' ...
	       num2str(tr_corner(ii) - 1) ']'];
    end
  end
  
  % Get the hyperslab of data. Note that we allow for the vector of variable
  % values to be at different "depths" in the structure. It seems that a
  % variable that is also a dimension will be at values_struct.varid but
  % otherwise it will be at values_struct.varid.varid. We also use
  % method_of_call == 2 which, according to the documentation, does the same
  % thing as method_of_call == 1 but would be expected to be slower. However,
  % this is necessary for loaddods to work properly on Windows boxes when
  % an array is filled with characters. I have no idea what the bug might be
  % in loaddods. There is also a hint that the same problem can occur with
  % loaddap on some linux boxes.
  
  method_of_call = 2;
  switch method_of_call
   case 1
    switch mex_name
     case 'loaddap'
      values_struct = loaddap('+v', [full_name '?' varid str_h]);
     case 'loaddods'
      values_struct = loaddods('+v', [full_name '?' varid str_h]);
    end
   case 2
    switch mex_name
     case 'loaddap'
      loaddap('+v', [full_name '?' varid str_h]);
     case 'loaddods'
      loaddods('+v', [full_name '?' varid str_h]);
    end
    eval(['values_struct = ' varid ';'])
  end
  if isstruct(values_struct) 
    %xx = values_struct.(varid);
    xx = getfield(values_struct, varid);
    if isstruct(xx)
      %values = xx.(varid);
      values = getfield(xx, varid);
    else
      values = xx;
    end
  else
    values = values_struct;
  end
  
  % If we are getting an array of character data then loaddap returns a row
  % vector regardless of the original format. The following kluge is designed
  % to recreate the original array structure.
  
  kluge_char = 1;
  if kluge_char & ischar(values)
    dim_new = zeros(1, nvdims);
    for kk = 1:nvdims
      dim_new(nvdims - kk + 1) = length(bl_corner(kk):stride(kk):tr_corner(kk));
    end
    if nvdims == 1
      values = reshape(values, [1 dim_new]);
    elseif nvdims == 2
      values = reshape(values, dim_new);
    else
      values = reshape(values, dim_new);
      %error('Can''t handle code with a character array of more than 2 dimensions')
    end
  end
  
  % Squeeze the array if required.
  
  if squeeze_it
    values = squeeze(values);
  end

  % Now deal with permutation problems in the data.
  
  size_values = size(values);
  num_dim = length(size_values);
  if num_dim == 2
    if min(size_values) == 1
      % After squeezing a vector may be a row or column vector and so
      % turn any row vector into a column vector for consistency.
      
      values = values(:);
    else
      % Apply the order option to transpose the matrix
      if kluge_char & ischar(values)
	if order == -1
	  values = values';
	end
      else
	if order == -2
	  values = values';
	end
      end
    end
  elseif num_dim > 2
    
    % Permute the array as required.  There are 3 issues here.
    % 1) loaddap returns multi-dimensional numerical arrays in a strange order -
    % half-fortran and half-C conventions. First explain this conceptually. We
    % see the variable dimensions in a certain order when we look at the DDS -
    % the order is the same using inqnc or ncdump on a netcdf file. After
    % loaddap puts the data into a matlab variable it looks as though there
    % have been 2 transformations. First the order is totally reversed and
    % then the order of the first 2 dimensions are reversed
    % again. Schematically, suppose that a variable has N dimensions labelled
    % 1, 2, 3, ..., N-4, N-3, N-2, N-1, N in the DDS. The matlab order is then
    % N-1, N, N-2, N-3, N-4, ..., 3, 2, 1.
    % 2) loaddap returns multi-dimensional character arrays as vectors. It
    % simply streams out the data in the same order that it is stored on
    % disk.
    % 3) We must permute the array according to the value of the argument
    % "order" also.
    
    if kluge_char & ischar(values)
      if order == -1
	vec_permute = [num_dim:-1:1];
      elseif order == -2
	vec_permute = [1:num_dim];
      end
    else
      if order == -1
	vec_permute = [num_dim:-1:3 1 2];
      elseif order == -2
	vec_permute = [2 1 3:num_dim];
      end
    end
    values = permute(values, vec_permute);
  end

  % Do possible byte correction.
    
  if strcmp(lower(var_type), 'byte')
    ff = find(values > 127);
    if ~isempty(ff)
      values(ff) = values(ff) - 256;
    end
  end

  % Possibly do a rescale of the variable.
  
  if rescale_var
    if ~isempty(scalef)
      values = values*scalef;
    end
    if ~isempty(addoff)
      values = values + addoff;
    end
  end

  % check for missing values.  Note that a
  % missing value is taken to be one less than valid_min, greater than
  % valid_max or 'close to' _FillValue or missing_value.
  % Note 1: valid_min and valid_max may be specified by the attribute
  % valid_range and if valid_range exists than the existence of
  % valid_min and valid_max is not checked.
  % Note 2: a missing value must be OUTSIDE the valid range to be
  % recognised.
  % Note 3: a range does not make sense for character arrays.
  % Note 4: By 'close to' _FillValue I mean that an integer or character
  % must equal _FillValue and a real must be in the range
  % 0.99999*_FillValue tp 1.00001*_FillValue.  This allows real*8 
  % rounding errors in moving the data from the netcdf file to matlab;
  % these errors do occur although I don't know why given that matlab
  % works in double precision.
  % Note 5: An earlier version of this software checked for an attribute
  % named missing_value.  This check was taken out because,
  % although in common use, missing_value was not given in the netCDF
  % manual list of attribute conventions.  Since it has now appeared in
  % the netCDF manual I have put the check back in.
  
  % The indices of the data points containing missing value indicators
  % will be stored separately in index_miss_low, index_miss_up, 
  % index_missing_value and index_fillvalue.
  
  index_miss_low = [];
  index_miss_up = [];
  index_fillvalue = [];
  index_missing_value = [];
  miss_low_orig = [];
  miss_up_orig = [];
  fill_value_orig = [];
  if ~isempty(findstr(lower(var_type), 'float')) | ...
	~isempty(findstr(lower(var_type), 'int'))
    var_is_number = 1;
  else
    var_is_number = 0;
  end

  % First find the indices of the data points that are outside the valid
  % range.
    
  if ~isempty(valid_range)
    
    % Check that valid_range is a 2 element vector.
    
    if length(valid_range) == 2
      miss = valid_range;
    else
      error(['The valid_range attribute must be a 2 element vector'])
    end
    
    miss_low = miss(1);
    miss_up = miss(2);
    miss_low_orig = miss_low;
    miss_up_orig = miss_up;
    
    index_miss_low = find (values < miss_low);
    index_miss_up = find (values > miss_up);
  else
    if ~isempty(valid_min)
      miss_low = valid_min;

      % Check that valid_min is a scalar
	
      if length(miss_low) ~= 1
	error(['The valid_min attribute must be a scalar'])
      end
      index_miss_low = find ( values < miss_low );
    end
    
    if ~isempty(valid_max)
      miss_up = valid_max;
  
      % Check that valid_max is a scalar
	
      if length(miss_up) ~= 1
	error(['The valid_max attribute must be a scalar'])
      end
      index_miss_up = find ( values > miss_up );
    end
  end
  
  % Now find the indices of the data points that are the same as the
  % _FillValue.
  
  if ~isempty(fillvalue)
    miss_val = fillvalue;
    
    % Check that _FillValue is a scalar
    
    if length(miss_val) ~= 1
      error(['The _FillValue attribute must be a scalar'])
    end
 
    % Check whether _FillValue is outside the valid range to decide
    % whether to keep going.
    
    keep_going = 1;
    if ~isempty(miss_low_orig)
      if (miss_val < miss_low_orig)
	keep_going = 0;
      end
    end
    if ~isempty(miss_up_orig)
      if (miss_val > miss_up_orig)
	keep_going = 0;
      end
    end
    
    if keep_going
      if var_is_number
	if miss_val < 0
	  miss_val_low = 1.00001*miss_val;
	  miss_val_up = 0.99999*miss_val;
	else
	  miss_val_low = 0.99999*miss_val;
	  miss_val_up = 1.00001*miss_val;
	end
	index_fillvalue = find((miss_val_low <= values) & ...
				(values <= miss_val_up));
      else
	index_fillvalue = find(values == miss_val);
      end
    end
  end

  % Now find the indices of the data points that are the same as the
  % missing_value.
  
  if ~isempty(missing_value)
    miss_val = missing_value;
    
    % Check that missing_value is a scalar
    
    if length(miss_val) ~= 1
      error(['The missing_value attribute must be a scalar'])
    end
 
    % Check whether missing_value is outside the valid range to decide
    % whether to keep going. Also, don't bother doing the search if
    % _FillValue is identical to missing_value.
    
      keep_going = 1;
      if ~isempty(miss_low_orig)
	if (miss_val < miss_low_orig )
	  keep_going = 0;
	end
      end
      if ~isempty(miss_up_orig)
	if (miss_val > miss_up_orig )
	  keep_going = 0;
	end
      end
    if ~isempty(fillvalue)
      if (missing_value == fillvalue)
	keep_going = 0;
      end
    end
    
    if keep_going
      if var_is_number
	if miss_val < 0
	  miss_val_low = 1.00001*miss_val;
	  miss_val_up = 0.99999*miss_val;
	else
	  miss_val_low = 0.99999*miss_val;
	  miss_val_up = 1.00001*miss_val;
	end
	index_missing_value = find((miss_val_low <= values) & ...
				(values <= miss_val_up));
      else
	index_missing_value = find(values == miss_val);
      end
    end
  end
  
  % Combine the arrays of missing value indices into one unordered array.
  % Note that for real numbers the range of the _FillValue and
  % missing_value may intersect both the valid and invalid range and so
  % some indices may appear twice; this does not cause any inaccuracy,
  % although it will result in some inefficiency.
  
  index_miss = [ index_miss_low(:); index_fillvalue(:); ...
		 index_missing_value(:); index_miss_up(:) ];
  len_index_miss = length(index_miss);
  
  if len_index_miss > 0
    if change_miss == 2
      if isnumeric(values)
	values(index_miss) = NaN;
      else
	values(index_miss) = char(0);
      end
    elseif change_miss == 3
      if isnumeric(values)
	values(index_miss) = new_miss;
      else
	if isnumeric(new_miss)
	  values(index_miss) = char(new_miss);
	else
	  values(index_miss) = new_miss;
	end
      end
    elseif change_miss ~= 1
      s = [ 'getnc_s was passed change_miss = ' int2str(change_miss) ];
      error(s)
    end
  end
 case 'java'
  
  % Translate rescale_opts to a 2 element vector.
  
  if length(rescale_opts) == 1
    if rescale_opts == -1
      rescale_opts = [1 1];
    else
      error('rescale_opts must be a 2 element vector or -1')
    end
  end
      
  % Open dataset so that this script will have to do the
  % scale/offset/missing value stuff. This is because the java
  % interface does not deal with missing values properly for character
  % arrays and rescale_opts may force some things to be rescaled and others
  % not. Note that the 2nd argument relates to enhancing
  % the data, i.e., automatically doing scale/offset/missing stuff. The 3rd
  % argument relates to error handling and hopefully is the way to pass a
  % null to the constructor.

  try
    enhance = false;
    ncdJ = ucar.nc2.dataset.NetcdfDataset.openDataset(file, enhance, []);
  catch
    ss = lasterror;
    mess_str = ss.message;
    rcode = -1000000;
    values = error_handle('java', mess_str, rcode, err_opt);
    return
  end
  
  % Get the variable object. Note that this is different according to whether
  % whether varid is numeric or not.
  
  if isnumeric(varid)
    
    % Check that varid is a scalar.
    
    if ndims(varid) > 2
      wrong_size = 1;
    else     
      if all(size(varid) == [1 1])
	wrong_size = 0;
      else
	wrong_size = 1;
      end
    end
    if wrong_size
      mess_str = 'varid must be either a string or a scalar';
      rcode = -1000000;
      values = error_handle(ncdJ, mess_str, rcode, err_opt);
    end
    
    % Check that varid is in the correct range.
    
    variables_list = ncdJ.getVariables();
    nvars = variables_list.size();
    if (varid < 1) || (varid > nvars)
      mess_str = ['Failed trying to access variable number: ' varid];
      rcode = -1000000;
      values = error_handle(ncdJ, mess_str, rcode, err_opt);
      return
    end
    
    % Get the variable object and reset varid to be the name of the variable.
    
    varid_old = varid;
    varJ = variables_list.get(varid_old - 1);
    varid = varJ.getName();
  else
    
    % Try to get the variable object and then check that this succeeded.
    
    varJ = ncdJ.findVariable(varid);
    if isempty(varJ)
      mess_str = ['Failed trying to access variable: ' varid];
      rcode = -1000000;
      values = error_handle(ncdJ, mess_str, rcode, err_opt);
      return
    end
  end
  
  % There are problems translating from netcdf to opendap formats when
  % dealing with character arrays. If varJ.getDataType returns 'char' then an
  % array of char is returned and it can be handled pretty much the same as a
  % numeric array although there are issues with missing values and switching
  % between matlabs numeric and char types. However the opendap data may be
  % seen as an array of strings. In this case varJ.getDataType will be
  % 'String' and things must be handled differently. This is complicated even
  % more by the issue of what is the fastest changing index - hence calls to
  % fliplr and permute.
  
  var_type_str = varJ.getDataType;
  
  % Get the data. We may either get all of the data or only some subset.
  
  if all(bl_corner == -1) && all(tr_corner == -1) && all(stride == -1)
    
    % Get all of variable.
    
    if strcmp(var_type_str, 'String')
      
      % Get the character data as a 1D array. Note all of the fiddling with
      % val_shape_row_orig and val_shape_row_reversed to deal with the
      % difference of fastest changing index between matlab and C.
      
      val_size = varJ.getSize();
      val_rank = varJ.getRank();
      val_shape = varJ.getShape();
      val_shape_row_orig = val_shape(:)';
      val_shape_row_reversed = fliplr(val_shape_row_orig);
      val_int = copyTo1DJavaArray(varJ.read());
      switch class(val_int(1))
       case 'opendap.dap.DString'
	
	% Read the string data allowing for weirdness of the returned
        % value. In particular it may call a null character an empty
        % string. Thus we must read the data one character at a time.
	
	len_val_int = length(val_int);
	values = repmat(' ', len_val_int, 1);
	for ii = 1:len_val_int
	  xx = val_int(ii).getValue();
	  if length(xx) == 0
	    values(ii) = char(0);
	  else
	    values(ii) = char(xx);
	  end
	end
	
	% Turn values into a properly dimensioned array (fastest changing
        % index problem).
	
	if val_rank > 1
	  values = reshape(values, val_shape_row_reversed);
	  values = permute(values, val_rank:-1:1);
	end
	
       otherwise
	error('Cannot handle object')
      end
    else
      values = copyToNDJavaArray(varJ.read());
    end
  else
    
    % Change bl_corner, tr_corner and stride so that they contain the
    % hyperslab description, i.e., replace the -1 values with proper indices
    % and then generate the string readSpec. This contains the fortran-like 
    % information readSpec = '0:7:2,3:8:3';
    
    varRank = varJ.getRank();
    varShape(1:varRank) = varJ.getShape();
    if varRank > 1
      if length(bl_corner) == 1
	bl_corner = repmat(bl_corner(1), varRank, 1);
      end
      if length(tr_corner) == 1
	tr_corner = repmat(tr_corner(1), varRank, 1);
      end
      if length(stride) == 1
	stride = repmat(stride(1), varRank, 1);
      end
    end
    for ii = 1:varRank
      if bl_corner(ii) < 0
	bl_corner(ii) = 1;
	stride(ii) = 1;
	tr_corner(ii) = varShape(ii);
      elseif (bl_corner(ii) == 0) || (bl_corner(ii) > varShape(ii))
	mess_str = '1: hyperslab is badly specified';
	rcode = -1000000;
	values = error_handle(ncdJ, mess_str, rcode, err_opt);
	return
      end
      if tr_corner(ii) == -1
	tr_corner(ii) = varShape(ii);
      elseif (tr_corner(ii) < bl_corner(ii)) || (tr_corner(ii) > varShape(ii))
	mess_str = '2: hyperslab is badly specified';
	rcode = -1000000;
	values = error_handle(ncdJ, mess_str, rcode, err_opt);
	return
      end
    end
    if any(stride < -1) || any(stride == 0)
      mess_str = 'bad stride';
      rcode = -1000000;
      values = error_handle(ncdJ, mess_str, rcode, err_opt);
      return
    end
    stride(find(stride == -1)) = 1;
      
    % Construct the string that specifies the hyperslab and then get the data.
    readSpec = '';
    for ii = 1:varRank
      readSpec = [readSpec num2str(bl_corner(ii) - 1) ':' ...
		  num2str(tr_corner(ii) - 1) ':' ...
		  num2str(stride(ii)) ','];
    end
    readSpec = readSpec(1:(end - 1));
    if strcmp(var_type_str, 'String')
      val_int = copyTo1DJavaArray(varJ.read(readSpec));
      switch class(val_int(1))
       case 'opendap.dap.DString'
	
	% Read the string data allowing for weirdness of the returned
        % value. In particular it may call a null character an empty
        % string. Thus we must read the data one character at a time.
	
	len_val_int = length(val_int);
	values = repmat(' ', len_val_int, 1);
	for ii = 1:len_val_int
	  xx = val_int(ii).getValue();
	  if length(xx) == 0
	    values(ii) = char(0);
	  else
	    values(ii) = char(xx);
	  end
	end
	
	% Turn into properly dimensioned array (fastest changing index
        % problem).
	
	val_rank = varJ.getRank();
	lims_dim = floor((tr_corner(:) - bl_corner(:))./stride(:)) + 1;
	lims_dim = fliplr(lims_dim');
	if val_rank > 1
	  values = reshape(values, lims_dim);
	  values = permute(values, val_rank:-1:1);
	end
	
       otherwise
	error('Cannot handle object')
      end
    else
      values = copyToNDJavaArray(varJ.read(readSpec));
    end
  end
  
  % Make sure that all numeric arrays are returned as double.
  
  if isnumeric(values)
    values = double(values);
  end

  % Possibly squeeze the array.
  
  if squeeze_it
    values = squeeze(values);
  end
  
  % Possibly change the order of the elements according to the value of the
  % variable "order".

  size_val = size(values);
  nd = length(size_val);
  if all(size(order) == 1)
    % order is a scalar
    if order == -2
      values = permute(values, (nd:-1:1));
    end
  else
    if length(order) == nd
      values = permute(values, order);
    else
      error('The order vector does not match the rank of the variable')
    end
  end
  
  % Make sure that all vectors are column vectors.
  
  size_val = size(values);
  nd = length(size_val);
  if (nd == 2) && any(size_val == 1)
    values = values(:);
  end
  
  % Deal with rescaling and change_miss stuff now.
  
  % Get all of the attributes that are required and store them in suitable
  % variables.
  
  att_list = varJ.getAttributes();
  num_atts = att_list.size();
  miss_val_list = [];
  min_val = -Inf;
  max_val = Inf;
  scale_factor = [];
  add_offset = [];
  for ii = 0:(num_atts - 1)
    att = att_list.get(ii);
    att_name =  char(att.getName());
    switch att_name
     case '_FillValue'
      fillvalue = get_attribute_value(att);
      miss_val_list = [miss_val_list fillvalue];
     case 'missing_value'
      missing_value = get_attribute_value(att);
      miss_val_list = [miss_val_list missing_value];
     case 'valid_range'
      valid_range = get_attribute_value(att);
      min_val = max([min_val valid_range(1)]);
      max_val = min([max_val valid_range(2)]);
     case 'valid_min'
      valid_min = get_attribute_value(att);
      min_val = max([min_val valid_min]);
     case 'valid_max'
      valid_max = get_attribute_value(att);
      max_val = min([max_val valid_max]);
     case 'scale_factor'
      scale_factor = get_attribute_value(att);
     case 'add_offset'
      add_offset = get_attribute_value(att);
    end
  end
  
  % Deal with rescaling first
  
  if rescale_opts(1)
    % Rescale the variable.

    if ~isempty(scale_factor)
      values = values*scale_factor;
    end
    if ~isempty(add_offset)
      values = values + add_offset;
    end
  end
  
  if rescale_opts(2)
    % Rescale the variable.

    if ~isempty(scale_factor)
      min_val = min_val*scale_factor;
      max_val = max_val*scale_factor;
      miss_val_list = miss_val_list*scale_factor;
    end
    if ~isempty(add_offset)
      min_val = min_val + add_offset;
      max_val = max_val + add_offset;
      miss_val_list = miss_val_list + add_offset;
    end
  end

  % Now do the missing value calculations.
  
  if change_miss == 1
    replace_missing_values = 0;
  elseif change_miss == 3
    new_miss_val = new_miss;
    replace_missing_values = 1;
  else
    new_miss_val = NaN;
    replace_missing_values = 1;
  end

  if replace_missing_values
    if isnumeric(values)
      for ii = 1:length(miss_val_list)
	ff = find(abs(values - miss_val_list(ii)) < 1e-5);
	values(ff) = new_miss_val;
      end
      values(find(values < min_val)) = new_miss_val;
      values(find(values > max_val)) = new_miss_val;
    else
      if isnumeric(new_miss_val)
	if isnan(new_miss_val)
	  new_miss_val = char(0);
	else
	  new_miss_val = char(new_miss_val);
	end
      end
      si = size(values);
      values = values(:)';
      for ii = 1:length(miss_val_list)
	ff = strfind(values, miss_val_list(ii));
	values(ff) = new_miss_val;
      end
      values = reshape(values, si);
    end
  end
  
  % Make sure that all numeric arrays are returned as double.
  
  if isnumeric(values)
    values = double(values);
  end

  % Close the file object

  ncdJ.close();
 case 'none'
  error(['Couldn''t find a suitable mex-file for reading ' file])
end

function att_val = get_attribute_value(att)
% Get the value of an attribute object.
  
  if att.isString()
    att_val = char(att.getStringValue());
  else
    num = double(att.getLength);
    att_val = zeros(1, num);
    for mm = 1:num
      att_val(mm) = double(att.getNumericValue(mm - 1));
    end
  end

function new_val = remove_quotes(val)
% If val is a string then we replace any extraneous quote marks.
  if ischar(val)
    s = abs(val);
    if (s(1) == s(end)) & ((s(1) == 34) | (s(1) == 39))
      s = s(2:(end - 1));
    end
    new_val = char(s);
  else
    new_val = val;
  end

function values = error_handle(fid, mess_str, rcode, err_opt)

% error_handle is called after a mexnc or java call has failed. It ensures
% that an open netcdf file is closed. The value of err_opt determines what
% else is done. For a mexnc call fid is cdfid, the handle to the open
% file. For a java call fid is the opened file object.
%    err_opt == 0 prints an error message and then aborts
%            == 1 prints a warning message and then returns an empty
%                 array. This is the default.
%            == 2 returns an empty array. This is a very dangerous option and
%                 should only be used with caution. It might be used when
%                 getnc_s is called in a loop and you want to do your own
%                 error handling without being bothered by warning messages.

% Decide what part of the code made the call.
  
  if isempty(fid)
    called_by = 'loadd';
  else
    if isnumeric(fid)
      called_by = 'mexnc';
    else
      called_by = 'java';
    end
  end

% Close an open netcdf of java file.
  
  switch called_by
   case 'mexnc'
    if fid >= 0
      [rcode_sub] = mexnc('ncclose', fid);
    end
   case 'java'
    if isjava(fid)
      fid.close();
    end
  end
  
  % Handle the errors according to the value of err_opt. If rcode is empty
  % then this is probably because loaddap or loaddods was called.
  
  values = [];
  switch err_opt
   case 0
    if isempty(rcode)
      str = ['ERROR: ' mess_str];
    else
      str = ['ERROR: ' mess_str ' : rcode = ' num2str(rcode)];
    end
    error(str)
   case 1
    if isempty(rcode)
      str = ['WARNING: ' mess_str];
    else
      str = ['WARNING: ' mess_str ' : rcode = ' num2str(rcode)];
    end
    disp(str)
   case 2
    return
   otherwise
    error(['error_handle was called with err_opt = ' num2str(err_opt)])
  end

function [file, varid, bl_corner, tr_corner, stride, order, change_miss, ...
	  new_miss, squeeze_it, rescale_opts, err_opt] = parse_args(cell_args);

  % parse_args receives the cell that was passed to getnc via the varargin
  % method. It parses it to get all of the output arguments. It also sets the
  % defaults for these.
  
  % Specify all of the default values.
  
  file_def = []; % 1
  varid_def = []; % 2
  bl_corner_def = -1; % 3
  tr_corner_def = -1; % 4
  stride_def = -1; % 5
  order_def = -1; % 6
  change_miss_def = 2; % 7
  new_miss_def = 0; % 8
  squeeze_it_def = 1; % 9
  rescale_opts_def = -1; % 10
  err_opt_def = 1; % 11
  
  % Set the variables to their default values.
  
  file = file_def;
  varid = varid_def;
  bl_corner = bl_corner_def;
  tr_corner = tr_corner_def;
  stride = stride_def;
  order = order_def;
  change_miss = change_miss_def;
  new_miss = new_miss_def;
  squeeze_it = squeeze_it_def;
  rescale_opts = rescale_opts_def;
  err_opt = err_opt_def;
  
  % Check that we do not have too many arguments passed.
  
  len_cell_args = length(cell_args);
  if len_cell_args > 11
    error([num2str(len_cell_args) ' input arguments (too many) were passed'])
  end
  
  % Work through cell_args, overwriting variables as the input arguments are
  % read.
  % 1) The last element of cell_args may be a structure containing multiple
  % cases of values to be assigned to variables.
  % 2) If the passed value is an empty array or -1 then it is assumed that
  % the user wanted to pass the default. The exception to this is new_miss
  % where -1 would be an acceptable value to use.
  
  for ii = 1:len_cell_args
    if isstruct(cell_args{ii})
      if ii == len_cell_args
	fname = fieldnames(cell_args{ii});
	for jj = 1:length(fname)
	  val = getfield(cell_args{ii}, fname{jj});
	  switch fname{jj}
	   case 'file'
	    file = substitute_new(val, order_def, 1);
	   case 'varid'
	    varid = substitute_new(val, order_def, 1);
	   case 'bl_corner'
	    bl_corner = substitute_new(val, order_def, 1);
	   case 'tr_corner'
	    tr_corner = substitute_new(val, order_def, 1);
	   case 'stride'
	    stride = substitute_new(val, order_def, 1);
	   case 'order'
	    order = substitute_new(val, order_def, 1);
	   case 'change_miss'
	    change_miss = substitute_new(val, change_miss_def, 1);
	   case 'new_miss'
	    new_miss = substitute_new(val, new_miss_def, 0);
	   case 'squeeze_it'
	    squeeze_it = substitute_new(val, squeeze_it_def, 1);
	   case 'rescale_opts'
	    rescale_opts = substitute_new(val, rescale_opts_def, 1);
	   case 'err_opt'
	    err_opt = substitute_new(val, err_opt_def, 1);
	   otherwise
	    error(['Don''t recognise the field ' fname{jj}])
	  end
	end
      else
	error('Only the last argument in getnc can be a structure')
      end
    else
      switch ii
       case 1
	file = substitute_new(cell_args{ii}, order_def, 1);
       case 2
	varid = substitute_new(cell_args{ii}, order_def, 1);
       case 3
	bl_corner = substitute_new(cell_args{ii}, order_def, 1);
       case 4
	tr_corner = substitute_new(cell_args{ii}, order_def, 1);
       case 5
	stride = substitute_new(cell_args{ii}, order_def, 1);
       case 6
	order = substitute_new(cell_args{ii}, order_def, 1);
       case 7
	change_miss = substitute_new(cell_args{ii}, change_miss_def, 1);
       case 8
	new_miss = substitute_new(cell_args{ii}, new_miss_def, 0);
       case 9
	squeeze_it = substitute_new(cell_args{ii}, squeeze_it_def, 1);
       case 10
	rescale_opts = substitute_new(cell_args{ii}, rescale_opts_def, 1);
       case 11
	err_opt = substitute_new(cell_args{ii}, err_opt_def, 1);
       otherwise
	error('There can only be 11 arguments passed to getnc')
      end	
    end
  end
  
  % Flag errors if some necessary argument has not been passed or now
  % contains an unacceptable value.
  
  % 1: file
  
  if ~ischar(file)
    error(' FILE is not a string');
  end

  % 2: varid
  
  if isnumeric(varid)
    size_var = size(varid);
    if prod(size(varid)) == 1
      if varid < 0
	error('ERROR: varid is less than zero');
      end
    else
      error('varid must be a scalar or a string')
    end
  elseif ~ischar(varid)
    error('varid must be a scalar or a string')
  end
  
  % 7: change_miss
  
  is_bad = 1;
  if isnumeric(change_miss)
    if length(change_miss) == 1
      if ((change_miss >= 1) & (change_miss <= 3)) | (change_miss <= -1)
	is_bad = 0;
      end
    end
  end
  if is_bad
    error('change_miss must be a scalar')
  end
    
  % 8: new_miss
  
  is_bad = 1;
  if isnumeric(new_miss)
    if length(new_miss) == 1
	is_bad = 0;
    end
  end
  if is_bad
    error('new_miss must be a scalar')
  end
    
  % 9: squeeze_it
  
  is_bad = 1;
  if isnumeric(squeeze_it)
    if length(squeeze_it) == 1
	is_bad = 0;
    end
  end
  if is_bad
    error('squeeze_it must be a scalar')
  end
    
  % 11: err_opt
  
  is_bad = 1;
  if isnumeric(err_opt)
    if length(err_opt) == 1
      if (err_opt >= -1) & (err_opt <= 2)
	is_bad = 0;
      end
    end
  end
  if is_bad
    error('err_opt must be a scalar == 0, 1 or 2')
  end
  % disp(['a:bl_corner = ' num2str(bl_corner)])
  % disp(['a:tr_corner = ' num2str(tr_corner)])
  % disp(['a:stride = ' num2str(stride)])

function return_val = substitute_new(new_val, default_val, check_neg1)
  
% If new_val is empty set return_val == default_val.
% If new_val == -1 and check_neg1 is .true. set return_val = default_val.
% Otherwise set return_val = new_val.
  
  if isempty(new_val)
    return_val = default_val;
  elseif (check_neg1 & all(new_val == -1))
    return_val = default_val;
  else
    return_val = new_val;
  end
