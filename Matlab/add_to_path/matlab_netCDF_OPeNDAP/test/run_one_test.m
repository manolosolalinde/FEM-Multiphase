% This is called from test_general_all. Different values of the variable test_no
% determining exactly what reading test will be carried out.

% $Id: run_one_test.m,v 1.3 2007/08/06 10:08:10 man133 Exp man133 $
% Copyright J. V. Mansbridge, CSIRO, Fri Jun 16 17:43:42 EST 2006

if (test_no == 0)
  
  % Now test ddsnc

  try
    desc_dds_test = ddsnc(file);
  catch
    disp('A call to ddsnc failed with the following message:')
    disp(' ')
    disp(lasterr)
    disp(' ')
    error(['Either the file ' file ' is inaccessible or the interface is ' ...
	   'not working at all.'])
  end
  
  % Compare the layout of the two structures and check that all of the
  % expected fields are there.

  if ~isstruct(desc_dds_ref)
    error('desc_dds_ref is not a structure')
  end
  if ~isstruct(desc_dds_test)
    error('desc_dds_test is not a structure')
  end

  if ~isfield(desc_dds_test, 'variable')
    error('desc_dds_test does not have a field named "variable"')
  end
  if ~isfield(desc_dds_test.variable, 'type')
    error('desc_dds_test.variable does not have a field named "type"')
  end
  if ~isfield(desc_dds_test.variable, 'name')
    error('desc_dds_test.variable does not have a field named "name"')
  end
  if ~isfield(desc_dds_test.variable, 'dim_statement')
    error('desc_dds_test.variable does not have a field named "dim_statement"')
  end
  if ~isfield(desc_dds_test.variable, 'dim_idents')
    error('desc_dds_test.variable does not have a field named "dim_idents"')
  end

  if ~isfield(desc_dds_test, 'dimension')
    error('desc_dds_test does not have a field named "dimension"')
  end
  if ~isfield(desc_dds_test.dimension, 'name')
    error('desc_dds_test.dimension does not have a field named "name"')
  end
  if ~isfield(desc_dds_test.dimension, 'length')
    error('desc_dds_test.dimension does not have a field named "length"')
  end

  % Check that the dimensions have the same names and lengths.

  num_elements = length(desc_dds_ref.dimension);
  for ii = 1:num_elements
    val = desc_dds_ref.dimension(ii).name;
    jj = look_for_val(desc_dds_test, 'dimension', 'name', val, 1);
    if desc_dds_ref.dimension(ii).length ~= desc_dds_test.dimension(jj).length
      error(['length of ref dimension ' val ' == ' ...
	     num2str(desc_dds_ref.dimension(ii).length) ' but is ' ...
	     num2str(desc_dds_test.dimension(jj).length) ' for test case'])
    end
  end

  % Check that the variables have the same names and dimensions. The types
  % must be mapped to common names since the call to ddsnc can give different
  % names to the types according to whether mexnc or java calls are made.

  num_elements = length(desc_dds_ref.variable);
  for ii = 1:num_elements
    val = desc_dds_ref.variable(ii).name;
    jj = look_for_val(desc_dds_test, 'variable', 'name', val, 1);
    include_type_test = 0;
    if include_type_test
      type_ref = lower(desc_dds_ref.variable(ii).type);
      switch type_ref
       case 'char'
	type_ref = 'string';
       case 'short'
	type_ref = 'int16';
       case 'int'
	type_ref = 'int32';
       case 'long'
	type_ref = 'int64';
       case 'float'
	type_ref = 'float32';
       case 'double'
	type_ref = 'float64';
      end
      type_test = lower(desc_dds_test.variable(jj).type);    
      switch type_test
       case 'char'
	type_test = 'string';
       case 'short'
	type_test = 'int16';
       case 'int'
	type_test = 'int32';
       case 'long'
	type_test = 'int64';
       case 'float'
	type_test = 'float32';
       case 'double'
	type_test = 'float64';
      end
      if ~strcmp(type_ref, type_test)
	error(['type of ref variable ' val ' is ' ...
	       desc_dds_ref.variable(ii).type ' but is ' ...
	       desc_dds_test.variable(jj).type ' for test case'])
      end
    end
    if ~strcmp(desc_dds_ref.variable(ii).dim_statement, ...
	       desc_dds_test.variable(jj).dim_statement)
      error(['dim_statement of ref variable ' val ' is ' ...
	     desc_dds_ref.variable(ii).dim_statement ' but is ' ...
	     desc_dds_test.variable(jj).dim_statement ' for test case'])
    end
  end
  
elseif (test_no >= 1) & (test_no <= 4)

  % First test timenc
  
  k = test_no;
  switch k
   case 1
    [gt, st, gb, sb, si, stj, sbj] = timenc(file);
    ra = 1:3;
   case 2
    [gt, st, gb, sb, si, stj, sbj] = timenc(file, 'time');
    ra = 1:3;
   case 3
    [gt, st, gb, sb, si, stj, sbj] = timenc(file, 'time', -1, -1);
    ra = 1:3;
   case 4
    [gt, st, gb, sb, si, stj, sbj] = timenc(file, 'time', 2, 3);
    ra = 2:3;
  end

  if size(gt) ~= size(gregorian_time(ra, :))
    error('!! timenc failed: produces wrong size of gregorian_time')
  else
    if size(gt) > 0
      dd = abs(gt - gregorian_time(ra, :)) + ...
	   abs(isnan(gt) - isnan(gregorian_time(ra, :)));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of gregorian_time')
      end
    end
  end

  if size(st) ~= size(serial_time(ra))
    error('!! timenc failed: produces wrong size of serial_time')
  else
    if size(st) > 0
      dd = abs(st - serial_time(ra)) + abs(isnan(st) - ...
					    isnan(serial_time(ra)));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of serial_time')
      end
    end
  end

  if size(gb) ~= size(gregorian_base)
    error('!! timenc failed: produces wrong size of gregorian_base')
  else
    if size(gb) > 0
      dd = abs(gb - gregorian_base) + abs(isnan(gb) - isnan(gregorian_base));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of gregorian_base')
      end
    end
  end

  if size(sb) ~= size(serial_base)
    error('!! timenc failed: produces wrong size of serial_base')
  else
    if size(sb) > 0
      dd = abs(sb - serial_base) + abs(isnan(sb) - isnan(serial_base));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of serial_base')
      end
    end
  end

  if size(si) ~= size(sizem)
    error('!! timenc failed: produces wrong size of sizem')
  else
    if size(si) > 0
      dd = abs(si - sizem) + abs(isnan(si) - isnan(sizem));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of sizem')
      end
    end
  end

  stj = stj(:);
  serial_time_jd_ra = serial_time_jd(ra);
  serial_time_jd_ra = serial_time_jd_ra(:);
  if any(size(stj) ~= size(serial_time_jd_ra))
    error('!! timenc failed: produces wrong size of serial_time_jd')
  else
    if size(stj) > 0
      dd = abs(stj - serial_time_jd_ra) + ...
	   abs(isnan(stj) - isnan(serial_time_jd_ra));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of serial_time_jd')
      end
    end
  end

  if size(sbj) ~= size(serial_base_jd)
    error('!! timenc failed: produces wrong size of serial_base_jd')
  else
    if size(sbj) > 0
      dd = abs(sbj - serial_base_jd) + abs(isnan(sbj) - isnan(serial_base_jd));
      if max(dd(:)) > 0
	error('!! timenc failed: produces wrong value of serial_base_jd')
      end
    end
  end
 
elseif (test_no == 5)
  
  % Now test the global attributes

  [avg, anlg] = attnc(file, 'global');
  if size(avg) ~= size(att_val_global)
    error('!! attnc failed: produces wrong number of global attributes')
  else
    for ii = 1:length(avg)
      if ischar(avg{ii}) & ischar(att_val_global{ii})
	if ~strcmp(avg{ii}, att_val_global{ii})
	  error(['!! attnc failed: 1: global attribute ' num2str(ii) ' wrong'])
	end
      elseif isnumeric(avg{ii}) & isnumeric(att_val_global{ii})
	if size(avg{ii}) == size(att_val_global{ii})
	  dd = abs(avg{ii} - att_val_global{ii}) + ...
	       abs(isnan(avg{ii}) - isnan(att_val_global{ii}))
	  if max(dd(:)) > 0
	    error(['!! attnc failed: 2: global attribute ' num2str(ii) ' wrong'])
	  end
	else
	  error(['!! attnc failed: 3: global attribute ' num2str(ii) ' wrong'])
	end
      else
	error(['!! attnc failed: 4: global attribute ' num2str(ii) ' wrong'])
      end
    end
  end
  
elseif (test_no >= 6) & (test_no <= 19)

  % Now test the attributes of each of the variables using a 2 argument call
  % to attnc inside a single loop.

  ii = test_no - 5;
  [av, anl] = attnc(file, var{ii});
  cmd = ['att_val = att_val_' var{ii} ';'];
  eval(cmd)
  cmd = ['att_name_list = att_name_list_' var{ii} ';'];
  eval(cmd)
  if length(av) < length(att_val)
    % Note that there may be extra attributes returned by attnc because the
    % opendap setver could add extra information.
    error(['!! attnc failed: finds wrong number of attributes for ' var{ii}])
  else
    for jj = 1:length(att_name_list)
      found_match = 0;
      for kk = 1:length(av)
	if strcmp(att_name_list{jj}, anl{kk})
	  found_match = 1;
	  break
	end
      end
      if ~found_match
	error(['!! attnc failed: 0.5: could not find attribute ' ...
	       att_name_list{jj}])
      end
      if ischar(av{kk}) & ischar(att_val{jj})
	if ~strcmp(av{kk}, att_val{jj})
	  error(['!! attnc failed: 1: attribute ' num2str(jj) ' of ' ...
		 var{ii} ' wrong'])
	end
      elseif isnumeric(av{kk}) & isnumeric(att_val{jj})
	if size(av{kk}) == size(att_val{jj})
	  switch test_type
	    case 'netcdf'
	     dd = abs(av{kk} - att_val{jj}) + ...
		  abs(isnan(av{kk}) - isnan(att_val{jj}));
	     if max(dd(:)) > 0
	       error(['!! attnc failed: 2: attribute ' num2str(jj) ' of ' ...
		      var{ii} ' wrong'])
	     end
	   case 'opendap' % allow for possible rounding error.
	     dd1 = abs(isnan(av{kk}) - isnan(att_val{jj}));
	     if max(dd1(:)) > 0
	       error(['!! attnc failed: 2a: attribute ' num2str(jj) ' of ' ...
		      var{ii} ' wrong (NaNs don''t match'])
	     end
	     % Check for a large relative error. This is necessary because of
             % loss of precision in when using opendap.
	     dd2 = abs(av{kk} - att_val{jj});
	     if max(dd2(:)) > 0
	       vv = att_val{jj};
	       ff = find(abs(vv) > 0);
	       dd2(ff) = dd2(ff)./vv(ff);
	       if max(dd2(:)) > 1e-5
		 error(['!! attnc failed: 2b: attribute ' num2str(jj) ' of ' ...
			var{ii} ' wrong'])
	       end
	     end	    
	  end
	else
	  error(['!! attnc failed: 3: attribute ' num2str(jj) ' of ' ...
		 var{ii} ' wrong'])
	end
      else
	error(['!! attnc failed: 4: attribute ' num2str(jj) ' of ' ...
		 var{ii} ' wrong'])
      end
    end
  end

elseif (test_no >= 20) & (test_no <= 33)

  % Now test the attributes of each of the variables using a 3 argument call
  % to attnc inside a double loop.

  ii = test_no - 19;
  cmd = ['att_val = att_val_' var{ii} ';'];
  eval(cmd)
  cmd = ['att_name_list = att_name_list_' var{ii} ';'];
  eval(cmd)
  for jj = 1:length(att_name_list)
    [av, anl] = attnc(file, var{ii}, att_name_list{jj});

    if ischar(av) & ischar(att_val{jj})
      if ~strcmp(av, att_val{jj})
	error(['!! attnc failed: 1: attribute ' num2str(jj) ' of ' ...
	       var{ii} ' wrong'])
      end
    elseif isnumeric(av) & isnumeric(att_val{jj})
      if size(av) == size(att_val{jj})
	switch test_type
	 case 'netcdf'
	  dd = abs(av - att_val{jj}) + ...
	       abs(isnan(av) - isnan(att_val{jj}));
	  if max(dd(:)) > 0
	    error(['!! attnc failed: 2: attribute ' num2str(jj) ' of ' ...
		   var{ii} ' wrong'])
	  end
	 case 'opendap' % allow for possible rounding error.
	  dd1 = abs(isnan(av) - isnan(att_val{jj}));
	  if max(dd1(:)) > 0
	    error(['!! attnc failed: 2a: attribute ' num2str(jj) ' of ' ...
		   var{ii} ' wrong (NaNs don''t match'])
	  end
	  % Check for a large relative error. This is necessary because of
	  % loss of precision in when using opendap.
	  dd2 = abs(av - att_val{jj});
	  if max(dd2(:)) > 0
	    vv = att_val{jj};
	    ff = find(abs(vv) > 0);
	    dd2(ff) = dd2(ff)./vv(ff);
	    if max(dd2(:)) > 1e-5
	      error(['!! attnc failed: 2b: attribute ' num2str(jj) ' of ' ...
		     var{ii} ' wrong'])
	    end
	  end	    
	end
      else
	error(['!! attnc failed: 3: attribute ' num2str(jj) ' of ' ...
	       var{ii} ' wrong'])
      end
    else
      error(['!! attnc failed: 4: attribute ' num2str(jj) ' of ' ...
	     var{ii} ' wrong'])
    end
  end

elseif (test_no >= 34) & (test_no <= 47)

  % Now step through each of the variables. When the arrays are found replace
  % appropriate values by NaNs. (Values which match the missing_value or
  % _FillValue or lie outside the valid_range, valid_min or valid_max.)

  ii = test_no - 33;
  cmd = ['val = ' var{ii} '_nan;'];
  eval(cmd)
  
  % Get the entire array, with missing values replaced by NaNs, and test it. The
  % array is retrieved in different ways including sending a structure as an
  % argument to getnc.
  
  switch mod(ii, 3)
   case 0
    xx = getnc(file, var{ii});
   case 1
    clear struc
    struc.varid = var{ii};
    xx = getnc(file, struc);
   case 2
    clear struc
    struc.file = file;
    struc.varid = var{ii};
    xx = getnc(struc);
  end
  compare_mats_getnc(val, xx, var{ii}, 'nan')

elseif (test_no >= 48) & (test_no <= 61)

  % Now step through each of the variables. When the arrays are found do not
  % replace any values by NaNs.

  ii = test_no - 47;
  cmd = ['val = ' var{ii} '_nonan;'];
  eval(cmd)
  
  % Get the entire array, without any NaNs, and test it. The array is retrieved
  % in different ways including sending a structure as an argument to getnc.
  
  switch mod(ii, 4)
   case 0
    xx = getnc(file, var{ii}, -1, -1, -1, -1, 1);
   case 1
    clear struc
    struc.varid = var{ii};
    struc.change_miss = 1;
    xx = getnc(file, struc);
   case 2
    clear struc
    struc.change_miss = 1;
    xx = getnc(file, var{ii}, struc);
   case 3
    clear struc
    struc.file = file;
    struc.varid = var{ii};
    struc.change_miss = 1;
   xx = getnc(struc);
  end
  compare_mats_getnc(val, xx, var{ii}, 'nonan')

elseif (test_no >= 62) & (test_no <= 75)

  % Now step through each of the variables. When the arrays are found replace
  % appropriate values by 64. (Values which match the missing_value or
  % _FillValue or lie outside the valid_range, valid_min or valid_max.)

  ii = test_no - 61;
   cmd = ['val = ' var{ii} '_mi;'];
  eval(cmd)
  
  % Get the entire array, with missing values replaced by 64 and test it.
  
  switch mod(ii, 4)
   case 0
    xx = getnc(file, var{ii}, -1, -1, -1, -1, 3, 64);
   case 1
    clear struc
    struc.varid = var{ii};
    struc.new_miss = 64;
    xx = getnc(file, var{ii}, -1, -1, -1, -1, 3, struc);
   case 2
    clear struc
    struc.varid = var{ii};
    struc.change_miss = 3;
    struc.new_miss = 64;
    xx = getnc(file, var{ii}, -1, -1, -1, -1, 3, struc);
   case 3
    clear struc
    struc.file = file;
    struc.varid = var{ii};
    struc.change_miss = 3;
    struc.new_miss = 64;
    xx = getnc(struc);
  end
  compare_mats_getnc(val, xx, var{ii}, 'mi')

elseif (test_no >= 76) & (test_no <= 89)  

  % Now get a subset of the data.

  ii = test_no - 75;
  cmd = ['val = ' var{ii} '_nonan;'];
  eval(cmd)
  
  % Specify a subset of the array. Note the need to fiddle the case of array
  % being a vector. Note also that we slip in -1 a few times just to show
  % that it can be done.
  
  si_full = size(val);
  ref_array = zeros(si_full);
  if ndims(val) == 2
    if si_full(1) == 1
      si_full = si_full(2);
    elseif si_full(2) == 1
      si_full = si_full(1);
    end
  end
  corner = max(floor(si_full.*0.3), 1);
  end_point = ceil(si_full.*0.5);
  stride = ceil(si_full.*0.2);
  switch length(si_full)
   case 1
    val_sub = val(corner(1):stride(1):end_point(1));
   case 2
    val_sub = val(corner(1):stride(1):end_point(1), ...
		  corner(2):stride(2):end_point(2));
   case 3
    val_sub = val(corner(1):stride(1):end_point(1), ...
		  corner(2):end_point(2), ...
		  corner(3):stride(3):end_point(3));
    stride(2) = -1;
   case 4
    val_sub = val(corner(1):stride(1):end_point(1), ...
		  corner(2):stride(2):end, ...
		  corner(3):stride(3):end_point(3), ...
		  corner(4):stride(4):end_point(4));
    end_point(2) = -1;
   case 5
    val_sub = val(1:1:end, ...
		  corner(2):stride(2):end_point(2), ...
		  corner(3):stride(3):end_point(3), ...
		  corner(4):stride(4):end_point(4), ...
		  corner(5):stride(5):end_point(5));
    corner(1) = -1;
    end_point(1) = -1;
  end
  
  % Now get the subset of the array, without NaNs, and test it.

  switch mod(ii, 3)
   case 0
    xx = getnc(file, var{ii}, corner, end_point, stride, -1, 1);
   case 1
    clear struc
    struc.stride = stride;
    struc.change_miss = 1;
    xx = getnc(file, var{ii}, corner, end_point, struc);
   case 2
    clear struc
    struc.file = file;
    struc.varid = var{ii};
    struc.bl_corner = corner;
    struc.tr_corner = end_point;
    struc.stride = stride;
    struc.change_miss = 1;
    xx = getnc(struc);
  end
  compare_mats_getnc(val_sub, xx, var{ii}, 'sub')
  
elseif (test_no >= 90) & (test_no <= 103)  

  % Now step through each of the variables. When the arrays are found do not
  % replace any values by NaNs. Permute the arrays using the order_array
  % cell. For opendap access only test those occasions for which the order is a
  % scalar since this is the only version supported.

  ii = test_no - 89;
  order = order_array{ii};
  if strcmp(test_type, 'netcdf') | ...
	(strcmp(test_type, 'opendap') & (length(order) == 1))
    % disp(['Doing ' var{ii}])
    cmd = ['val = ' var{ii} '_nonan;'];
    eval(cmd)
    if length(order) == 1
      if (order == -2)
	si = size(val);
	nd = length(si);
	if (nd > 2) | ((si(1) > 1) & (si(2) > 1))
	  val = permute(val, length(size(val)):-1:1);
	end
      end
    else
      val = permute(val, order);
    end
    
    % First get the entire array, without any NaNs, and test it.
    
  switch mod(ii, 3)
   case 0
    xx = getnc(file, var{ii}, -1, -1, -1, order, 1);
   case 1
    clear struc
    struc.stride = stride;
    struc.change_miss = 1;
    xx = getnc(file, var{ii}, -1, -1, [], order, struc);
   case 2
    clear struc
    struc.file = file;
    struc.varid = var{ii};
    struc.order = order;
    struc.change_miss = 1;
    xx = getnc(struc);
  end
    compare_mats_getnc(val, xx, var{ii}, 'order')
  end

elseif (test_no >= 104) & (test_no <= 117)  

  % Now get a subset of the data. This is the same as the 76 to 89 case
  % except that we use different factors for the corner, end_point and
  % stride.

  ii = test_no - 103;
  cmd = ['val = ' var{ii} '_nonan;'];
  eval(cmd)
  
  % Specify a subset of the array. Note the need to fiddle the case of array
  % being a vector. Note also that we slip in -1 a few times just to show
  % that it can be done.
  
  si_full = size(val);
  ref_array = zeros(si_full);
  if ndims(val) == 2
    if si_full(1) == 1
      si_full = si_full(2);
    elseif si_full(2) == 1
      si_full = si_full(1);
    end
  end
  corner = max(floor(si_full.*0.3), 1);
  end_point = ceil(si_full.*0.9);
  stride = ceil(si_full.*0.3);
  switch length(si_full)
   case 1
    val_sub = val(corner(1):stride(1):end_point(1));
   case 2
    val_sub = val(corner(1):stride(1):end_point(1), ...
		  corner(2):stride(2):end_point(2));
   case 3
    val_sub = val(corner(1):stride(1):end_point(1), ...
		  corner(2):end_point(2), ...
		  corner(3):stride(3):end_point(3));
    stride(2) = -1;
   case 4
    val_sub = val(corner(1):stride(1):end_point(1), ...
		  corner(2):stride(2):end, ...
		  corner(3):stride(3):end_point(3), ...
		  corner(4):stride(4):end_point(4));
    end_point(2) = -1;
   case 5
    val_sub = val(1:1:end, ...
		  corner(2):stride(2):end_point(2), ...
		  corner(3):stride(3):end_point(3), ...
		  corner(4):stride(4):end_point(4), ...
		  corner(5):stride(5):end_point(5));
    corner(1) = -1;
    end_point(1) = -1;
  end
  
  % Allow for the possibility that val_sub is a row vector when our standard
  % specifies that only column vectors are to be returned.
  
  si_val_sub = size(val_sub);
  if length(si_val_sub) == 2
    if si_val_sub(1) == 1
      val_sub = val_sub(:);
    end
  end
  
  % Now get the subset of the array, without NaNs, and test it.

  switch mod(ii, 3)
   case 0
    xx = getnc(file, var{ii}, corner, end_point, stride, -1, 1);
   case 1
    clear struc
    struc.stride = stride;
    struc.change_miss = 1;
    xx = getnc(file, var{ii}, corner, end_point, struc);
   case 2
    clear struc
    struc.file = file;
    struc.varid = var{ii};
    struc.bl_corner = corner;
    struc.tr_corner = end_point;
    struc.stride = stride;
    struc.change_miss = 1;
    xx = getnc(struc);
  end
  compare_mats_getnc(val_sub, xx, var{ii}, 'sub')
  
end % end of giant if statement
