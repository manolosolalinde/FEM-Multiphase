function [gregorian_time, serial_time, gregorian_base, serial_base, ...
	  sizem, serial_time_jd, serial_base_jd] = ...
    timenc(file, time_var, bl_corner, tr_corner, calendar)
% TIMENC   Returns time information for a file that meets COARDS standards.
%
% function [gregorian_time, serial_time, gregorian_base, serial_base, ...
% 	  sizem, serial_time_jd, serial_base_jd] = ...
%         timenc(file, time_var, bl_corner, tr_corner, calendar);
%
% DESCRIPTION:
%
% timenc finds the time vector and the corresponding base date for a
% netcdf file or DODS/OPEnDAP dataset that meets COARDS standards. In
% practice this means that time-like variable should have a units attribute
% of a certain form. An example is:
%      'seconds since 1992-10-8 15:15:42.5 -6:00'
% This indicates seconds since October 8th, 1992 at 3 hours, 15
% minutes and 42.5 seconds in the afternoon in the time zone
% which is six hours to the west of Coordinated Universal Time
% (i.e. Mountain Daylight Time). The time zone specification can
% also be written without a colon using one or two-digits
% (indicating hours) or three or four digits (indicating hours
% and minutes).  Instead of 'seconds' the string may contain 'minutes',
% 'hours', 'days' and 'weeks' and all of these may be singular or plural
% and have capital first letters.  The letters 'UTC' or 'UT' are allowed at
% the end of the string, but these are ignored.
%
% The parsing of the units attribute is done by the function PARSETNC. This
% function can be called by itself to test that a string is sytactically
% correct.
%
% It is possible to have many different types of calendars but timenc only
% implements five at present. These are necessary because there is some
% confusion with dates before October 15 1582 when the Gregorian calendar was
% introduced. A problem also arises when the reference date in the units
% attribute is before this. timenc deals with this by recognising some of the
% CF conventions and returns different answers depending on the value of the
% calendar attribute of the time-like variable. Also, some numerical models
% like to pretend that every year has the same number of days - 365, 366 and
% 360 are all used.
%    1) calendar == 'standard', 'gregorian' (this is the default).
% In this case the relevant calculations are done using the functions
% get_calendar_date and get_julian_day which know about both the Julian and
% Gregorian calendars and so work back to julian day 0 in the year
% -4712. This means that the day after 4 October 1582 is 15 October 1582 as
% decreed by Pope Gregory XIII. This is the calendar almost universally used
% today and what udunits works with today. timenc has worked this way since
% revision 1.10 in 2000.
%    2) calendar == 'proleptic_gregorian'.
% In this case the relevant calculations are done using the matlab functions
% datenum and datevec which simply extend the way our present calendar works
% backwards into the past. This follows international standard ISO 8601. This
% means that dates are continuous, i.e., the day after 4 October 1582 is 5
% October 1582, which does NOT correspond to historical time anywhere. As well
% there is a year zero. timenc used to work this way before revision 1.10 in the
% year 2000 and I believe that udunits also did at some stage in the past.
%    3) calendar == 'noleap', '365_day'
% Here it is assumed that every year has 365 days.
%    4) calendar == 'all_leap', '366_day'
% Here it is assumed that every year has 366 days.
%    5) calendar == '360_day'
% Here it is assumed that every year has 360 days and every month has 30 days.
%
% Note that other values of the calendar attribute produce an error
% message. This can usually be avoided by the user specifying the calendar
% explicitly in the call to timenc.
%
% INPUT:
%  file may be the name of a netCDF file with or without the .cdf or .nc
%    extent. file may also be the URL of a DODS/OPEnDAP dataset.
% time_var: the name of the 'time' variable in the netCDF file or DODS/OPEnDAP
%    dataset.  If this argument is missing then it is assumed that variable
%    name is 'time'. If time_var is multi-dimensional then it will be handled
%    as if it had been reshaped as one 'giant' vector.
% bl_corner: a vector of length n specifying the hyperslab corner
%    with the lowest index values (the bottom left-hand corner in a
%    2-space).  The corners refer to the dimensions in the same
%    order that these dimensions are listed in the relevant questions
%    in getnc.m and in the inqnc.m description of the variable.  A
%    negative element means that all values in that direction will be
%    returned.  If this argument is missing or a negative scalar or empty
%    array is used this means that all of the elements in the array will be
%    returned.
%  tr_corner is a vector of length n specifying the hyperslab corner
%    with the highest index values (the top right-hand corner in a
%    2-space).  The corners refer to the dimensions in the same order
%    that these dimensions are listed in the relevant questions in
%    getnc.m and in the inqnc.m description of the variable. A negative
%    element means that the returned hyperslab should run to the highest
%    possible index. Note, however, that the value of an element in the
%    tr_corner vector will be ignored if the corresponding element in the
%    corner vector is negative.
%  calendar is a string determining the type of calendar to be used and is
%    discussed above.
%
% OUTPUT:
% gregorian_time: an Mx6 matrix where the rows refer to the M times
%    specified in the 'time' variable in the netCDF file.  The columns
%    are the year, month, day, hour, minute, second in that order, UT.
% serial_time: an M vector giving the serial times (in UT) specified in the
%    'time' variable in the netCDF file. Serial times are used by datestr,
%    datevec & datenum. Thus gregorian_time = datevec(serial_time). Note that
%    the 'time' variable actually contains the serial time relative to a
%    base time.
% gregorian_base: a 6-vector giving the year, month, day, hour, minute,
%    second of the base time as specified in the 'units' attribute of
%    the 'time' variable. This is in UT.
% serial_base: the serial time of the base time, in UT, as determined by
%    matlab's datenum function. Thus gregorian_base = datevec(serial_base).
%    serial_base will be a NaN for times before October 15 1582, when the
%    Gregorian calendar was adopted, since datenum is not meaningful in this
%    case.
% sizem: the size of the 'time' variable in the netCDF file.
% serial_time_jd: an M vector giving the julian day number (in UT) specified
%    in the 'time' variable in the netCDF file. (julian day numbers are used
%    by get_julian_day and get_calendar_date. Thus gregorian_time =
%    get_calendar_date(serial_time_jd).
% serial_base_jd: the Julian day number of the base time, in UT, as
%    determined by get_julian_day. Thus gregorian_base =
%    get_calendar_date(serial_base_jd).
%
% EXAMPLES:
%       1)
%   [gregorian_time, serial_time] = timenc('file.nc', 't');
%
% This gives you the time information in gregorian and 'matlab' serial
% time for a netcdf file named file.nc containing a time variable t.
%
%       2)
%   [gregorian_time, serial_time] = timenc('file.nc', 't', 3, 8);
%
% This gives you the time information in gregorian and 'matlab' serial
% time for elements 3 to 8 of the time variable t.
%
%       3)
%   [gregorian_time, serial_time] = timenc('file.nc', 't', -1, -1, 'gregorian');
%
% The use of the 5th (calendar) input allows you to override the calendar
% attribute in the original netcdf file. This attribute may be invalid or
% simply not supported by timenc.
%
%     Copyright J. V. Mansbridge, CSIRO, Tue May  9 11:36:06 EST 1995
%     (with the ability to handle a multi-dimensional time variable
%     added by Rose O'Connor).

% This function calls: check_nc.m, get_calendar_date.m, loaddap or
%                      loaddods, mexnc, parsetnc.m, pos_cds.m

%$Id: timenc.m Mon, Wed Nov 29 16:35:09 EDT 2006 $

% Process the input arguments

if ~ismember(nargin , [1 2 4 5])
  error('timenc takes either 1, 2, 4 or 5 input arguments (try "help timenc")')
end

if nargin == 1
  time_var = 'time';
end

if nargin <= 2
  get_all = 1;
elseif nargin >= 4
  if isempty(bl_corner) | isempty(tr_corner)
    get_all = 1;
  else
    if (bl_corner < 1)
      get_all = 1;
    else
      get_all = 0;
    end
  end
end

if nargin < 5
  calendar = [];
end

% Do some initialisation.

CSIRO_add_jar_file_maybe;
[mex_name, full_name, desc_das, file_status, exe_name] = ...
    choose_mexnc_opendap(file);

% We also use method_of_call == 2 which, according to the documentation, does
% the same thing as method_of_call == 1 but would be expected to be
% slower. However, this is necessary for loaddods to work properly on Windows
% boxes when an array is filled with characters. I have no idea what the bug
% might be in loaddods. There is also a hint that the same problem can occur
% with loaddap on some linux boxes.
  
method_of_call = 2;

switch mex_name
 case 'mexnc'
  
  %Open the netcdf file.

  [cdfid, rcode] = mexnc('OPEN', full_name, 'NC_NOWRITE');
  if rcode == -1
    error(['** ERROR ** ncopen: rcode = ' num2str(rcode)])
  end

  %Suppress all error messages from netCDF 

  [rcode] = mexnc('setopts', 0);

  %Get the id number of the variable 'time' and find out the info. about
  %the variable 'time'.

  [varid, rcode] = mexnc('varid', cdfid, time_var);
  if rcode == -1
    error(['** ERROR ** ncvarid: time variable = ''' time_var ''' not found'])
  end

  if get_all == 1
    [varnam, vartypv, nvdims, vdims, nvatts, rcode] = ...
	mexnc('ncvarinq', cdfid, varid);
    if rcode == -1
      error(['** ERROR ** ncvarinq: rcode = ' num2str(rcode) ...
	     ', time variable = ''' time_var ''''])
    end
    
    %Find out the size of the dimension 'time' which is the same as the
    %variable 'time'.

    %[name, sizem, rcode] = mexnc('ncdiminq', cdfid, vdims(1));
    %ROC
    sizem = zeros(1,nvdims);
    st    = zeros(1,nvdims);
    for i = 1:nvdims
      [name, sizem(i), rcode] = mexnc('ncdiminq', cdfid, vdims(i));
      if rcode == -1
	error(['** ERROR ** ncdiminq: rcode = ' num2str(rcode)])
      end
    end
    
    %Retrieve the elements of the variable 'time'.  These are the serial day
    %relative to the base date.

    if any(sizem == 0)
      serial_rel = []; % Presumably time is an unlimited dimension of length 0.
      disp(['Warning: There are apparently no ''time'' records'])
    else
      %[serial_rel, rcode] = mexnc('ncvarget', cdfid, varid, 0, sizem);
      [serial_rel, rcode] = mexnc('ncvarget', cdfid, varid, st, sizem);
      if rcode == -1
	error(['** ERROR ** ncvarget: rcode = ' num2str(rcode)])
      end
    end
  else
    [varnam, vartypv, nvdims, vdims, nvatts, rcode] = ...
	mexnc('ncvarinq', cdfid, varid);
    if rcode == -1
      error(['** ERROR ** ncvarinq: rcode = ' num2str(rcode) ...
	     ', time variable = ''' time_var ''''])
    end
    
    sizem = zeros(1,nvdims);
    for i = 1:nvdims
      [name, sizem(i), rcode] = mexnc('ncdiminq', cdfid, vdims(i));
      if rcode == -1
	error(['** ERROR ** ncdiminq: rcode = ' num2str(rcode)])
      end
      if tr_corner(i) < 0
	tr_corner(i) = sizem(i);
      end
    end
    
    if any(sizem == 0)
      serial_rel = []; % Presumably time is an unlimited dimension of length 0.
      disp(['Warning: There are apparently no ''time'' records'])
    elseif any(bl_corner < 1)
      error('bl_corner values are too small')
    elseif any(tr_corner > sizem)
      error('tr_corner values are too large')
    else
      [serial_rel, rcode] = mexnc('ncvarget', cdfid, varid, ...
				  bl_corner-1, tr_corner-bl_corner+1);
      if rcode == -1
	error(['** ERROR ** ncvarget: rcode = ' num2str(rcode)])
      end
    end
  end

  %Get the string describing the base date.

  [base_str, rcode] = mexnc('attget', cdfid, varid, 'units');
  if rcode == -1
    error(['** ERROR ** ncattget: rcode = ' num2str(rcode)])
  end

  %Close the netcdf file.

  [rcode] = mexnc('ncclose', cdfid);
  if rcode == -1
    error(['** ERROR ** ncclose: rcode = ' num2str(rcode)])
  end
  
 case {'loaddap', 'loaddods'}

  % Dealing with a DODS file.
  
  % Check the value of time_var.
  
 if ischar(time_var)
    if ~isfield(desc_das, time_var)
      error([time_var ' is not a variable in ' full_name])
    end
  else
    error('For dods data varid must be a string')
  end
  
  % Get the structure of time_var and then use this to extract required
  % information about attributes. We also strip off extraneous quote marks.
  
  %desc_varid = desc_das.(time_var);
  desc_varid = getfield(desc_das, time_var);
  if isfield(desc_varid, 'units')
    %base_str = desc_varid.('units');
    base_str = getfield(desc_varid, 'units');
    if ischar(base_str)
      base_str = base_str(2:(end - 1));
    else
      error(['The units attribute of ' time_var ' is not a string'])
    end
  else
    error([time_var ' does not have a units attribute']);
  end
  if isfield(desc_varid, 'DODS_ML_Size')
    %sizem = desc_varid.('DODS_ML_Size');
    sizem = getfield(desc_varid, 'DODS_ML_Size');
  else
    %xx = desc_varid.(time_var);
    xx = getfield(desc_varid, time_var);
    %sizem = xx.('DODS_ML_Size');
    sizem = getfield(xx, 'DODS_ML_Size');
  end
  num_dim = length(sizem);

  % Get the time vector. Note that we allow for the vector of variable values
  % to be at different "depths" in the structure. It seems that a variable
  % that is also a dimension will be at values_struct.time but otherwise it
  % will be at values_struct.time.time.
  
  if get_all
    switch method_of_call
     case 1
      switch mex_name
       case 'loaddap'
	values_struct = loaddap('+v', [full_name '?' time_var]);
       case 'loaddods'
	values_struct = loaddods('+v', [full_name '?' time_var]);
      end
     case 2
      switch mex_name
       case 'loaddap'
	loaddap('+v', [full_name '?' time_var]);
       case 'loaddods'
	loaddods('+v', [full_name '?' time_var]);
      end
      eval(['values_struct = ' time_var ';'])
    end
  else
    str_h = [];
    for ii = 1:num_dim
      if bl_corner(ii) < 0
	bl_corner(ii) = 1;
      elseif bl_corner(ii) > sizem(ii)
	error(['bl_corner so big that it is outside the hyperslab'])
      end
      if tr_corner(ii) < 0
	tr_corner(ii) = sizem(ii);
      elseif tr_corner(ii) > sizem(ii)
	error(['tr_corner so big that it is outside the hyperslab'])
      elseif tr_corner(ii) < bl_corner(ii)
	error(['tr_corner is less than bl_corner'])
      end      
      str_h = [str_h '[' num2str(bl_corner(ii) - 1) ':' ...
	       num2str(tr_corner(ii) - 1) ']'];
    end
    
    switch method_of_call
     case 1
      switch mex_name
       case 'loaddap'
	values_struct = loaddap('+v', [full_name '?' time_var str_h]);
       case 'loaddods'
	values_struct = loaddods('+v', [full_name '?' time_var str_h]);
      end
     case 2
      switch mex_name
       case 'loaddap'
	loaddap('+v', [full_name '?' time_var str_h]);
       case 'loaddods'
	loaddods('+v', [full_name '?' time_var str_h]);
      end
      eval(['values_struct = ' time_var ';'])
    end
  end
  
  % Allow for the various weird ways that the data may be returned.
  
  if isstruct(values_struct)
    %xx = values_struct.(time_var);
    xx = getfield(values_struct, time_var);
    if isstruct(xx)
      %serial_rel = xx.(time_var);
      serial_rel = getfield(xx, time_var);
    else
      serial_rel = xx;
    end
  else
    serial_rel = values_struct;
  end

  % Allow for serial_rel to be multi-dimensional in which case we have to
  % work around the weird way that loaddap returns multi-dimensional
  % data. Then make it into a column vector for backwards compatibility with
  % the code that Rose put in for a multi-dimensional time variable in netcdf
  % files.
  
  if num_dim > 2
    vec_permute = [num_dim:-1:3 1 2];
    serial_rel = permute(serial_rel, vec_permute);
  end
  serial_rel = serial_rel(:)';
 case 'java'
  
  % Get the file object
  
  try
    ncdJ = ucar.nc2.dataset.NetcdfDataset.openDataset(file);
  catch
    ss = lasterror;
    mess_str = ['Failed trying to open file: ' file ss.message];      
    rcode = -1000000;
    [gregorian_time, serial_time, gregorian_base, serial_base, ...
     sizem, serial_time_jd, serial_base_jd] = error_handle('java', ...
						  mess_str, rcode, err_opt);
    return
  end
  
  % Get the variable object
  
  try
    varJ = ncdJ.findVariable(time_var);
  catch
    ss = lasterror;
    mess_str = ['Failed getting ' time_var ' in file: ' file ss.message];      
    rcode = -1000000;
    [gregorian_time, serial_time, gregorian_base, serial_base, ...
     sizem, serial_time_jd, serial_base_jd] = error_handle('java', ...
						  mess_str, rcode, err_opt);
    return
  end
  
  % Get the actual data as serial_rel.
  
  sizem = varJ.getShape();
  if any(sizem == 0)
    serial_rel = []; % Presumably time is an unlimited dimension of length 0.
    disp(['Warning: There are apparently no ''time'' records'])
  else
    if get_all == 1
      serial_rel = squeeze(copyToNDJavaArray(varJ.read()));
    else
      % Check that the bl_corner points are acceptable.
      
      varRank = varJ.getRank;
      for ii = 1:varRank
	if (bl_corner(ii) < 1) || (bl_corner(ii) > tr_corner(ii)) || ...
	      (tr_corner(ii) > sizem(ii))
	  mess_str = 'hyperslab is badly specified';
	  rcode = -1000000;
	  [gregorian_time, serial_time, gregorian_base, serial_base, ...
	   sizem, serial_time_jd, serial_base_jd] = error_handle('java', ...
						  mess_str, rcode, err_opt);
	  return
	end
      end
      
      % Construct the string that specifies the hyperslab and then get the data.
      readSpec = '';
      for ii = 1:varRank
	readSpec = [readSpec num2str(bl_corner(ii) - 1) ':' ...
		    num2str(tr_corner(ii) - 1) ':1,'];
      end
      readSpec = readSpec(1:(end - 1));
      serial_rel = squeeze(copyToNDJavaArray(varJ.read(readSpec)));
    end
  end
    
  % Get the string describing the base date.
  
  try
    att = varJ.findAttribute('units');
    base_str = char(att.getStringValue());
  catch
    ss = lasterror;
    mess_str = ['Failed getting units attribute in file: ' file ss.message];
    rcode = -1000000;
    [gregorian_time, serial_time, gregorian_base, serial_base, ...
     sizem, serial_time_jd, serial_base_jd] = error_handle('java', ...
						  mess_str, rcode, err_opt);
    return
  end  

  % Close the file object
  
  ncdJ.close();
 case 'none'
  error(['Couldn''t find a suitable mex-file for reading ' file])
end

% Find out what calendar we are using. If there is no calendar attribute then
% set it to 'gregorian';

if isempty(calendar)
  [att_val, att_name_list] = attnc(full_name, time_var);
  calendar = 'gregorian';
  for ii = 1:length(att_name_list)
    if strcmp(lower(att_name_list{ii}), 'calendar')
      calendar = att_val{ii};
      break;
    end
  end
end

% Parse the string containing the base date to get its constituents and
% then find its serial and gregorian dates. Also rescale the relative serial
% time vector to turn it into days since the base time.

[gregorian_base, rescale_serial_rel, serial_base_jd, serial_base] = ...
    parsetnc(base_str);
if rescale_serial_rel ~= 1
  serial_rel = rescale_serial_rel*serial_rel;
end

% Find the absolute serial date and resultant gregorian date of the time
% vector.

serial_time_jd = serial_rel + serial_base_jd;
if isempty(serial_time_jd)
  gregorian_time = [];
  serial_time = [];
else
  switch lower(calendar)
   case {'standard', 'gregorian'}
    gregorian_time = get_calendar_date(serial_time_jd);
    serial_time = datenum(gregorian_time(:, 1), gregorian_time(:, 2), ...
			  gregorian_time(:, 3), gregorian_time(:, 4), ...
			  gregorian_time(:, 5), gregorian_time(:, 6));

   case 'proleptic_gregorian'
    serial_time = serial_rel + serial_base;
    serial_time = serial_time(:);
    gregorian_time = datevec(serial_time);
   case {'noleap', '365_day'}
    % We use serial_base to give us a proper starting time and work from
    % there in steps of 365 days per year.
    days_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];
    days_ref = [0 cumsum(days_per_month)];
    [year_b, month_b, day_b, hour_b, minute_b, sec_b] = datevec(serial_base);
    days_from_year_base = days_ref(month_b) + day_b - 1 + hour_b/24 + ...
       minute_b/1440 + sec_b/86400;
    day_full = serial_rel + days_from_year_base;
    year_rel = floor(day_full/365);
    year_abs = year_b + year_rel;
    rem_1 = day_full - year_rel*365;
    month_abs = zeros(1, sizem);
    for ii = 1:sizem
       ff = find(days_ref <= rem_1(ii));
       month_abs(ii) = ff(end);
       rem_2(ii) = rem_1(ii) - days_ref(month_abs(ii));
    end
    day_rel = floor(rem_2);
    day_abs = day_rel + 1;
    rem_3 = (rem_2 - day_rel)*24;
    hour_abs = floor(rem_3);
    rem_4 = (rem_3 - hour_abs)*60;
    minute_abs = floor(rem_4);
    second_abs = (rem_4 - minute_abs)*60;
    gregorian_time = [year_abs(:) month_abs(:) day_abs(:) hour_abs(:) ...
       minute_abs(:) second_abs(:)];
    serial_time = datenum(gregorian_time);
   case {'all_leap', '366_day'}
    % We use serial_base to give us a proper starting time and work from
    % there in steps of 366 days per year.
    days_per_month = [31 29 31 30 31 30 31 31 30 31 30 31];
    days_ref = [0 cumsum(days_per_month)];
    [year_b, month_b, day_b, hour_b, minute_b, sec_b] = datevec(serial_base);
    days_from_year_base = days_ref(month_b) + day_b - 1 + hour_b/24 + ...
       minute_b/1440 + sec_b/86400;
    day_full = serial_rel + days_from_year_base;
    year_rel = floor(day_full/366);
    year_abs = year_b + year_rel;
    rem_1 = day_full - year_rel*366;
    month_abs = zeros(1, sizem);
    for ii = 1:sizem
       ff = find(days_ref <= rem_1(ii));
       month_abs(ii) = ff(end);
       rem_2(ii) = rem_1(ii) - days_ref(month_abs(ii));
    end
    day_rel = floor(rem_2);
    day_abs = day_rel + 1;
    rem_3 = (rem_2 - day_rel)*24;
    hour_abs = floor(rem_3);
    rem_4 = (rem_3 - hour_abs)*60;
    minute_abs = floor(rem_4);
    second_abs = (rem_4 - minute_abs)*60;
    gregorian_time = [year_abs(:) month_abs(:) day_abs(:) hour_abs(:) ...
       minute_abs(:) second_abs(:)];
    serial_time = datenum(gregorian_time);
   case '360_day'
    % We use serial_base to give us a proper starting time and work from
    % there in steps of 366 days per year.
    [year_b, month_b, day_b, hour_b, minute_b, sec_b] = datevec(serial_base);
    days_from_year_base = 30*(month_b - 1) + day_b - 1 + hour_b/24 + ...
       minute_b/1440 + sec_b/86400;
    day_full = serial_rel + days_from_year_base;
    year_rel = floor(day_full/360);
    year_abs = year_b + year_rel;
    rem_1 = day_full - year_rel*360;
    month_abs = floor(rem_1/30) + 1;
    rem_2 = rem_1 - (month_abs - 1)*30;
    day_rel = floor(rem_2);
    day_abs = day_rel + 1;
    rem_3 = (rem_2 - day_rel)*24;
    hour_abs = floor(rem_3);
    rem_4 = (rem_3 - hour_abs)*60;
    minute_abs = floor(rem_4);
    second_abs = (rem_4 - minute_abs)*60;
    gregorian_time = [year_abs(:) month_abs(:) day_abs(:) hour_abs(:) ...
       minute_abs(:) second_abs(:)];
    serial_time = datenum(gregorian_time);
   otherwise
    disp(['!! timenc cannot handle the calendar attribute **' calendar '**'])
    disp('!! which may have been found in the original netCDF file. The help')
    disp('!! message for timenc tells you how to specify a different calendar')
    error('strange calendar')
  end
end

function [gregorian_time, serial_time, gregorian_base, serial_base, ...
	  sizem, serial_time_jd, serial_base_jd] = error_handle(fid, ...
						  mess_str, rcode, err_opt)

% error_handle is called after a mexnc or java call has failed. It ensures
% that an open netcdf file is closed. The value of err_opt determines what
% else is done. For a mexnc call fid is cdfid, the handle to the open
% file. For a java call fid is the opened file object.
%    err_opt == 1 prints an error message and then aborts
%            == 2 prints a warning message and then returns an empty
%                 array. This is the default.
%            == 3 returns an empty array. This is a very dangerous option and
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
  
  if ~exist('gregorian_time')
    gregorian_time = [];
  end
  if ~exist('serial_time')
    serial_time = [];
  end
  if ~exist('gregorian_base')
    gregorian_base = [];
  end
  if ~exist('serial_base')
    serial_base = [];
  end
  if ~exist('sizem')
    sizem = [];
  end
  if ~exist('serial_time_jd')
    serial_time_jd = [];
  end
  if ~exist('serial_base_jd')
    serial_base_jd = [];
  end
  
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

