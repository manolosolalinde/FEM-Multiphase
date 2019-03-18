function serial_time = get_serial_time(relative_time, time_units_str)
% get_serial_time: Find the absolute time value given a relative time. The
% string time_units_str specifies how to go from relative to absolute time.
%
%     serial_time = get_serial_time(relative_time, time_units_str)
%
%        INPUT:
% relative_time: Elapsed time since the starting time specified in
%    time_units_str; this may be an array.
% time_units_str: The string stored as the units attribute of the time
%    variable in a COARDS standard netcdf file. It is like:
%    'seconds since 1992-10-8 15:15:42.5 -6:00'
%    OUTPUT:
% serial_time: An array giving serial times (in UT) as used by datestr,
%    datevec & datenum. Thus gregorian_time = datevec(serial_time). Note that
%    the 'time' variable in a COARDS standard netcdf file actually contains
%    the serial time relative to a base time.
%
% NOTES:
% In a COARDS standard netcdf file a time-like variable may be specified by
% the units attribute being a string like 
%      'seconds since 1992-10-8 15:15:42.5 -6:00'
% This indicates seconds since October 8th, 1992 at 3 hours, 15
% minutes and 42.5 seconds in the afternoon in the time zone
% which is six hours to the west of Coordinated Universal Time
% (i.e. Mountain Daylight Time). The time zone specification can
% also be written without a colon using one or two-digits
% (indicating hours) or three or four digits (indicating hours
% and minutes).  Instead of 'seconds' the string may contain 'minutes',
% 'hours', 'days' and 'weeks' and all of these may be singular or plural
% and have capital first letters.  I also allow the letters 'UTC' or
% 'UT' at the end of the string, but these are ignored.
%
% This function allows us to find the absolute time value given a relative
% time specified as above. Thus we could have the example below:
%  >> relative_time = [0 1 2];
%  >> time_units_str = 'weeks since 1990-10-8 13:00:00';
%  >> serial_time = get_serial_time(relative_time, time_units_str);
%  >> datestr(serial_time)
%  ans =
%  08-Oct-1990 13:00:00
%  15-Oct-1990 13:00:00
%  22-Oct-1990 13:00:00
%
% If the data is actually in a netcdf file then it is more convenient to use
% the function timenc which carries out the equivalent operation while doing
% extra checks and handling missing values as well.

% $Id: get_serial_time.m Mon, 03 Jul 2006 17:16:40 $
% Copyright J. V. Mansbridge, CSIRO, Tue Aug 12 15:30:49 EST 2003

[gregorian_base, rescale_serial_rel, serial_base_jd, serial_base] = ...
    parsetnc(time_units_str);
if rescale_serial_rel ~= 1
  relative_time = rescale_serial_rel*relative_time;
end
serial_time = relative_time + serial_base;
