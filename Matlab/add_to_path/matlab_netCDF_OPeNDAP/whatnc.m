function whatnc

% WHATNC lists all of the netCDF files in the current directory
%--------------------------------------------------------------------
%     Copyright (C) J. V. Mansbridge, CSIRO, january 24 1992
%     Revision $Revision: 1.7 $
%
% DESCRIPTION:
% whatnc lists all of the netCDF files (including compressed ones) in
% the current directory.  It also lists all of the netcdf files in the
% common data set.
%
% Note: the path for the common data set is found by a call to pos_cds.
%
% INPUT:
% none
%
% OUTPUT:
% messages to the user's terminal
%
% EXAMPLE:
% Simply type whatnc at the matlab prompt.
%
% This function calls: pos_cds, message_for_whatnc
% This function is called by:
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

%     Copyright (C), J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: whatnc.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

% 1) current directory netCDF files

dir_name = '.';
dir_list = dir(dir_name);

string_list = {'.nc', '.cdf'};
num_strings = length(string_list);
for jj = 1:num_strings
  len_string(jj) = length(string_list{jj});
end
names = {};
count = 0;
len = zeros(length(dir_list), 1);
for ii = 1:length(dir_list)
  if ~dir_list(ii).isdir
    for jj = 1:num_strings
      if length(dir_list(ii).name) > len_string(jj)
	if strcmp(dir_list(ii).name((end - len_string(jj) + 1):end), ...
		  string_list{jj})
	  count = count + 1;
	  names{count} = dir_list(ii).name;
	  len(count) = length(names{count});
	end
      end
    end
  end
end
maxlen = max(len);

disp(' ')
disp('-----  current directory netCDF files  -----')
print_list(names, count, maxlen)

% 2) current directory compressed netCDF files

string_list = {'.nc.Z', '.cdf.Z', '.nc.gz', '.cdf.gz'};
num_strings = length(string_list);
for jj = 1:num_strings
  len_string(jj) = length(string_list{jj});
end
names = {};
count = 0;
len = zeros(length(dir_list), 1);
for ii = 1:length(dir_list)
  if ~dir_list(ii).isdir
    for jj = 1:num_strings
      if length(dir_list(ii).name) > len_string(jj)
	if strcmp(dir_list(ii).name((end - len_string(jj) + 1):end), ...
		  string_list{jj})
	  count = count + 1;
	  names{count} = dir_list(ii).name;
	  len(count) = length(names{count});
	end
      end
    end
  end
end
maxlen = max(len);

disp(' ')
disp('-----  current directory compressed netCDF files  -----')
print_list(names, count, maxlen)

% 3) common data set of netCDF files

dir_name = pos_cds;
if exist(dir_name, 'dir')
  dir_list = dir(dir_name);

  string_list = {'.nc', '.cdf'};
  num_strings = length(string_list);
  for jj = 1:num_strings
    len_string(jj) = length(string_list{jj});
  end
  names = {};
  count = 0;
  len = zeros(length(dir_list), 1);
  for ii = 1:length(dir_list)
    if ~dir_list(ii).isdir
      for jj = 1:num_strings
	if length(dir_list(ii).name) > len_string(jj)
	  if strcmp(dir_list(ii).name((end - len_string(jj) + 1):end), ...
		    string_list{jj})
	    count = count + 1;
	    names{count} = dir_list(ii).name;
	    len(count) = length(names{count});
	  end
	end
      end
    end
  end
  maxlen = max(len);

  disp(' ')
  disp('-----  common data set of netCDF files -----')
  print_list(names, count, maxlen)
end

% 4) look for a message to be added to which.

disp(' ')
disp('--------------------------------------------')
try
  mess = message_for_whatnc;
  if ischar(mess)
    disp(mess)
  elseif iscell(mess)
    for ii = 1:length(mess)
      if ischar(mess{ii})
	disp(mess{ii})
      end
    end
  end
end

function print_list(names, count, maxlen)
  if count > 0
    str_blank = repmat(' ', 1, 50);
    if maxlen < 10
      lim = 10;
      ncols = 7;
    elseif maxlen < 12
      lim = 12;
      ncols = 6;
    elseif maxlen < 15
      lim = 15;
      ncols = 5;
    elseif maxlen < 18
      lim = 18;
      ncols = 4;
    elseif maxlen < 25
      lim = 25;
      ncols = 3;
    elseif maxlen < 37
      lim = 37;
      ncols = 2;
    else
      lim = -1;
    end
    if lim > 0
      nrows = ceil(count/ncols);
      str = repmat(' ', nrows, ncols*lim);
      xx = zeros(ncols, nrows);
      xx(1:count) = NaN;
      xx = xx';
      ff = find(isnan(xx));
      [ii, jj] = ind2sub([nrows ncols], ff); 
      for kk = 1:count
	jst = (jj(kk) - 1)*lim + 1;
	jfi = jst + length(names{kk}) - 1;
	str(ii(kk), jst:jfi) = names{kk};
      end
      disp(str)
    else
      for ii = 1:count
	disp(names{ii})
      end
    end
  else
    disp('EMPTY')
  end
function print_list_opposite(names, count, maxlen)
  if count > 0
    str_blank = repmat(' ', 1, 50);
    if maxlen < 10
      lim = 10;
      str = [];
      for ii = 1:count
	str = [str names{ii} str_blank(1:(lim - length(names{ii})))];
	if rem(ii, 7) == 0
	  str = [str '\n'];
	end
      end
      disp(str)
    else
      for ii = 1:count
	disp(names{ii})
      end
    end
  else
    disp('EMPTY')
  end

