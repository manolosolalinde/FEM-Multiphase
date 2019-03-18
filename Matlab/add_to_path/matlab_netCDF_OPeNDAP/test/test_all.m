% test_all tests the CSIRO matlab/netcdf interface. The user is prompted to
% test its ability to read both netcdf and opendap files. The retrieved
% values are compared with the "correct" ones stored in a mat file.

% $Id: test_all.m,v 1.1 2005/10/27 05:10:02 man133 Exp man133 $
% Copyright J. V. Mansbridge, CSIRO, Mon Aug 29 11:53:06 EST 2005

clear all

var = {'time', 'u', 'ureverse', 'uchar1', 'uchar2', 'uchar3', 'ushort', ...
       'ulong', 'udouble', 'no_atts', 'big_var', 'depth1', 'depth2', 'dim3'};
order_array = {1, [2 1], [1 2], -1, [1 2], [2 3 1], [1 3 2], ...
	       [3 2 1 4], [2 1], -2, [2 3 5 1 4], -1, -2, [1]};
index_pos_char = [23 37 51 65 79 93];
index_list_char = [index_pos_char (index_pos_char + 1) (index_pos_char + 2)];

load test_1.mat

% Read a netcdf file.

while 1
  reply = input('Do you want to test the netcdf installation? (y/n): ', 's');
  reply = lower(reply);
  if strcmp(reply, 'y')
    if ~exist('mexnc')
      error('!! Didn''t find mexnc (from the netcdf interface)')
    end
    file = './test_1.nc';
    test_type = 'netcdf';
    t1 = now;
    test_general_all
    t2 = now;
    disp(' ')
    disp(['netcdf test took ' num2str((t2 - t1)*86400) ' seconds'])
    break
  elseif strcmp(reply, 'n')
    break
  else
    disp('Answer must be y or n')
  end
end

% Read an opendap file.

while 1
  reply = input('Do you want to test the opendap installation? (y/n): ', 's');
  reply = lower(reply);
  if strcmp(reply, 'y')
    if ~exist('mexnc')
      error('!! Didn''t find mexnc (from the netcdf interface)')
    end
    file = 'http://www.marine.csiro.au/dods/nph-dods/dods-data/test_data/test_1.nc';
    test_type = 'opendap';
    t1 = now;
    test_general_all
    t2 = now;
    disp(' ')
    disp(['opendap test took ' num2str((t2 - t1)*86400) ' seconds'])
    break
  elseif strcmp(reply, 'n')
    break
  else
    disp('Answer must be y or n')
  end
end


