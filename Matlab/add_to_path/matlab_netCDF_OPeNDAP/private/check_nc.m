function status = check_nc(cdf)

%  CHECK_NC checks whether a netcdf or opendap file is accessible.
%--------------------------------------------------------------------
%     Copyright J. V. Mansbridge, CSIRO, january 24 1992
%     Revision $Revision: 1.4 $
% CHANGE
%
% status = check_nc(cdf)
%
% DESCRIPTION:
%  This checks whether the file cdf is accessible.  It returns a
%  status flag (see below) which describes the accessibility of the
%  file.
%  
% INPUT:
%  cdf is the name of a netCDF file, including the .cdf or .nc extent,
%  but without any .Z extent or it may be a url.
%
% OUTPUT:
%  status; a status flag; 
%  = 0 if the netCDF file is in the current directory.
%  = 1 if the file cannot be found anywhere.
%  = 2 if the file is in the directory specified by a call to the
%      m-function pos_cds.
%  = 3 if a compressed version of the file is in the current directory.
%  = 4 if this is a url.
%
% EXAMPLE:
%  status = check_nc('fred.nc')
%
% Note: check_nc only checks whether a file of the given name exists or
%       whether cdl looks like a url. It does not actually read the file or
%       try to access the url.
%
% This function calls: pos_cds.m
% This function is called by: attnc.m, getnc.m, getnc_s.m, inqnc.m, timenc.m,
%
% AUTHOR:   J. V. Mansbridge, CSIRO

%     Copyright (C), 1992, J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: check_nc.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

% Check for a url

if strncmp(cdf, 'http://', 7)
  status = 4;
  return
end

% Try to find the netcdf file.

if exist(cdf) == 2
  status = 0;
else
  % temp = getenv('NETCDFDATAPATH');
  temp = pos_cds;
  data_set = [ temp cdf ];
  if exist(data_set) == 2
    status = 2;
  else
    comp = [cdf '.Z'];
    if exist(comp) == 2
      status = 3;
    else
      status = 1;
    end
  end
end

