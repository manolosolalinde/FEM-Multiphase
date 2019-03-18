% TMEXNC Test of mexnc Mex-file interface to NetCDF.

% Copyright (C) 1993 Dr. Charles R. Denham, ZYDECO.
% All Rights Reserved.

status = mexnc('setopts', 0);

for i = 0:10
   status = mexnc('close', i);
end

% Phase 1: Simple creation and fill.

cdf = mexnc('create', 'mycdf.cdf', 'clobber');

xdim = mexnc('dimdef', cdf, 'x', 10);
ydim = mexnc('dimdef', cdf, 'y', 10);

avar = mexnc('vardef', cdf, 'a', 'double', 1, [xdim]);
bvar = mexnc('vardef', cdf, 'b', 'double', 1, [ydim]);

status = mexnc('close', cdf);

cdf = mexnc('open', 'mycdf.cdf', 'write');

status = mexnc('redef', cdf);

cvar = mexnc('vardef', cdf, 'c', 'double', 2, [xdim ydim]);

status = mexnc('endef', cdf);

a = rand(10, 1);
b = ones(10, 1) * pi;

status = mexnc('varput', cdf, 'a', [0], [10], a);
status = mexnc('varput', cdf, 'b', [0], [10], b);

status = mexnc('close', cdf);

disp('TMEXNC: Phase 1 done.')

% Phase 2: foo.cdf example, NetCDF User's Guide, Chapter 10.

cdf = mexnc('create', 'foo.cdf', 'clobber');

dlat = mexnc('dimdef', cdf, 'lat', 10);
dlon = mexnc('dimdef', cdf, 'lon', 5);
dtime = mexnc('dimdef', cdf, 'time', 'unlimited');

vlat = mexnc('vardef', cdf, 'lat', 'long', 1, dlat);
vlon = mexnc('vardef', cdf, 'lon', 'long', 1, dlon);
vtime = mexnc('vardef', cdf, 'time', 'long', 3, [dtime dlat dlon]);
vz = mexnc('vardef', cdf, 'z', 'float', 3, [dtime dlat dlon]);
vt = mexnc('vardef', cdf, 't', 'float', 3, [dtime dlat dlon]);
vp = mexnc('vardef', cdf, 'p', 'double', 3, [dtime dlat dlon]);
vrh = mexnc('vardef', cdf, 'rh', 'long', 3, [dtime dlat dlon]);

status = mexnc('endef', cdf);

status = mexnc('attput', cdf, vlat, 'units', 'char', -1, 'degrees_north');
status = mexnc('attput', cdf, vlon, 'units', 'char', -1, 'degrees_east');
status = mexnc('attput', cdf, vtime, 'units', 'char', -1, 'seconds');
status = mexnc('attput', cdf, vz, 'units', 'char', -1, 'meters');
status = mexnc('attput', cdf, vz, 'valid_range', 'float', -1, [0 5000]);
status = mexnc('attput', cdf, vp, '_FillValue', 'double', -1, -9999);
status = mexnc('attput', cdf, vrh, '_FillValue', 'long', -1, -1);

lat = [0 10 20 30 40 50 60 70 80 90];
lon = [-140 -118 -96 -84 -52];

status = mexnc('varput', cdf, vlat, 0, 10, lat);
status = mexnc('varput', cdf, vlon, 0, 5, lon);

status = mexnc('close', cdf);

disp('TMEXNC: Phase 2 done.')

% Phase 3: bar.cdf for NC_UNLIMITED data.

cdf = mexnc('create', 'bar.cdf', 'clobber');

dimid = mexnc('dimdef', cdf, 'i', 'unlimited');

varid = mexnc('vardef', cdf, 'x', 'double', 1, dimid);

status = mexnc('endef', cdf);

status = mexnc('close', cdf);

[x, j] = sort(rand(11, 1));
j = j(:).' - 1;
disp(j)

for i = j
   cdf = mexnc('open', 'bar.cdf', 'write');
   varid = mexnc('varid', cdf, 'x');
   status = mexnc('varput', cdf, varid, i, 1, i);
   x = mexnc('varget', cdf, 'x', i, 1);
   if x ~= i
      disp(['Bad put/get: ', int2str(i), num2str(x)]);
   end
   status = mexnc('close', cdf);
end

disp('TMEXNC: Phase 3 done.')

