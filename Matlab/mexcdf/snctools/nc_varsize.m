function varsize = nc_varsize(ncfile, varname)
% NC_VARSIZE:  return the size of the requested netncfile variable
%
% VARSIZE = NC_VARSIZE(NCFILE,NCVAR) returns the size of the netCDF variable 
% NCVAR in the netCDF file NCFILE.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_varsize.m 2659 2009-04-01 17:38:36Z johnevans007 $
% $LastChangedDate: 2009-04-01 13:38:36 -0400 (Wed, 01 Apr 2009) $
% $LastChangedRevision: 2659 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(2,2,nargin,'struct'));
error(nargoutchk(1,1,nargout,'struct'));

if ~ischar(ncfile)
	error ( 'SNCTOOLS:NC_VARSIZE:badInputType', 'The input filename must be a string.' );
end
if ~ischar(varname)
	error ( 'SNCTOOLS:NC_VARSIZE:badInputType', 'The input variable name must be a string.' );
end


v = nc_getvarinfo ( ncfile, varname );

varsize = v.Size;

return

