function bool = nc_isunlimitedvar ( ncfile, varname )
% NC_ISUNLIMITEDVAR:  determines if a variable has an unlimited dimension
%
% BOOL = NC_ISUNLIMITEDVAR ( NCFILE, VARNAME ) returns TRUE if the netCDF
% variable VARNAME in the netCDF file NCFILE has an unlimited dimension, 
% and FALSE otherwise.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_isunlimitedvar.m 2659 2009-04-01 17:38:36Z johnevans007 $
% $LastChangedDate: 2009-04-01 13:38:36 -0400 (Wed, 01 Apr 2009) $
% $LastChangedRevision: 2659 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


error(nargchk(2,2,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

try
    DataSet = nc_getvarinfo ( ncfile, varname );
catch
    e = lasterror;
    switch ( e.identifier )
        case { 'SNCTOOLS:NC_GETVARINFO:badVariableName', ...
               'SNCTOOLS:NC_VARGET:MEXNC:INQ_VARID', ...
	       'MATLAB:netcdf:inqVarID:variableNotFound' }
            bool = false;
            return
        otherwise
            error('SNCTOOLS:NC_ISUNLIMITEDVAR:unhandledCondition', e.message );
    end
end

bool = DataSet.Unlimited;

return;
