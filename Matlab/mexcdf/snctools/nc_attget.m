function values = nc_attget(ncfile, varname, attribute_name )
% NC_ATTGET: Get the values of a NetCDF attribute.
%
% USAGE:  att_value = nc_attget(ncfile, varname, attribute_name);
%
% PARAMETERS:
% Input:
%   ncfile:  
%       name of netcdf file in question
%   varname:  
%       name of variable in question.  In order to retrieve a global
%       attribute, use NC_GLOBAL for the variable name argument.  Do
%       Do NOT use 'global'!
%   attribute_name:  
%       name of attribute in question
% Output:    
%   values:  
%       value of attribute asked for.  Returns the empty matrix 
%       in case of an error.  There is an ambiguity in the case of 
%       NC_BYTE data, so it is always retrieved as an int8 datatype.
%       If you wanted uint8, then you must cast it yourself.
%
% You can specify that java be used instead of the mex-file by setting
% the appropriate preference, i.e.
%     >> setpref('SNCTOOLS','USE_JAVA',true);
%
% Example:
%    >> values = nc_attget('foo.nc', 'x', 'scale_factor')
%
% Example:  retrieving a global attribute.  Note we don't use 
%    'nc_global' or 'global'. 
% 
%    >> history = nc_attget('foo.nc', nc_global, 'history')
%
% SEE ALSO:  NC_GLOBAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_attget.m 2688 2009-04-30 18:58:17Z johnevans007 $
% $LastChangedDate: 2009-04-30 14:58:17 -0400 (Thu, 30 Apr 2009) $
% $LastChangedRevision: 2688 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(nargchk(3,3,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

backend = snc_read_backend(ncfile);
switch(backend)
	case 'tmw'
		values = nc_attget_tmw(ncfile,varname,attribute_name);
	case 'java'
		values = nc_attget_java(ncfile,varname,attribute_name);
	case 'mexnc'
		values = nc_attget_mexnc(ncfile,varname,attribute_name);
	otherwise
		error('SNCTOOLS:NC_ATTGET:unhandledBackend', ...
		      '%s is not a recognized backend for SNCTOOLS.', ...
			  backend);
end


return



%--------------------------------------------------------------------------
function values = nc_attget_mexnc(ncfile, varname, attribute_name )

[ncid, status] =mexnc('open', ncfile, nc_nowrite_mode );
if ( status ~= 0 )
    ncerror = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_ATTGET:MEXNC:OPEN', ncerror );
end

switch class(varname)
case { 'double' }
    varid = varname;

case 'char'
    varid = figure_out_varid ( ncid, varname );

otherwise
    error ( 'SNCTOOLS:NC_ATTGET:badType', 'Must specify either a variable name or NC_GLOBAL' );

end


funcstr = determine_funcstr(ncid,varid,attribute_name);

%
% And finally, retrieve the attribute.
[values, status]=mexnc(funcstr,ncid,varid,attribute_name);
if ( status ~= 0 )
    ncerror = mexnc ( 'strerror', status );
    err_id = ['SNCTOOLS:NC_ATTGET:MEXNC:' funcstr ];
    error ( err_id, ncerror );
end

status = mexnc('close',ncid);
if ( status ~= 0 )
    ncerror = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_ATTGET:MEXNC:CLOSE', ncerror );
end


return;











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function funcstr = determine_funcstr(ncid,varid,attribute_name)
% This function is for the mex-file backend.  Determine which netCDF function
% string we invoke to retrieve the attribute value.

[dt, status]=mexnc('inq_atttype',ncid,varid,attribute_name);
if ( status ~= 0 )
    mexnc('close',ncid);
    ncerror = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_ATTGET:MEXNC:INQ_ATTTYPE', ncerror );
end

switch ( dt )
case nc_double
    funcstr = 'GET_ATT_DOUBLE';
case nc_float
    funcstr = 'GET_ATT_FLOAT';
case nc_int
    funcstr = 'GET_ATT_INT';
case nc_short
    funcstr = 'GET_ATT_SHORT';
case nc_byte
    funcstr = 'GET_ATT_SCHAR';
case nc_char
    funcstr = 'GET_ATT_TEXT';
otherwise
    mexnc('close',ncid);
    error ( 'SNCTOOLS:NC_ATTGET:badDatatype', 'Unhandled datatype ID %d', dt );
end

return





%===============================================================================
%
% Did the user do something really stupid like say 'global' when they meant
% NC_GLOBAL?
function varid = figure_out_varid ( ncid, varname )

if isempty(varname)
    varid = nc_global;
    return;
end

if ( strcmpi(varname,'global') )
    [varid, status] = mexnc ( 'inq_varid', ncid, varname ); %#ok<ASGLU>
    if status 
        %
        % Ok, the user meant NC_GLOBAL
        warning ( 'SNCTOOLS:nc_attget:doNotUseGlobalString', ...
                  'Please consider using the m-file NC_GLOBAL.M instead of the string ''%s''.', varname );
        varid = nc_global;
        return;
    end
end

[varid, status] = mexnc ( 'inq_varid', ncid, varname );
if ( status ~= 0 )
    mexnc('close',ncid);
    ncerror = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_ATTGET:MEXNC:INQ_VARID', ncerror );
end

