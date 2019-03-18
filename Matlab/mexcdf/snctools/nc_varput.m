function nc_varput( ncfile, varname, data, varargin )
% NC_VARPUT:  Writes data into a netCDF file.
%
% NC_VARPUT(NCFILE,VARNAME,DATA) writes the matlab variable DATA to
% the variable VARNAME in the netCDF file NCFILE.  The main requirement
% here is that DATA have the same dimensions as the netCDF variable.
%
% NC_VARPUT(NCFILE,VARNAME,DATA,START,COUNT) writes DATA contiguously, 
% starting at the zero-based index START and with extents given by
% COUNT.
%
% NC_VARPUT(NCFILE,VARNAME,DATA,START,COUNT,STRIDE) writes DATA  
% starting at the zero-based index START with extents given by
% COUNT, but this time with strides given by STRIDE.  If STRIDE is not
% given, then it is assumes that all data is contiguous.
%
% Setting the preference 'PRESERVE_FVD' to true will compel MATLAB to 
% display the dimensions in the opposite order from what the C utility 
% ncdump displays.  
% 
% EXAMPLES:
%    Suppose you have a netcdf variable called 'x' of size 6x4.  If you 
%    have an array of data called 'mydata' that is 6x4, then you can 
%    write to the entire variable with 
% 
%        >> nc_varput ( 'foo.nc', 'x', mydata );
%
%    If you wish to only write to the first 2 rows and three columns,
%    you could do the following
%
%        >> subdata = mydata(1:2,1:3);
%        >> nc_varput ( 'foo.nc', 'x', subdata, [0 0], [2 3] );
%
%
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: nc_varput.m 2691 2009-05-01 14:11:15Z johnevans007 $
% $LastChangedDate: 2009-05-01 10:11:15 -0400 (Fri, 01 May 2009) $
% $LastChangedRevision: 2691 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


error(nargchk(3,6,nargin));
error(nargoutchk(0,0,nargout));

[start, count, stride] = parse_and_validate_args(ncfile,varname,varargin{:});

backend = snc_write_backend(ncfile);
switch(backend)
	case 'tmw'
		nc_varput_tmw(ncfile,varname,data,start,count,stride);
	case 'java'
		nc_varput_java(ncfile,varname,data,start,count,stride);
	case 'mexnc'
		nc_varput_mexnc(ncfile,varname,data,start,count,stride);
end

return





%-----------------------------------------------------------------------
function nc_varput_mexnc( ncfile, varname, data, start,count,stride )


[ncid, status] = mexnc('open', ncfile, nc_write_mode);
if (status ~= 0)
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARPUT:MEXNC:OPEN', ncerr );
end




%
% check to see if the variable already exists.  
[varid, status] = mexnc('INQ_VARID', ncid, varname );
if ( status ~= 0 )
    mexnc ( 'close', ncid );
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARPUT:MEXNC:INQ_VARID', ncerr );
end


[dud,var_type,nvdims,var_dim,dud, status]=mexnc('INQ_VAR',ncid,varid); %#ok<ASGLU>
if status ~= 0 
    mexnc ( 'close', ncid );
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARPUT:MEXNC:INQ_VAR', ncerr );
end


v = nc_getvarinfo ( ncfile, varname );
nc_count = v.Size;


[start, count] = nc_varput_validate_indexing (ncid,nvdims,data,start,count,stride);

%
% check that the length of the start argument matches the rank of the variable.
if length(start) ~= length(nc_count)
    mexnc ( 'close', ncid );
    error ( 'SNCTOOLS:NC_VARPUT:badIndexing', ...
              'Length of START index (%d) does not make sense with a variable rank of %d.\n', ...
              length(start), length(nc_count) );
end



%
% Figure out which write routine we will use.  If the target variable is a singleton, then we must use
% VARPUT1.  If a stride was given, we must use VARPUTG.  Otherwise just use VARPUT.
if isempty(start) || (nvdims == 0)
    write_op = 'put_var';
elseif isempty(count)
    write_op = 'put_var1';
elseif isempty(stride)
    write_op = 'put_vara';
else
    write_op = 'put_vars';
end



data = handle_fill_value ( ncid, varid, data );
data = handle_scaling(ncid,varid,data);

preserve_fvd = getpref('SNCTOOLS','PRESERVE_FVD',false);
if preserve_fvd
    start = fliplr(start);
    count = fliplr(count);
    stride = fliplr(stride);
else
    data = permute(data,fliplr(1:ndims(data)));
end

write_the_data(ncid,varid,start,count,stride,write_op,data);


status = mexnc ( 'close', ncid );
if ( status ~= 0 )
    error ( 'SNCTOOLS:nc_varput:close', mexnc('STRERROR',status));
end


return




%-------------------------------------------------------------------------
function [start, count, stride] = parse_and_validate_args(ncfile,varname,varargin)

%
% Set up default outputs.
start = [];
count = [];
stride = [];


switch length(varargin)
case 2
    start = varargin{1};
    count = varargin{2};

case 3
    start = varargin{1};
    count = varargin{2};
    stride = varargin{3};

end



%
% Error checking on the inputs.
if ~ischar(ncfile)
    error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the filename must be character.' );
end
if ~ischar(varname)
    error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the variable name must be character.' );
end

if ~isnumeric ( start )
    error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the ''start'' argument must be numeric.' );
end
if ~isnumeric ( count )
    error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the ''count'' argument must be numeric.' );
end
if ~isnumeric ( stride )
    error ( 'SNCTOOLS:NC_VARPUT:badInput', 'the ''stride'' argument must be numeric.' );
end


return








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_scaling(ncid,varid,data)
% HANDLE_MEX_SCALING
%     If there is a scale factor and/or  add_offset attribute, convert the data
%     to double precision and apply the scaling.
%

[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'scale_factor' ); %#ok<ASGLU>
if ( status == 0 )
    have_scale_factor = 1;
else
    have_scale_factor = 0;
end
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, 'add_offset' ); %#ok<ASGLU>
if ( status == 0 )
    have_add_offset = 1;
else
    have_add_offset = 0;
end

%
% Return early if we don't have either one.
if ~(have_scale_factor || have_add_offset)
    return;
end

scale_factor = 1.0;
add_offset = 0.0;


if have_scale_factor
    [scale_factor, status] = mexnc ( 'get_att_double', ncid, varid, 'scale_factor' );
    if ( status ~= 0 )
        mexnc ( 'close', ncid );
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARPUT:MEXNC:GET_ATT_DOUBLE', ncerr );
    end
end

if have_add_offset
    [add_offset, status] = mexnc ( 'get_att_double', ncid, varid, 'add_offset' );
    if ( status ~= 0 )
        mexnc ( 'close', ncid );
        ncerr = mexnc('strerror', status);
        error ( 'SNCTOOLS:NC_VARPUT:MEXNC:GET_ATT_DOUBLE', ncerr );
    end
end

[var_type,status]=mexnc('INQ_VARTYPE',ncid,varid);
if status ~= 0 
    mexnc ( 'close', ncid );
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARPUT:MEXNC:INQ_VARTYPE', ncerr );
end

data = (double(data) - add_offset) / scale_factor;

%
% When scaling to an integer, we should add 0.5 to the data.  Otherwise
% there is a tiny loss in precision, e.g. 82.7 should round to 83, not 
% .
switch var_type
    case { nc_int, nc_short, nc_byte, nc_char }
        data = round(data);
end


return






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_fill_value(ncid,varid,data)

%
% Handle the fill value.  We do this by changing any NaNs into
% the _FillValue.  That way the netcdf library will recognize it.
[dud, dud, status] = mexnc('INQ_ATT', ncid, varid, '_FillValue' ); %#ok<ASGLU>
if ( status == 0 )

    switch ( class(data) )
    case 'double'
        funcstr = 'get_att_double';
    case 'single'
        funcstr = 'get_att_float';
    case 'int32'
        funcstr = 'get_att_int';
    case 'int16'
        funcstr = 'get_att_short';
    case 'int8'
        funcstr = 'get_att_schar';
    case 'uint8'
        funcstr = 'get_att_uchar';
    case 'char'
        funcstr = 'get_att_text';
    otherwise
        mexnc ( 'close', ncid );
        error ( 'SNCTOOLS:NC_VARPUT:unhandledDatatype', ...
            'Unhandled datatype for fill value, ''%s''.', ...
            class(data) );
    end

    [fill_value, status] = mexnc(funcstr,ncid,varid,'_FillValue' );
    if ( status ~= 0 )
        mexnc ( 'close', ncid );
        ncerr = mexnc('strerror', status);
        err_id = [ 'SNCTOOLS:NC_VARPUT:MEXNC:' funcstr ];
        error ( err_id, ncerr );
    end


    data(isnan(data)) = fill_value;

end

    













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_the_data(ncid,varid,start,count,stride,write_op,pdata)

%
% write the data
switch ( write_op )

    case 'put_var1'
        switch ( class(pdata) )
        case 'double'
            funcstr = 'put_var1_double';
        case 'single'
            funcstr = 'put_var1_float';
        case 'int32'
            funcstr = 'put_var1_int';
        case 'int16'
            funcstr = 'put_var1_short';
        case 'int8'
            funcstr = 'put_var1_schar';
        case 'uint8'
            funcstr = 'put_var1_uchar';
        case 'char'
            funcstr = 'put_var1_text';
        otherwise
            mexnc('close',ncid);
            error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', ...
                'unhandled data class %s\n', ...
                class(pdata));
        end
        status = mexnc (funcstr, ncid, varid, start, pdata );

    case 'put_var'
        switch ( class(pdata) )
        case 'double'
            funcstr = 'put_var_double';
        case 'single'
            funcstr = 'put_var_float';
        case 'int32'
            funcstr = 'put_var_int';
        case 'int16'
            funcstr = 'put_var_short';
        case 'int8'
            funcstr = 'put_var_schar';
        case 'uint8'
            funcstr = 'put_var_uchar';
        case 'char'
            funcstr = 'put_var_text';
        otherwise
            mexnc('close',ncid);
            error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', ...
                'unhandled data class %s\n', class(pdata)  );
        end
        status = mexnc (funcstr, ncid, varid, pdata );
    
    case 'put_vara'
        switch ( class(pdata) )
        case 'double'
            funcstr = 'put_vara_double';
        case 'single'
            funcstr = 'put_vara_float';
        case 'int32'
            funcstr = 'put_vara_int';
        case 'int16'
            funcstr = 'put_vara_short';
        case 'int8'
            funcstr = 'put_vara_schar';
        case 'uint8'
            funcstr = 'put_vara_uchar';
        case 'char'
            funcstr = 'put_vara_text';
        otherwise
            mexnc('close',ncid);
            error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType',...
                'unhandled data class %s\n', class(pdata) );
        end
        status = mexnc (funcstr, ncid, varid, start, count, pdata );

    case 'put_vars'
        switch ( class(pdata) )
        case 'double'
            funcstr = 'put_vars_double';
        case 'single'
            funcstr = 'put_vars_float';
        case 'int32'
            funcstr = 'put_vars_int';
        case 'int16'
            funcstr = 'put_vars_short';
        case 'int8'
            funcstr = 'put_vars_schar';
        case 'uint8'
            funcstr = 'put_vars_uchar';
        case 'char'
            funcstr = 'put_vars_text';
        otherwise
            mexnc('close',ncid);
            error ( 'SNCTOOLS:NC_VARPUT:unhandledMatlabType', ...
                'unhandled data class %s\n', class(pdata) );
        end
        status = mexnc (funcstr, ncid, varid, start, count, stride, pdata );

    otherwise 
        mexnc ( 'close', ncid );
        error ( 'SNCTOOLS:NC_VARPUT:unhandledWriteOp', ...
            'unknown write operation''%s''.\n', write_op );


end

if ( status ~= 0 )
    mexnc ( 'close', ncid );
    ncerr = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
        'write operation ''%s'' failed with error ''%s''.', ...
        write_op, ncerr);
end

return
