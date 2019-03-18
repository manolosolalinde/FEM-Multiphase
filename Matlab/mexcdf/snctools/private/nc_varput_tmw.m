function nc_varput_tmw( ncfile, varname, data, start,count,stride )

try
    ncid = netcdf.open(ncfile, nc_write_mode);
    varid = netcdf.inqVarID(ncid, varname );
    [dud,var_type,var_dim]=netcdf.inqVar(ncid,varid); %#ok<ASGLU>
    nvdims = numel(var_dim);
    
    v = nc_getvarinfo(ncfile,varname);
    nc_count = v.Size;
    
    [start, count] = nc_varput_validate_indexing(ncid,nvdims,data,start,count,stride);
    
    %
    % check that the length of the start argument matches the rank of the variable.
    if numel(start) ~= numel(nc_count)
        netcdf.close( ncid );
        error ( 'SNCTOOLS:NC_VARPUT:badIndexing', ...
                'Length of START index (%d) does not make sense with a variable rank of %d.\n', ...
                 numel(start), numel(count) );
    end
    
    data = handle_fill_value_tmw ( ncid, varid, data );
    data = handle_scaling_tmw(ncid,varid,data);

	preserve_fvd = getpref('SNCTOOLS','PRESERVE_FVD',false);
    if ~preserve_fvd
        data = permute(data,fliplr(1:ndims(data)));
        start = fliplr(start);
        count = fliplr(count);
        stride = fliplr(stride);
    end
    
    if isempty(start) || (nvdims == 0)
        netcdf.putVar(ncid,varid,data);
    elseif isempty(count)
        netcdf.putVar(ncid,varid,start,data);
    elseif isempty(stride)
        netcdf.putVar(ncid,varid,start,count,data);
    else
        netcdf.putVar(ncid,varid,start,count,stride,data);
    end
    
    
    netcdf.close(ncid);

catch myException
    if exist('ncid','var')
        netcdf.close(ncid);
        rethrow(myException);
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_scaling_tmw(ncid,varid,data)
% HANDLE_SCALING_TMW
%     If there is a scale factor and/or  add_offset attribute, convert the data
%     to double precision and apply the scaling.
%

have_scale_factor = 0;
have_add_offset = 0;

try
    netcdf.inqAtt(ncid, varid, 'scale_factor' );
    have_scale_factor = 1;
catch %#ok<CTCH>
    
end

try
    netcdf.inqAtt(ncid, varid, 'add_offset' );
    have_add_offset = 1;
catch %#ok<CTCH>
    
end

%
% Return early if we don't have either one.
if ~(have_scale_factor || have_add_offset)
    return;
end

scale_factor = 1.0;
add_offset = 0.0;

try

    if have_scale_factor
        scale_factor = netcdf.getAtt(ncid, varid, 'scale_factor','double' );
    end
    
    if have_add_offset
        add_offset = netcdf.getAtt(ncid, varid, 'add_offset','double' );
    end
    
    data = (double(data) - add_offset) / scale_factor;
    
    %
    % When scaling to an integer, we should add 0.5 to the data.  Otherwise
    % there is a tiny loss in precision, e.g. 82.7 should round to 83, not 
    % 82.
    [varname,xtype] = netcdf.inqVar(ncid,varid);  %#ok<ASGLU>
    switch xtype
        case { nc_int, nc_short, nc_byte, nc_char }
            data = round(data);
    end

catch myException
    netcdf.close(ncid);
    rethrow(myException);
end

return






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = handle_fill_value_tmw(ncid,varid,data)

%
% Handle the fill value.  We do this by changing any NaNs into
% the _FillValue.  That way the netcdf library will recognize it.
try
    switch ( class(data) )
    case 'double'
        myClass = 'double';
    case 'single'
        myClass = 'float';
    case 'int32'
        myClass = 'int';
    case 'int16'
        myClass = 'short';
    case 'int8'
        myClass = 'schar';
    case 'uint8'
        myClass = 'uchar';
    case 'char'
        myClass = 'text';
    otherwise
        netcdf.close(ncid);
        error ( 'SNCTOOLS:NC_VARPUT:unhandledDatatype', ...
                'Unhandled datatype for fill value, ''%s''.', ...
                class(data) );
    end

    fill_value  = netcdf.getAtt(ncid,varid,'_FillValue',myClass);

    data(isnan(data)) = fill_value;

catch myException %#ok<NASGU>
    return
end


    











