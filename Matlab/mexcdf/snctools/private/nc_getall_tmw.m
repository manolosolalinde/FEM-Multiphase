function data = nc_getall_tmw ( ncfile )
data = [];

cdfid=netcdf.open(ncfile,'NOWRITE');


[dud,nvars,ngatts] = netcdf.inq(cdfid); %#ok<ASGLU>

for varid=0:nvars-1

    varstruct = [];

    [varname, datatype, dims, natts] = netcdf.inqVar(cdfid, varid); %#ok<ASGLU>
    ndims = numel(dims);


    %
    % If ndims is zero, then it must be a singleton variable.  Don't bother trying
    % to retrieve the data, there is none.
    if ( ndims == 0 )
        varstruct.data = [];
    else
        values = nc_varget(ncfile, varname);
        varstruct.data = values;
    end



    %
    % get all the attributes
    for attnum = 0:natts-1

        attname = netcdf.inqAttName(cdfid, varid, attnum);
        try
            attval = nc_attget(ncfile, varname, attname);
        catch ME
            netcdf.close(cdfid);
            error ( 'SNCTOOLS:nc_getall:attributeRetrievalFailed', ...
                'nc_attget failed, ''%s''.\n', ME.message );
        end
        
        %
        % Matlab structures don't like the leading '_'
        if strcmp(attname,'_FillValue' )
            attname = 'FillValue';
        end


        sanitized_attname = matlab_sanitize_attname ( attname );


        %
        % this puts the attribute into the variable structure
        varstruct.(sanitized_attname) = attval;


    end


    %
    % Add this variable to the entire file structure
    data.(varname) = varstruct;

end


%
% Do the same for the global attributes
%
% get all the attributes
global_atts = [];
for attnum = 0:ngatts-1

    attname = netcdf.inqAttName(cdfid, nc_global, attnum);
    try
        attval = nc_attget(ncfile, nc_global, attname);
    catch ME
        netcdf.close(cdfid);
        error ( 'MATLAB:nc_getall:globalAttributeRetrievalFailed', ...
            'Failed to retrieve attribute #%d, ''%s''', ...
            attnum, ME.message);
    end
    
    sanitized_attname = matlab_sanitize_attname ( attname );


    %
    % this puts the attribute into the variable structure
    global_atts.(sanitized_attname) = attval;



end

if ~isempty ( global_atts )
    data.global_atts = global_atts;
end

netcdf.close(cdfid);

if isempty(data)
    data = struct([]);
end

return



