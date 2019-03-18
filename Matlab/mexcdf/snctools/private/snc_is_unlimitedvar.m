function tf = snc_is_unlimitedvar(ncfile,varname)     

% Make use of the error structure when using 2007b and higher.
try
    DataSet = nc_getvarinfo ( ncfile, varname );
catch me
    
    switch ( me.identifier )
        case { 'SNCTOOLS:NC_GETVARINFO:badVariableName', ...
                'SNCTOOLS:NC_VARGET:MEXNC:INQ_VARID', ...
                'MATLAB:netcdf:inqVarID:variableNotFound', ...
                'MATLAB:netcdf:inqVarID:enotvar:variableNotFound' }
            tf = false;
            return
        otherwise
            error('SNCTOOLS:NC_ISUNLIMITEDVAR:unhandledCondition', me.message );
    end
end

tf = DataSet.Unlimited;