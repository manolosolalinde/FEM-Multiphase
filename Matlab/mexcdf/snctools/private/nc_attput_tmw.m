%-----------------------------------------------------------------------
function nc_attput_tmw ( ncfile, varname, attribute_name, attval )

ncid  =netcdf.open(ncfile, nc_write_mode );

try
    netcdf.redef(ncid);

    if isnumeric(varname)
        varid = varname;
    else
        varid = netcdf.inqVarID(ncid, varname );
    end
    
    netcdf.putAtt(ncid,varid,attribute_name,attval);
    netcdf.endDef(ncid);
    netcdf.close(ncid);

catch myException
    netcdf.close(ncid);
    rethrow(myException);
end


return;
