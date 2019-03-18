function bool = nc_isvar_tmw ( ncfile, varname )

ncid = netcdf.open(ncfile, nc_nowrite_mode );
try
	varid = netcdf.inqVarID(ncid,varname);
	bool = true;
catch myException
	bool = false;
end

netcdf.close(ncid);
return









