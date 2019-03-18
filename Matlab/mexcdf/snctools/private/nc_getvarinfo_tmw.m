function Dataset = nc_getvarinfo_tmw ( arg1, arg2 )

%
% If we are here, then we must have been given something local.
if ischar(arg1) && ischar(arg2)

	ncfile = arg1;
	varname = arg2;


	try
		ncid=netcdf.open(ncfile,nc_nowrite_mode);
		varid = netcdf.inqVarid(ncid, varname);
		Dataset = get_varinfo_tmw ( ncid,  varid );
		netcdf.close(ncid);
	catch myException
		if exist('ncid','var')
			netcdf.close(ncid);
			rethrow(myException);
		end
	end

elseif isnumeric ( arg1 ) && isnumeric ( arg2 )

	ncid = arg1;
	varid = arg2;

	Dataset = get_varinfo_tmw ( ncid,  varid );

else
	error ( 'SNCTOOLS:NC_GETVARINFO:tmw:badTypes', ...
	        'Must have either both character inputs, or both numeric.' );
end


return




