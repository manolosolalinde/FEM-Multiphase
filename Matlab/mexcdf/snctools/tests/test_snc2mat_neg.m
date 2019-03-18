function test_snc2mat_neg ( ncfile )

if nargin == 0
	ncfile = 'foo.nc';
end

test_file_does_not_exist ( ncfile );

return










%--------------------------------------------------------------------------
function test_file_does_not_exist ( ncfile )

% netcdf file does not exist.
matfile_name = [ ncfile '.mat' ];
try
	snc2mat ( 'bad.nc', matfile_name );
catch me %#ok<NASGU>
    %  'MATLAB:netcdf:open:noSuchFile'
    return
end








