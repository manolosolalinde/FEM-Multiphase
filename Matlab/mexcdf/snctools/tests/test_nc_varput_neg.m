function test_nc_varput_neg (  )

test_netcdf3;
test_hdf4;
test_netcdf4;

return



%--------------------------------------------------------------------------
function test_netcdf3()
fprintf('\tRunning netcdf-3 tests...  ' );
testroot = fileparts(mfilename('fullpath'));

ncfile = fullfile(testroot,'testdata/varput.nc');
run_generic_tests(ncfile);
run_singleton_tests(ncfile);

% This doesn't work for nc4 or hdf4
test_neg_2d_to_singleton ( ncfile );

fprintf('OK\n');
return

%--------------------------------------------------------------------------
function test_hdf4()
fprintf('\tRunning hdf4 tests...  ' );
testroot = fileparts(mfilename('fullpath'));

ncfile = fullfile(testroot,'testdata/varput.hdf');
run_generic_tests(ncfile);

fprintf('OK\n');
return


%--------------------------------------------------------------------------
function test_netcdf4()

if ~netcdf4_capable
    fprintf('\tmexnc (netcdf-4) backend testing filtered out on ');
    fprintf('configurations where the library version < 4.\n');
    return
end

fprintf('\tRunning netcdf-4 tests...' );
testroot = fileparts(mfilename('fullpath'));
test_no_input_arguments;

ncfile = fullfile(testroot,'testdata/empty-4.nc');
test_only_one_argument ( ncfile );
test_only_two_arguments ( ncfile );
test_bad_filename ('i_do_not_exist.nc');

ncfile = fullfile(testroot,'testdata/varput4.nc');
run_generic_tests(ncfile);
run_singleton_tests(ncfile);

fprintf('OK\n');
return


%--------------------------------------------------------------------------
function run_singleton_tests(input_ncfile)

ncfile = 'foo.nc';
copyfile(input_ncfile,ncfile);

test_singleton_bad_start ( ncfile );
test_singleton_bad_count ( ncfile );
test_singleton_with_stride_which_is_bad ( ncfile );

return

%--------------------------------------------------------------------------
function run_generic_tests(input_ncfile)

ncfile = 'foo.nc';
copyfile(input_ncfile,ncfile);



test_write_1D_size_mismatch ( ncfile );
test_write_1D_bad_stride ( ncfile );


test_write_2D_strided_bad_start ( ncfile );
test_write_2D_chunk_bad_count ( ncfile );
test_write_2D_bad_stride ( ncfile );


return






















%--------------------------------------------------------------------------
function test_neg_2d_to_singleton ( ncfile )
% Don't test on 2007a or earlier

global ignore_eids
try
    nc_varput ( ncfile, 'test_singleton', [2 1] );
catch me 
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'MATLAB:netcdf:putVar:dataSizeMismatch', ...
                'SNCTOOLS:NC_VARPUT:MEXNC:varput:dataSizeMismatch', ...
                'MATLAB:netcdf:open:notANetcdfFile'}
            return
        otherwise
            rethrow(me);
    end
end
error('nc_varput succeeded when it should not have.');








%--------------------------------------------------------------------------
function test_write_1D_size_mismatch ( ncfile )


global ignore_eids

input_data = 3.14159;
try
    nc_varput ( ncfile, 'test_1D', input_data, 4, 2 );
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'MATLAB:netcdf:putVara:dataSizeMismatch', ...
                'SNCTOOLS:NC_VARPUT:MEXNC:putVara:dataSizeMismatch', ...
                'SNCTOOLS:varput:hdf4:writedataFailed'}
            return
        otherwise
            rethrow(me);
    end
end        


return







%--------------------------------------------------------------------------
function test_write_1D_bad_stride ( ncfile )


global ignore_eids;

input_data = [3.14159; 2];
try
    nc_varput ( ncfile, 'test_1D', input_data, 0, 2, 8 );
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'MATLAB:netcdf:putVars:indexExceedsDimensionBound', ...
                'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
                'SNCTOOLS:varput:hdf4:writedataFailed'}
            return
        otherwise
            rethrow(me);
    end
end
    
error('nc_varput succeeded when it should have failed.');






%--------------------------------------------------------------------------
function test_singleton_bad_start ( ncfile )

global ignore_eids;
input_data = 3.14159;
try
    nc_varput ( ncfile, 'test_singleton', input_data, 4, 1 );
catch me

    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:badIndexing' }
            return
        otherwise
            rethrow(me);
    end

end






%--------------------------------------------------------------------------
function test_singleton_bad_count ( ncfile )

global ignore_eids;
input_data = 3.14159;
try
    nc_varput ( ncfile, 'test_singleton', input_data, 0, 2 );
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:badIndexing' }
            return
        otherwise
            rethrow(me);
    end    

end










%--------------------------------------------------------------------------
function test_singleton_with_stride_which_is_bad ( ncfile )

global ignore_eids;

input_data = 3.14159;
try
    nc_varput ( ncfile, 'test_singleton', input_data, 0, 1, 1 );
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:badIndexing' }
            return
        otherwise
            rethrow(me);
    end   

end

return







%--------------------------------------------------------------------------
function test_write_2D_too_much_with_putvar ( ncfile )

global ignore_eids;

input_data = 1:49;
input_data = reshape(input_data,7,7);
try
    nc_varput ( ncfile, 'test_2D', input_data );
catch me
  
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
                'MATLAB:netcdf:putVara:startPlusCountExceedsDimensionBound', ...
                'SNCTOOLS:varput:hdf4:writedataFailed' }
            return
        otherwise
            rethrow(me);
    end   

end
error('failed');







%--------------------------------------------------------------------------
function test_write_2D_chunk_bad_offset ( ncfile )
% write with a bad offset

global ignore_eids;

sz = nc_varsize(ncfile,'test_2D');
start = [1 1];
count = sz;

input_data = 1:prod(count);
input_data = reshape(input_data,count);
try
    nc_varput ( ncfile, 'test_2D', input_data, start, count );
catch me
    
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
                'MATLAB:netcdf:putVara:startPlusCountExceedsDimensionBound', ...
                'SNCTOOLS:varput:hdf4:writedataFailed'}
            return
        otherwise
            rethrow(me);
    end  
end
error('failed');








%--------------------------------------------------------------------------
function test_write_2D_strided_bad_start ( ncfile )
% write using put_vars with a bad offset

global ignore_eids;

sz = nc_varsize(ncfile,'test_2D');
start = [2 1];
count = sz/2;
stride = [2 2];

input_data = (1:prod(count)) + 3.14159;
input_data = reshape(input_data,count);

try
    nc_varput ( ncfile, 'test_2D', input_data, start, count, stride);
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
                'MATLAB:netcdf:putVars:indexExceedsDimensionBound', ...
                'SNCTOOLS:varput:hdf4:writedataFailed'}
            return
        otherwise
            rethrow(me);
    end  
end
error('failed');







%--------------------------------------------------------------------------
function test_write_2D_chunk_bad_count ( ncfile )
% vara with bad count

global ignore_eids;

sz = nc_varsize(ncfile,'test_2D');
start = [0 0];
count = sz+1;

input_data = (1:prod(count)) + 3.14159;
input_data = reshape(input_data,count);
try
    nc_varput ( ncfile, 'test_2D', input_data, start, count );
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
                'MATLAB:netcdf:putVara:startPlusCountExceedsDimensionBound', ...
                'SNCTOOLS:varput:hdf4:writedataFailed'}
            return
        otherwise
            rethrow(me);
    end 
end
error('failed');







%--------------------------------------------------------------------------
function test_write_2D_bad_stride ( ncfile )

global ignore_eids;

sz = nc_varsize(ncfile,'test_2D');
start = [0 0];
count = sz/2;
stride = [3 3];

input_data = (1:prod(count)) + 3.14159;
input_data = reshape(input_data,count);
try
    nc_varput ( ncfile, 'test_2D', input_data, start, count, stride);
catch me
    if ignore_eids
        return
    end
    switch(me.identifier)
        case { 'SNCTOOLS:NC_VARPUT:writeOperationFailed', ...
                'MATLAB:netcdf:putVars:indexExceedsDimensionBound', ...
                'SNCTOOLS:varput:hdf4:writedataFailed'}
            return
        otherwise
            rethrow(me);
    end  
end
error('failed');




%--------------------------------------------------------------------------
function test_no_input_arguments()

try
    nc_varput;
catch me %#ok<NASGU>
    %  'MATLAB:nargchk:notEnoughInputs'
	return
end
error('nc_varput succeeded when it should not have.');


%--------------------------------------------------------------------------
function test_only_one_argument ( ncfile )
try
    nc_varput ( ncfile );
catch %#ok<CTCH>
	return
end
error('nc_varput succeeded when it should not have.');

%--------------------------------------------------------------------------
function test_only_two_arguments ( ncfile )

try
    nc_varput ( ncfile, 'test_2d' );
catch %#ok<CTCH>
	return
end
error('nc_varput succeeded when it should not have.');

%--------------------------------------------------------------------------
function test_bad_filename ( ncfile )

try
    nc_varput ( ncfile, 'test_2d', rand(5,5) );
catch me  %#ok<NASGU>
    % 'MATLAB:netcdf:open:noSuchFile'
	return
end
error('nc_varput succeeded when it should not have.');





