function test_snctools()
% TEST_SNCTOOLS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Id: test_snctools.m 2661 2009-04-10 19:05:08Z johnevans007 $
% $LastChangedDate: 2009-04-10 15:05:08 -0400 (Fri, 10 Apr 2009) $
% $LastChangedRevision: 2661 $
% $LastChangedBy: johnevans007 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% switch off some warnings
mver = version('-release');
switch mver
    case {'11', '12'}
        error ( 'This version of MATLAB is too old, SNCTOOLS will not run.' );
    case {'13'}
        error ( 'R13 is not supported in this release of SNCTOOLS');
    otherwise
        warning('off', 'SNCTOOLS:nc_archive_buffer:deprecated' );
        warning('off', 'SNCTOOLS:nc_datatype_string:deprecated' );
        warning('off', 'SNCTOOLS:nc_diff:deprecated' );
        warning('off', 'SNCTOOLS:nc_getall:deprecated' );
        warning('off', 'SNCTOOLS:snc2mat:deprecated' );
end


run_backend_neutral_tests;
run_backend_mex_tests;

fprintf ( 1, '\nAll  possible tests for your configuration have been run.  Bye.\n\n' );

return




%----------------------------------------------------------------------
function run_mexnc_tests()

% Is mexnc ok?
mexnc_loc = which ( 'mexnc' );
mexnc_ok = ~isempty(which('mexnc'));

pause_duration = 3;
if ~mexnc_ok
    fprintf ( 1, 'MEXNC was not found, so the tests requiring mexnc\n' );
    fprintf ( 1, 'will not be run.\n\n' );
    return
end

fprintf ( 1, '\n' );
fprintf ( 1, 'Ok, we found mexnc.  ' );
fprintf ( 1, 'Remote OPeNDAP/mexnc tests ' );
if getpref('SNCTOOLS','TEST_REMOTE_MEXNC',false)
    fprintf ( 1, 'will ' );
    setpref('SNCTOOLS','TEST_REMOTE',true)
else
    fprintf ( 1, 'will NOT ' );
    setpref('SNCTOOLS','TEST_REMOTE',false)
end
fprintf ( 1, 'be run.\n  Starting tests in ' );
for j = 1:pause_duration
    fprintf ( 1, '%d... ', pause_duration - j + 1 );
    pause(1);
end
fprintf ( 1, '\n' );

run_backend_neutral_tests;
run_backend_mexnc_tests;


return


%----------------------------------------------------------------------
function run_all_tests()

fprintf ( 1, 'Ok, about to start testing in  ' );
pause_duration = 3;
for j = 1:pause_duration
    fprintf ( 1, '%d... ', pause_duration - j + 1 );
    pause(1);
end
fprintf ( 1, '\n' );

test_nc_attget;
test_nc_datatype_string;
test_nc_iscoordvar;
test_nc_isunlimitedvar;
test_nc_dump;
test_nc_getlast;
test_nc_isvar;
test_nc_varsize;
test_nc_getvarinfo;
test_nc_info;
test_nc_getbuffer;
test_nc_varget;
test_nc_getdiminfo;

test_nc_varput           ( 'test.nc' );
test_nc_add_dimension    ( 'test.nc' );
test_nc_addhist          ( 'test.nc' );
test_nc_addvar           ( 'test.nc' );
test_nc_attput           ( 'test.nc' );
test_nc_create_empty     ( 'test.nc' );
test_nc_varrename        ( 'test.nc' );
test_nc_addnewrecs       ( 'test.nc' );
test_nc_add_recs         ( 'test.nc' );
test_nc_archive_buffer   ( 'test.nc' );

test_snc2mat             ( 'test.nc' );
test_nc_getall           ( 'test.nc' );
test_nc_diff             ( 'test1.nc', 'test2.nc' );
test_nc_cat_a;


return




%----------------------------------------------------------------------
function run_tmw_tests()

fprintf ( 1, 'Ok, about to start TMW testing in  ' );
pause_duration = 3;
for j = 1:pause_duration
    fprintf ( 1, '%d... ', pause_duration - j + 1 );
    pause(1);
end
fprintf ( 1, '\n' );

run_backend_neutral_tests;
run_backend_mexnc_tests;

return





%----------------------------------------------------------------------
function run_backend_neutral_tests()

test_nc_attget;
test_nc_datatype_string;
test_nc_iscoordvar;
test_nc_isunlimitedvar;
test_nc_dump;
test_nc_getlast;
test_nc_isvar;
test_nc_varsize;
test_nc_getvarinfo;
test_nc_info;
test_nc_getbuffer;
test_nc_varget;
test_nc_getdiminfo;


return




%----------------------------------------------------------------------
function run_backend_mex_tests()

if ~(snctools_use_tmw || snctools_use_mexnc)
	fprintf ( 1, 'Cannot use native netcdf support or mexnc, no tests requiring netcdf output can be run.\n' );	
	return
end

test_nc_varput           ( 'test.nc' );
test_nc_add_dimension    ( 'test.nc' );
test_nc_addhist          ( 'test.nc' );
test_nc_addvar           ( 'test.nc' );
test_nc_attput           ( 'test.nc' );
test_nc_create_empty     ( 'test.nc' );
test_nc_varrename        ( 'test.nc' );
test_nc_addnewrecs       ( 'test.nc' );
test_nc_add_recs         ( 'test.nc' );
test_nc_archive_buffer   ( 'test.nc' );

test_snc2mat             ( 'test.nc' );
test_nc_getall           ( 'test.nc' );
test_nc_diff             ( 'test1.nc', 'test2.nc' );
test_nc_cat_a;



return

