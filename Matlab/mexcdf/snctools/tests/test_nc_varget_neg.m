function test_nc_varget_neg()

testroot = fileparts(mfilename('fullpath'));

run_opendap_tests;
run_nc4_java_tests(testroot);

return


%--------------------------------------------------------------------------
function run_nc4_java_tests(testroot)

	if ~getpref('SNCTOOLS','USE_JAVA',false)
		return
	end

	ncfile = fullfile(testroot,'testdata/tst_pres_temp_4D_netcdf4.nc');
	% Partial Retrievals
	test_readSingleValueFromNc4File(ncfile );
	test_readFullFloatVariableNc4(ncfile);
    

return

%--------------------------------------------------------------------------
function run_opendap_tests()
if getpref('SNCTOOLS','TEST_OPENDAP',false)

	% Regression
	test_regressionErrorMsgBadUrl;

end

return

%==========================================================================
function test_regressionErrorMsgBadUrl ()
% Regression test.  If the URL is wrong, then the error message must give 
% name of the wrong url.   01-04-2007
% 

    url = 'http://doesntexits:8080/thredds/dodsC/nexrad/composite/1km/agg';
    try
        nc_varget ( url, 'y', 0, 1 );
    catch me
        if ~strcmp(me.identifier,'SNCTOOLS:nc_varget_java:fileOpenFailure')
            error ( 'Error id ''%s'' was not expected.', id );
        end
        if ~findstr(me.message, url)
            error ( 'Error message did not contain the incorrect url.');
        end
        fprintf('\n\n\tThe above error message is expected, don''t freak out...\n\n');
    end
return








%--------------------------------------------------------------------
function test_readSingleValueFromNc4File ( ncfile )

switch ( version('-release') )
    case { '2007b', '2008a' }
        try 
            actData = nc_varget ( ncfile, 'latitude', 1, 1 );
        catch me
            
            switch(me.identifier)
                case 'MATLAB:Java:GenericException'
                    fprintf('Could not read a NC4 file with java, you ');
                    fprintf('need >= version 4.0 of toolsUI installed.\n');
                    return
                case 'SNCTOOLS:NC_VARGET:MEXNC:OPEN'
                    fprintf('\n\n\n' );
                    fprintf('Could not read a NC4 file with mexnc, you ');
                    fprintf('would need to compile the netcdf library ');
                    fprintf('version >= version 4.0 (no, I won''t do ');
                    fprintf('that for you).\n' );
                    fprintf('\n\n\n' );
                    pause(3);
                    return
                otherwise
                    error(eid,emsg);
            end


        end
    otherwise % assume >= R2008b
	    actData = test_readSingleValueFromNc4File_tmw(ncfile);
end

if ~isempty(actData) 
	test_readSingleValueFromNc4File_backend_neutral(actData);
end
    
return








%--------------------------------------------------------------------------
function test_readSingleValueFromNc4File_backend_neutral ( actData )

expData = 30;
if ndims(actData) ~= 2
    error ( 'rank of output data was not correct' );
end
if numel(actData) ~= 1
    error ( 'size of output data was not correct' );
end
ddiff = abs(expData(:) - actData(:));
if any( find(ddiff > eps) )
    error ( 'input data ~= output data ' );
end

return



%--------------------------------------------------------------------------
function test_readFullFloatVariableNc4 ( ncfile )

if snctools_use_java
    switch ( version('-release') )
        case { '2007b', '2008a' }
            try 
                actData = nc_varget ( ncfile, 'latitude' );
            catch me
 
                if strcmp(me.identifier,'MATLAB:Java:GenericException')
                    fprintf('\tCould not read an NC4 file, make sure ');
                    fprintf('you have >= version 4.0 of toolsUI ');
                    fprintf('installed.' );
                    return
                else
                    rethrow(me);
                end


            end
        otherwise % assume >= R2008b
			actData = test_readFullFloatVariableNc4_tmw ( ncfile );
    end
    
    test_readFullFloatVariable_backend_neutral (actData);
    
end
return






%--------------------------------------------------------------------------
function test_readFullFloatVariable_backend_neutral ( actData )

expData = [25 30 35 40 45 50];

if ndims(actData) ~= 2
    error ( 'rank of output data was not correct' );
end
if numel(actData) ~= 6
    error ( 'rank of output data was not correct' );
end
ddiff = abs(expData(:) - actData(:));
if any( find(ddiff > eps) )
    error ( 'input data ~= output data ' );
end

return




