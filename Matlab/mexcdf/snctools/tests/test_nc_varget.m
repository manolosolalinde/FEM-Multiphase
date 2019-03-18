function test_nc_varget( )


fprintf ( 1, 'NC_VARGET:  starting test suite...\n' );

%create_test_file ( ncfile );


% Partial Retrievals
test_readSingleValueFrom1dVariable ( 'testdata/varget.nc' );
test_readSingleValueFrom2dVariable ( 'testdata/varget.nc' );
test_read2x2hyperslabFrom2dVariable ( 'testdata/varget.nc' );
test_readSingleValueFromNc4File ( 'testdata/tst_pres_temp_4D_netcdf4.nc' );

% Full retrievals
test_readFullSingletonVariable ( 'testdata/varget.nc' );
test_readFullDoublePrecisionVariable ( 'testdata/varget.nc' );
test_readFullFloatVariableNc4 ( 'testdata/tst_pres_temp_4D_netcdf4.nc' );

% Strided retrievals
test_readStridedVariable ( 'testdata/varget.nc' );

% OPeNDAP
test_readOpendapVariable;

% Regression
test_regressionErrorMsgBadUrl;
test_regressionNegSize('testdata/varget.nc');
test_scaling('testdata/varget.nc');

return





function test_readOpendapVariable ()
if getpref('SNCTOOLS','TEST_OPENDAP',false) && getpref('SNCTOOLS','TEST_OPENDAP',false)
    url = 'http://motherlode.ucar.edu:8080/thredds/dodsC/nexrad/composite/1km/agg';
    fprintf ( 1, 'Testing remote URL access %s...\n', url );
    w = nc_varget ( url, 'y', [0], [1] );
else
	fprintf('Reading OPeNDAP data not tested.  Read the README for details.\n');
end
return



%==============================================================================
function test_regressionErrorMsgBadUrl ()
% Regression test.  If the URL is wrong, then the error message must give the
% name of the wrong url.   01-04-2007
% 
if snctools_use_java
    url = 'http://doesntexits:8080/thredds/dodsC/nexrad/composite/1km/agg';
    try
        w = nc_varget ( url, 'y', [0], [1] );
    catch
        [msg,id] = lasterr;
        if ~strcmp(id, 'SNCTOOLS:nc_varget_java:fileOpenFailure')
            error ( 'Error id ''%s'' was not expected.', id );
        end
        if ~findstr(msg, url)
            error ( 'Error message did not contain the incorrect url.');
        end
    end
end
return







function test_readSingleValueFrom1dVariable ( ncfile )

expData = 1.2;
actData = nc_varget ( ncfile, 'test_1D', 1, 1 );

ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    msg = sprintf ( 'input data ~= output data.' );
    error ( msg );
end

return








function test_readSingleValueFrom2dVariable ( ncfile )

expData = [1.5];
actData = nc_varget ( ncfile, 'test_2D', [2 2], [1 1] );

ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    msg = sprintf ( '%s:  input data ~= output data.\n', mfilename );
    error ( msg );
end

return




function test_read2x2hyperslabFrom2dVariable ( ncfile )

expData = [1.5 2.1; 1.6 2.2];
if getpref('SNCTOOLS','PRESERVE_FVD',false)
    expData = expData';
end
actData = nc_varget ( ncfile, 'test_2D', [2 2], [2 2] );

if ndims(actData) ~= 2
    error ( 'rank of output data was not correct' );
end
if numel(actData) ~= 4
    error ( 'rank of output data was not correct' );
end
ddiff = abs(expData(:) - actData(:));
if any( find(ddiff > eps) )
    error ( 'input data ~= output data ' );
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



%--------------------------------------------------------------------
function test_readSingleValueFromNc4File ( ncfile )

switch ( version('-release') )
    case { 'R14', '2006a', '2006b', '2007a', '2007b', '2008a' }
        try 
            actData = nc_varget ( ncfile, 'latitude', 1, 1 );
        catch
            [emsg,eid] = lasterr;
			switch(eid)
			case 'MATLAB:Java:GenericException'
                fprintf ( 1, 'Could not read a NC4 file with java, you need >= version 4.0 of toolsUI installed.\n' );
                return
			case 'SNCTOOLS:NC_VARGET:MEXNC:OPEN'
				fprintf ( 1, '\n\n\n' );
                fprintf ( 1, 'Could not read a NC4 file with mexnc, you would need to compile the netcdf library version >= version 4.0 (no, I won''t do that for you).\n' );
				fprintf ( 1, '\n\n\n' );
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





function test_readFullSingletonVariable ( ncfile )


expData = 3.14159;
actData = nc_varget ( ncfile, 'test_singleton' );

ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    error ( 'input data ~= output data.\n'  );
end

return



function test_readFullDoublePrecisionVariable ( ncfile )


expData = [1:24];
expData = reshape(expData,6,4) / 10;

if getpref('SNCTOOLS','PRESERVE_FVD',false)
    expData = expData';
end

actData = nc_varget ( ncfile, 'test_2D' );

ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    error ( 'input data ~= output data.\n'  );
end

return




%---------------------------------------------------------------------------
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



function test_readFullFloatVariableNc4 ( ncfile )

if snctools_use_java
    switch ( version('-release') )
        case { 'R14', '2006a', '2006b', '2007a', '2007b', '2008a' }
            try 
                actData = nc_varget ( ncfile, 'latitude' );
            catch
                [emsg,eid] = lasterr;
                if strcmp(eid,'MATLAB:Java:GenericException')
                    warning ('Could not read an NC4 file, make sure you have >= version 4.0 of toolsUI installed.\n' );
                    return
                else
                    error(eid,emsg);
                end


            end
        otherwise % assume >= R2008b
			actData = test_readFullFloatVariableNc4_tmw ( ncfile )
    end
    
    test_readFullFloatVariable_backend_neutral (actData);
    
end
return




function test_readStridedVariable ( ncfile )

expData = [1:24];
expData = reshape(expData,6,4) / 10;
expData = expData(1:2:3,1:2:3);
if getpref('SNCTOOLS','PRESERVE_FVD',false)
    expData = expData';
end

actData = nc_varget ( ncfile, 'test_2D', [0 0], [2 2], [2 2] );

ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    error ( 'input data ~= output data.\n'  );
end

return



function test_regressionNegSize ( ncfile )

expData = [1:24];
expData = reshape(expData,6,4) / 10;
sz = size(expData);
sz(2) = -1;
if getpref('SNCTOOLS','PRESERVE_FVD',false)
    expData = expData';
	sz = fliplr(sz);
end

actData = nc_varget ( ncfile, 'test_2D', [0 0], sz );

ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    error ( 'input data ~= output data.\n'  );
end

return

function test_scaling ( ncfile )

expData = [32 32 32 32; 50 50 50 50; 68 68 68 68; ...
           86 86 86 86; 104 104 104 104; 122 122 122 122];
expData = reshape(expData,6,4);
sz = size(expData);
sz(2) = -1;
if getpref('SNCTOOLS','PRESERVE_FVD',false)
    expData = expData';
	sz = fliplr(sz);
end

actData = nc_varget ( ncfile, 'temp' );

if ~isa(actData,'double')
    error ( 'short data was not converted to double');
end
ddiff = abs(expData - actData);
if any( find(ddiff > eps) )
    error ( 'input data ~= output data.\n'  );
end

return

