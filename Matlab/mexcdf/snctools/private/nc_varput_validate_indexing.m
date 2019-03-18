function [start, count] = nc_varget_validate_indexing(ncid,nvdims,data,start,count,stride)
% Check that any given start, count, and stride arguments actually make sense
% for this variable.  


%
% Singletons are a special case.  We need to set the start and count carefully.
if nvdims == 0

    if isempty(start) && isempty(count) && isempty(stride)

        %
        % This is the case of "nc_varput ( file, var, single_datum );"
        start = 0;
        count = 1;

    elseif ~isempty(start) && ~isempty(count) && ~isempty(stride)
        mexnc ( 'close', ncid );
        err_id = 'SNCTOOLS:NC_VARPUT:MEXNC:badIndexing';
        err_msg = 'Strides make no sense for a singleton variable.';
        error ( err_id, err_msg );
    end

    return;

end

% If START and COUNT not given, and if not a singleton variable, then START is [0,..] and COUNT is 
% the size of the data.  
if isempty(start) && isempty(count) && ( nvdims > 0 )
    start = zeros(1,nvdims);
    count = zeros(1,nvdims);
    for j = 1:nvdims
        count(j) = size(data,j);
    end
end


%
% Check that the start, count, and stride arguments have the same length.
if ( numel(start) ~= numel(count) )
    mexnc ( 'close', ncid );
    err_id = 'SNCTOOLS:NC_VARPUT:MEXNC:badIndexing';
    err_msg = 'START and COUNT arguments must have the same length.';
    error ( err_id, err_msg );
end
if ( ~isempty(stride) && (length(start) ~= length(stride)) )
    mexnc ( 'close', ncid );
    err_id = 'SNCTOOLS:NC_VARPUT:MEXNC:badIndexing';
    err_msg = 'START, COUNT, and STRIDE arguments must have the same length.';
    error ( err_id, err_msg );
end






