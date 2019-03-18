function retrieval_method = snc_read_backend(ncfile)

switch ( version('-release') )
	case { '11', '12', '13', '14', '2006a', '2006b', '2007a', '2007b', '2008a' }
		tmw_gt_r2008a = false;
	otherwise
		tmw_gt_r2008a = true;
end

file_is_nc3 = exist(ncfile,'file') && isnc3(ncfile);
use_java = getpref('SNCTOOLS','USE_JAVA',false);

if tmw_gt_r2008a && file_is_nc3
	% Use TMW for all local NC3 files when the version >= R2008b
	retrieval_method = 'tmw';
elseif file_is_nc3
	% Local NC3 files should rely on mexnc when the versiohn <= R2008a
	retrieval_method = 'mexnc';
elseif use_java
	% Can be a URL or NC4 file
	retrieval_method = 'java';
else
	% could be a URL where we prefer opendap-enabled mexnc
	% could be an NC4 file where we do not wish to use java
	retrieval_method = 'mexnc';
end

return



