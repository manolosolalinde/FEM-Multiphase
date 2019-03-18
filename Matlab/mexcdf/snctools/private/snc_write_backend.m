function backend = snc_write_backend(ncfile)

switch ( version('-release') )
	case { '11', '12', '13', '14', '2006a', '2006b', '2007a', '2007b', '2008a' }
		tmw_gt_r2008a = false;
	otherwise
		tmw_gt_r2008a = true;
end

file_is_nc3 = exist(ncfile,'file') && isnc3(ncfile);

if tmw_gt_r2008a && file_is_nc3
	% Use TMW for all local NC3 files when the version >= R2008b
	backend = 'tmw';
else
	% could be an NC4 file
	% could be a release less than R2008b
	backend = 'mexnc';
end

return




