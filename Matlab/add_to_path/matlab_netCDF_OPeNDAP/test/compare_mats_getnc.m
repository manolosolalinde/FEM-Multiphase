function compare_mats_getnc(val, xx, varname, desc)
% val is the array stored in the mat file
% xx is the array as returned by getnc
% varname is the name of the variable that was retrieved with getnc
% desc is a description string to be included in any error message

% $Id: startup-xemacs.el,v 1.2 2003/10/17 04:46:07 man133 Exp man133 $
% Copyright J. V. Mansbridge, CSIRO, Mon Jun 19 15:51:44 EST 2006
  
  if (ndims(xx) ~= ndims(val)) | (any(size(xx) ~= size(val)))
    if ischar(xx)
      val_1d = val';
      val_1d = val_1d(:);
      if size(xx) == size(val_1d)
	d1 = abs(xx - val_1d);
	d1(find(isnan(d1))) = 0;
	dd = d1 + abs(isnan(xx) - isnan(val_1d));
	if max(dd(:)) > 0
	  error(['!!a: ' desc ': getnc failed: got wrong value of ' ...
		 varname]);
	else
	  disp(['!!Character array data is scrambled for ' varname])
	end
      else
	error(['!!b: ' desc ': getnc failed: got wrong size of ' ...
	       varname]);
      end
    else
      if ischar(val) & ispc
	disp('!! Known windows problem with loaddods')
	error(['!! c: ' desc ': getnc got wrong value of character array ' ...
	       varname]);
      else
	error(['!!d: ' desc ': getnc failed: got wrong size of ' varname]);
      end
    end
  else
    if prod(size(xx)) > 0
      d1 = abs(xx - val);
      d1(find(isnan(d1))) = 0;
      dd = d1 + abs(isnan(xx) - isnan(val));
      if max(dd(:)) > 0
	if ischar(val) & ispc
	  disp('!! Known windows problem with loaddods')
	  error(['!! Got wrong value of character array ' var{ii} ':d_nan']);
	else
	  error(['!!e: ' desc ': getnc failed: got wrong value of ' varname]);
	end
      end
    end
  end
