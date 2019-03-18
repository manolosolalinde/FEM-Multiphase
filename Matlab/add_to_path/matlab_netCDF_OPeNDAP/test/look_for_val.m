function index = look_for_val(desc, field_name_1, field_name_2, val, ...
			      stop_on_error)
  index = [];
  f1 = getfield(desc, field_name_1);
  for ii = 1:length(f1)
    f2 = getfield(f1(ii), field_name_2);
    if strcmp(f2, val)
      index = ii;
      break
    end
  end
  if stop_on_error
    if isempty(index)
      error([field_name_1 '.' field_name_2 ' does not contain ' val])
    end
  end
