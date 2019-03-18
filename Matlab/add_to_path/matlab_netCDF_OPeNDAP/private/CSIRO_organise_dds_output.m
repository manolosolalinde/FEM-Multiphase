function desc = CSIRO_organise_dds_output(dds_text)
  
% Parse the text version, building the structure desc as we go.

% Divide text string into lines.

ff_eol = findstr(dds_text, char(10));
ff_start = [1 (ff_eol(1: (end-1)) + 1)];
ff_end = ff_eol - 1;

% Step through one line at a time (ignoring first and last lines).

in_grid = 0;
in_array = 0;
num_var = 0;
num_dim = 0;
for ii = 2:(length(ff_start) - 1)
  tt = dds_text(ff_start(ii): ff_end(ii));
  
  % Find out what type of line will be next. We only need to read one type of
  % line - the one that describes a variable, e.g.,
  %        Int16 sst[time = 427][lat = 180][lon = 360];
  
  if ~isempty(findstr(tt, '{'))
    in_grid = 1;
    continue
  elseif ~isempty(findstr(tt, '}'))
    in_grid = 0;
    continue
  elseif ~isempty(findstr(lower(tt), 'array:'))
    in_array = 1;
    continue
  elseif ~isempty(findstr(lower(tt), 'maps:'))
    in_array = 0;
    continue
  end
  
  % We have found the single type of line that we need to parse.
  
  if (in_grid & in_array) | ~in_grid
    
    % Get the variable name and type from the start of the line.
    
    num_var = num_var + 1;
    [t, r] = strtok(tt);
    desc.variable(num_var).type = t;
    [t, r] = strtok(r, '[');
    %desc.variable(num_var).name = strtrim(t);
    desc.variable(num_var).name = deblank(fliplr(deblank(fliplr(t))));
    
    % Get the dimension information from the latter part of the line. This is
    % messy because we have to figure out whether we have met this dimension
    % name previously.
    
    num_steps = length(findstr(r, '['));
    dim_idents = zeros(num_steps, 1);
    dim_st = {};
    for jj = 1:num_steps
      
      % Find the dimension name.
      
      [t, r] = strtok(r, ']');
      t_part = t(2:end);
      dim_st{jj} = t_part;
      [dim_name, remain] = strtok(t_part);
      
      % Find whether we have met this dimension name previously.
      
      [junk, remain] = strtok(remain);
      dim_length_str = strtok(remain);
      is_new_dim = 1;
      for kk = 1:num_dim
	if strcmp(dim_name, desc.dimension(kk).name)
	  is_new_dim = 0;
	  index_dim = kk;
	  break
	end
      end
      
      % If this is a new dimension name store information about it in desc.
      
      if is_new_dim
	num_dim = num_dim + 1;
	desc.dimension(num_dim).name = dim_name;
	desc.dimension(num_dim).length = str2num(dim_length_str);
	index_dim = num_dim;
      end
      dim_idents(jj) = index_dim;
    end
    
    % Store some dimension information in the desc.variable part of the
    % structure so that we can identify which dimensions are associated with
    % the given variable.
    
    desc.variable(num_var).dim_statement = dim_st;
    desc.variable(num_var).dim_idents = dim_idents;
  end
end

