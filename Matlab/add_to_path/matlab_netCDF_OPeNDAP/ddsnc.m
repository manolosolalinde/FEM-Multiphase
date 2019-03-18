function desc = ddsnc(file)
% DDSNC returns information about a netcdf file or DODS/OPEnDAP dataset
%
%    function desc = ddsnc(file)
%
% DESCRIPTION:
%  ddsnc is a non-interactive function that returns information about a
%  netcdf file or DODS dataset.
% INPUT:
%  file may be the name of a netCDF file with or it may be the URL of a
%  DODS/OPEnDAP dataset. 
%
% OUTPUT:
%  desc is a matlab structure. For an OPEnDAP data set desc will contain all
%  of the information in the DDS (Dataset Description Structure). For a
%  netCDF file desc will be almost identical.
%
% EXAMPLE:
% Information about the Reynolds data set can be found as follows:
% >> file = 'http://www.marine.csiro.au/dods/nph-dods/dods-data/climatology-netcdf/sst.wkmean.1981-1989.nc';
% >> desc = ddsnc(file)
%
% desc has 2 fields - variable and dimension. Looking at one element we see
% >> desc.variable(4) =
%
%              type: 'Int16'
%              name: 'sst'
%     dim_statement: {'time = 427'  'lat = 180'  'lon = 360'}
%        dim_idents: [3x1 double]
%
% The first 2 fields tell us that the variable is named 'sst' and is a 2 byte
% integer. The dim_statement field tells us that the sst variable has 3
% dimensions in the order given. For dim_idents we see
% >> desc.variable(4).dim_idents = 
% >>            [3 1 2]'.
%
% These integers refer to the dimensions of the sst array. Looking at
% >> desc.dimension(3) =
%       name: 'time'
%     length: 427
%
% we see that index 3 points us to the 3rd dimension, time and it has length
% 427. A generic program could then retrieve the information by setting:
% >> ii = desc.variable(4).dim_idents(1);
%
% and then referring to desc.dimension(ii).
%
% Notes
% 1) Different functions are called for netcdf and opendap files. Also the
% netcdf data model does not map exactly onto the opendap model. As a result the
% information in desc may be arranged differently according to whether a file is
% accessed as a local netcdf file or an opendap data set. Specifically, the
% dimensions may be in a different order in desc and so the values in dim_idents
% will be different.
% 2) If file is a URL then desc will simply be the structure returned by
% get_dods_dds.

% $Id: ddsnc.m Mon, 03 Jul 2006 17:16:40 $
% Copyright J. V. Mansbridge, CSIRO, Wed Dec 14 15:02:57 EST 2005

% Check the input arguments

if nargin < 1
   help ddsnc
   return
end

% Initialisation

CSIRO_add_jar_file_maybe;
[mex_name, full_name, desc_das, file_status, exe_name] = ...
    choose_mexnc_opendap(file);

switch mex_name
 case 'mexnc'
  
  % I make mexnc calls to find the integers that specify the attribute
  % types

  nc_byte = mexnc('parameter', 'nc_byte'); %1
  nc_char = mexnc('parameter', 'nc_char'); %2
  nc_short = mexnc('parameter', 'nc_short'); %3
  nc_long = mexnc('parameter', 'nc_long'); %4
  nc_float = mexnc('parameter', 'nc_float'); %5
  nc_double = mexnc('parameter', 'nc_double'); %6

  % Store the type information in a form that refers to the dods data types.
  
  data_type{nc_byte} = 'Byte';
  data_type{nc_char} = 'String';
  data_type{nc_short} = 'Int16';
  data_type{nc_long} = 'Int32';
  data_type{nc_float} = 'Float32';
  data_type{nc_double} = 'Float64';
  
  % Open the netcdf file.
  
  [cdfid, rcode] = mexnc('ncopen', full_name, 'nowrite');
  if rcode < 0
    error(['mexnc: ncopen: rcode = ' int2str(rcode)])
  end

  % Don't print out netcdf warning messages

  mexnc('setopts', 0);

  % Collect information about the cdf file.

  [num_dims, nvars, ngatts, recdim, rcode] =  mexnc('ncinquire', cdfid);
  if rcode < 0
    error([ 'mexnc: ncinquire: rcode = ' int2str(rcode) ])
  end
  
  % Initialise
  
  desc.variable = [];
  desc.dimension = [];
  
  % Step through each dimension adding the appropriate parts to the desc
  % structure
  
  for ii = 1:num_dims
    [dimnam, dimsiz, rcode] = mexnc('ncdiminq', cdfid, (ii - 1));
    desc.dimension(ii).name = dimnam;
    desc.dimension(ii).length = dimsiz;
  end

  % Step through each variable adding the appropriate parts to the desc
  % structure
  
  for ii = 1:nvars
    [varnam, vartyp, nvdims, vdims, nvatts, rcode] = mexnc('ncvarinq', ...
						  cdfid, (ii - 1));
    desc.variable(ii).type = data_type{vartyp};
    desc.variable(ii).name = varnam;
    if nvdims == 0
      desc.variable(ii).dim_statement = [];
    else
      for jj = 1:nvdims
	str = [desc.dimension(vdims(jj) + 1).name ' = ' ...
	       num2str(desc.dimension(vdims(jj) + 1).length)];
	desc.variable(ii).dim_statement{jj} = str;
      end
    end
    desc.variable(ii).dim_idents = vdims(:) + 1;
  end
  
  % Close the netcdf file.
  
  [rcode] = mexnc('ncclose', cdfid);
  if rcode < 0
    error(['** ERROR ** ncclose: rcode = ' num2str(rcode)])
  end
  
 case {'loaddap', 'loaddods'}
  
  [dds_text, desc] = get_dods_dds(full_name, exe_name);
  %[desc_text, status] = urlread([full_name '.dds']);
  %if ~status
  %  error(['couldn''t get the dds of ' full_name])
  %end
 case 'java'
  switch file_status
   case 4
    [dds_text, status] = urlread([full_name '.dds']);
    if ~status
      error(['couldn''t get the dds of ' file])
    end
    desc = CSIRO_organise_dds_output(dds_text);
   otherwise
    try
      ncdJ = ucar.nc2.dataset.NetcdfDataset.openDataset(full_name);
    catch
      error(['The file ' file ' is not accessible'])
    end

    % Initialise
  
    desc.variable = [];
    desc.dimension = [];
  
    % Step through each dimension adding the appropriate parts to the desc
    % structure
  
    dimensions_list = ncdJ.getDimensions();
    num_dims = dimensions_list.size();
    for ii = 1:num_dims
      dim = dimensions_list.get(ii - 1);
      desc.dimension(ii).name = char(dim.getName());
      desc.dimension(ii).length = dim.getLength();
    end
    
    % Step through each variable adding the appropriate parts to the desc
    % structure
  
    variables_list = ncdJ.getVariables();
    nvars = variables_list.size();
    for ii = 1:nvars
      var = variables_list.get(ii - 1);
      desc.variable(ii).type = char(var.getDataType());
      desc.variable(ii).name = char(var.getName());
      dim_list = var.getDimensions;
      nvdims = dim_list.size();
      vdims = zeros(1, nvdims);
      if nvdims == 0
	desc.variable(ii).dim_statement = [];
      else
	for jj = 1:nvdims
	  dim_val = dim_list.get(jj - 1);
	  name_comp = char(dim_val.getName);
	  for kk = 1:num_dims
	    if strcmp(desc.dimension(kk).name, name_comp)
	      vdims(jj) = kk;
	      break
	    end
	  end
	end
	
	for jj = 1:nvdims
	  str = [desc.dimension(vdims(jj)).name ' = ' ...
		 num2str(desc.dimension(vdims(jj)).length)];
	  desc.variable(ii).dim_statement{jj} = str;
	end
      end
      desc.variable(ii).dim_idents = vdims(:);
    end
    ncdJ.close();
  end
end


