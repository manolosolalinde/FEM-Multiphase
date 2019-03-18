function [dds_text, desc] = get_dods_dds(url, exe_name)
% get_dods_dds finds the Dataset Descriptor Structure (DDS) for a DODS data set
%
%  function [dds_text, desc] = get_dods_dds(url, exe_name)
%
%     INPUT:
% url: A DODS URL
% exe_name:
%   the name of the executable file associated with the mexfile. It contains
%   the full path name so that it can be called from other functions (notably
%   get_dods_dds). It may be writedap, writeval, writedap.exe or
%   writeval.exe. If the mexfile is mexnc then there is no associated
%   executable and so exe_name is 'none'.
%
%     OUTPUT:
% dds_text: A text version of the DODS Dataset descriptor structure. This is
%            the same thing that you would get by putting the .dds extension
%            on the url and feeding it into a web browser.
% desc: A matlab structure containing the same information as dds_text (in
%       fact it is derived from it). Being in a matlab structure it is then
%       handy for using in batch matlab code.
%
% EXAMPLE:
% The Reynolds data set can be found with the following url:
% url = 'http://www.marine.csiro.au/dods/nph-dods/dods-data/climatology-netcdf/sst.wkmean.1981-1989.nc';
%
% Using get_dods_dds gives:
% dds_text =
% Dataset {
%     Float32 lat[lat = 180];
%     Float32 lon[lon = 360];
%     Float64 time[time = 427];
%     Grid {
%      ARRAY:
%         Int16 sst[time = 427][lat = 180][lon = 360];
%      MAPS:
%         Float64 time[time = 427];
%         Float32 lat[lat = 180];
%         Float32 lon[lon = 360];
%     } sst;
%     Grid {
%      ARRAY:
%         Int16 mask[lat = 180][lon = 360];
%      MAPS:
%         Float32 lat[lat = 180];
%         Float32 lon[lon = 360];
%     } mask;
% } sst.wkmean.1981-1989.nc;
%
% This shows that there are 5 variables and 3 dimensions. The names lat, lon
% and time are used for both variables and dimensions. The lengths are 180,
% 360 and 427 repectively. The variable sst is 3 dimensional with the order
% of the dimensions as shown while the variable mask is 2 dimensional.
%
% desc has 2 fields - variable and dimension. Looking at one element we see
% >> desc.variable(4) =
%
%              type: 'Int16'
%              name: 'sst'
%     dim_statement: {'time = 427'  'lat = 180'  'lon = 360'}
%        dim_idents: [3x1 double]
%
% The meaning of the first 3 fields is obvious and
% desc.variable(4).dim_idents = [3 1 2]'. These integers refer to the
% dimensions of the sst array. Looking at
% >> desc.dimension(3) =
%       name: 'time'
%     length: 427
% we see that index 3 points us to the 3rd dimension, time and it has length
% 427.
%
% Copyright J. V. Mansbridge, CSIRO, Wed May 11 15:39:34 EST 2005

% This function calls: loaddap
% This function is called by: getnc_s.m, inqnc.m

% $Id: get_dods_dds.m Mon, 03 Jul 2006 17:16:40 $
% Copyright J. V. Mansbridge, CSIRO, Wed May 11 15:39:34 EST 2005

% Get the DDS in text form. Note that in the call to system we send double
% quote marks around the name of the executable. This is because it is
% possible under windows to have spaces in the directory names and we do not
% want these to be missinterpreted.

cmd = ['"' exe_name '" -D -- ', url];
[status, dds_text] = system(cmd);
if status == 0
  % It sometimes happens that when shelling out to the system some warning
  % messages are returned ahead the DDS. As a dangerous workaround to this we
  % get rid of all of the lines preceeding the line beginning "Dataset". This
  % is a horrible kluge but seems to be necessary in some cases.
  
  ff = findstr(dds_text, 'Dataset');
  if isempty(ff)
    error('Dodgy DDS returned by writedap')
  else
    if ff(1) > 1
      dds_text = dds_text(ff(1):end);
    end
  end
else
  error(['Couldn''t find the DDS for ' url])
end

desc = CSIRO_organise_dds_output(dds_text);

