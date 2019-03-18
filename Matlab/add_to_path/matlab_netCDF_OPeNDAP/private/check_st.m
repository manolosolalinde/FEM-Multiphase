function pos = check_st(name, string, n)

%  CHECK_ST checks if a character string is stored in a cell or character array.
%--------------------------------------------------------------------
%     Copyright (C) J. V. Mansbridge, CSIRO, january 23 1992
%     Revision $Revision: 1.3 $
% CHANGE   1.3 92/04/03
%
%  function pos = check_st(name, string, n)
%
% DESCRIPTION:
%  This function checks whether the character string stored in 'name' is stored
%  in the variable 'string'. 'string' may be an array of cells or a matrix of
%  type char. In the latter case it is assumed that each row in 'string' is
%  completed with blank fill if this is necessary.
% 
% INPUT:
%  name: is a string of characters that we are searching for in the
%       array 'string'.  It does not have blank fill at the end.
%  string: is an array of cells or a matrix where each row is a string which
%         may have blank fill at the end.
%  n: is the number of rows in the array 'string'. 
%
% OUTPUT:
%  pos is the index giving the position of name within the array.
%  pos = -1 if the name is not found within the array.
%
% EXAMPLE:
%  pos = check_st('fred', string, 4)
%  where string = [ 'jim   ' ; 'john  ' ; 'janet ' ; 'fred  ' ]
%
% This function calls: NONE
% This function is called by: getnc.m, getnc_s.m
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

% @(#)check_st.m   1.3   92/04/03
%     Copyright (C), J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: check_st.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

pos = -1;
if iscell(string)
  for ii = 1:n
    if strcmp(name, string{ii})
      pos = ii;
      return
    end
  end
else
  for ii = 1:n
    star = strtok(string(ii, :));
    if strcmp(name, star) == 1
      pos = ii;
      return
    end
  end
end