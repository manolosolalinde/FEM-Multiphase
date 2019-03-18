function y = return_v(s, x)

% RETURN_V prints out a string s to prompt the user to make a keyboard input
%--------------------------------------------------------------------
%     Copyright (C) J. V. Mansbridge, CSIRO, 
%     Revision $Revision: 1.4 $
% CHANGE   1.1 92/01/30
%
% function y = return_v(s, x)
%
% DESCRIPTION:
%
% return_v prints out the string s to prompt the user to make a
% keyboard input.  The type of x determines whether the input will be
% interpreted as a string or a matrix of numbers.  This input is
% returned by return_v unless the user hits the carriage return in
% which case the default (x) is returned.
% 
% INPUT:
% s is a string that is written to the user's terminal as a prompt.
% x is the default string or matrix of numbers that is returned if the
% user simply hits the carriage return.
%
% OUTPUT:
% y is the returned value.  It is either the input value or the
% default (x).
%
% EXAMPLE:
% y = return_v('Type in the ocean depth in metres [6000] ...', 6000);
%
% This function is called by: getnc.m
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

%     Copyright (C), J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: return_v.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

if ischar(x)
  xtemp = input(s, 's');
else
  xtemp = input(s);
end
if isempty(xtemp)
   y = x;
else
   y = xtemp;
end
