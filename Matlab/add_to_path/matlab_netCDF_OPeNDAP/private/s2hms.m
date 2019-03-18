function [hour,min,sec]=s2hms(secs)
% S2HMS      converts seconds to integer hour,minute,seconds
%
% This function is called by: gregorian.m

%     $Id: s2hms.m Mon, 03 Jul 2006 17:16:40 $

sec=round(secs);
hour=floor(sec/3600);
min=floor(rem(sec,3600)/60);
sec=round(rem(sec,60));
