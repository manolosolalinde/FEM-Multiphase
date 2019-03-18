function path_name = pos_cds() 

% POS_CDS returns the path to the common data set directory.
%--------------------------------------------------------------------
%     Copyright (C) J. V. Mansbridge, CSIRO, april 15 1992
%     Revision $Revision: 1.3 $
%
% DESCRIPTION:
% pos_cds returns a string containing the path to the common data set directory.
% This is the directory containing netcdf files accessible to all users. THIS M
% FILE SHOULD BE EDITED TO SUIT THE LOCAL INSTALLATION.
% 
% INPUT:
% none
%
% OUTPUT:
% path_name: a string containing the path to the common data set directory.
%
% EXAMPLE:
% Simply type path_name at the matlab prompt.
%
% This function is called by: attnc.m, check_nc.m, getnc.m, getnc_s.m, inqnc.m
%                             timenc.m, whatnc.m
%
% AUTHOR:   J. V. Mansbridge, CSIRO
%---------------------------------------------------------------------

%     Copyright (C), J.V. Mansbridge, 
%     Commonwealth Scientific and Industrial Research Organisation
%     $Id: pos_cds.m Mon, 03 Jul 2006 17:16:40 $
% 
%--------------------------------------------------------------------

path_name = [ ];

