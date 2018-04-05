function version = libra_version()
%  version = libra_version;
%  To get the versioning of the loaded LIBRA package.
%
%  Version info:
%  $Rev: 603 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2016-10-20 16:58:24 -0400 (Thu, 20 Oct 2016) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>


revision = '$Revision: 603 $';
revision = strtrim(strrep(revision,'$',''));
version = ['1.0.4 (' revision ')'];
