function [roih] = croi(I0, bdrind,coef)
%  Choose Region of Interest for Hough Transform
%  roih returns the region of interest for Hough transformation
%
%  Version info:
%  $Rev: 603 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2016-10-20 16:58:24 -0400 (Thu, 20 Oct 2016) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

[dimx, ~] = size(I0);
rub = round(bdrind(find(bdrind(10:end)>=3,1)+9)*coef);
lb = min(find(bdrind==max(bdrind(round(dimx*0.25):end)), 1, 'last'),round(dimx*0.50)); %skip the first few rows just in case 
roih = I0(1:lb,1:rub);
clear lb rub dimx dimy
