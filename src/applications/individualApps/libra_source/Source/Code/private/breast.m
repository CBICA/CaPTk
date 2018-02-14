function [mask] =  breast(I0,peind,bdrind,airthresh)
%  store the breast region
%
%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

%  muscle information/create a new grey level information
[dimx,dimy] = size(I0);
mask = zeros(dimx, dimy);
mini=airthresh;
for i = 1:dimx
    for j = 1:dimy
         mask(i,j) = (j > peind(i)&& j <= bdrind(i) && I0(i,j)>mini);
    end
end

tmp=bwlabel(mask);tmp(tmp==0)=NaN;
mask=double(tmp==mode(tmp(:)));

%now lets fill in any "islands" in the breast mask, can happen in processed
%mammograms where processing makes the intensity of breast pixels<=air
%pixels
mask_neg=~mask;
[L,NUM] = bwlabeln(mask_neg,8);

%Any 'background' region NOT adjacent to the image border needs to be
%filled in
%First, find labels of border regions
Lu=unique([L(1,:),L(:,1)',L(end,:),L(:,end)']);
%next, make a map of where to fill in the mask
fill_in=zeros(size(mask));
for i=1:NUM
    if ~max(Lu==i)
        fill_in(L==i)=1;
    end
end

%now fill the holes (if its 1 in either mask or fill_in, it's breast area)
mask=mask+fill_in;
     
