function [imgwind, outlined_image] = wind_and_outline(img,mask,mask2,lo,hi)
%  set windowing parameters
%
%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

if ~exist('lo','var')
    lo=min(min(img));
end
if ~exist('hi','var')
    hi=max(max(img));
end

%window input image
imgwind=uint8((double(img)-lo).*255./(hi-lo));

%get mask "outline"
%..strip off 10 pixel edge of mask
mask_small=mask;
for i=1:7
    mask_small=imerode(mask_small,strel('arbitrary',[0 1 0; 1 1 1; 0 1 0]));
end
mask_outline=logical(mask-mask_small);

mask2_small=mask2;
for i=1:7
    mask2_small=imerode(mask2_small,strel('arbitrary',[0 1 0; 1 1 1; 0 1 0]));
end
mask2_outline=logical(mask2-mask2_small);

%generate outlined image
outlined_image=imoverlay(imgwind,mask_outline,[0 1 0]);%[1 1 1]);
outlined_image=imoverlay(outlined_image,mask2_outline,[1 0 0]);%[1 1 1]);
