function [bdrind, airthresh] = con_breast(I0,tmp)
%  Find the contour of breast.
%  bdrind  returns the axis index of the boundary point
%
%  Version info:
%  $Rev: 603 $:  Revision of last commit
%  $LastChangedBy: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2016-10-20 16:58:24 -0400 (Thu, 20 Oct 2016) $:  Date of last commit
%
%  Developed by Dr. B.M. Keller
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

%First remove pixel rows from consideration of determining the threshold if
%the row is the basically the exact the same value
%I call them nonsense bars, because they can't possibly involve a breast/air
%interface, and we only want to use "central pixels"
row_range=max(I0,[],2)-min(I0,[],2);
row_range=row_range./max(row_range);
C1=find(row_range>0.001,1,'first');
C2=find(row_range>0.001,1,'last');
I0_center=tmp(C1:C2,:);

%I0 is the full image, tmp is the roi to evaluate airthresh in

[dimr] = size(I0,1);
%boundary between air and fat
bdrind = zeros(1,dimr);
%find "air-threshold"

img_pixs=double(I0_center(:));
x=min(img_pixs):(max(img_pixs)-min(img_pixs))/1000:max(img_pixs); % central location of bins in the histogram
n_elements = histc(img_pixs,x);
c_elements = [0 cumsum(n_elements)'];
dd=diff(c_elements,1);
dd=conv(dd,gausswin(25),'same');
peaklocation=find(dd>max(dd)*0.05,1);% note: old value for GE in med phys = 0.2 ie 20%

%old way - medphys style
%airthresh=x(peaklocation+find(dd(peaklocation:end)<=(dd(dd==max(dd))*0.05),1)+1);%med phys = 0.01 ie 1% Raised the threshold to 0.04 in response to the smoothing. --> 0.05 for better performance.

%new way - look at first - to + 0-crossing of empirical first derivative
%after the air peak - and have a max increased distance from peak so it can
%way over shoot
ddd=dd(2:end)-dd(1:end-1);
ddd_neg=peaklocation+find(ddd(peaklocation+1:end)<=0,1);
airthresh_opt1=ddd_neg+find(ddd(ddd_neg+1:end)>=0,1)-4; %-4 so as to be sure on the negative side of the 0 crossing, and in the smoother area of the breast edge
airthresh=min([x(airthresh_opt1),x(peaklocation+50)]); %the air peak is never THAT wide



%BMK-Aug28-2015: Deal with spacers/paddles - there should be air somewhere
%on the top or bottom row, if not (either a paddle or super huge breast),
%thus the threshold is too low, and we need to correct
%up it - only go up once though, should return a failure afterwards and
%skip the analysis of image - fix later
I0mask=I0>=airthresh;
I0mask_neg_cc = bwlabel(~I0mask,8);
%choose largest air region - %BMK Nov 2015 - this should be in terms of
%Y-extent, not in general
stats = regionprops(~I0mask,'BoundingBox'); %get y-extent of all regions
s=1;ys=stats(1).BoundingBox(end);
for i=2:size(stats,1)
    if stats(i).BoundingBox(end)>ys
        s=i;
        ys=stats(i).BoundingBox(end);
    end
end

airregion=I0mask_neg_cc==s;
r1=find(max(airregion,[],2), 1 );
r2=find(max(airregion,[],2), 1, 'last' );

%no air on either top or bottom row? fix (ie. choose next best option)
if r1~=1 && r2~=size(I0mask,1)
    disp('Bad mask. Might be presence of spacers or paddles. Redo.')  %%% TEMPORARY FLAG to find out bad cases.
    offset=peaklocation+find(dd(peaklocation:end)<=(dd(dd==max(dd))*0.05),1)+1;
    dd(1:offset)=[];
    peaklocation=find(dd>max(dd)*0.1,1);
    airthresh=x(offset+peaklocation+find(dd(peaklocation:end)<=(dd(dd==max(dd))*0.05),1)+1);
    I0mask=I0>=airthresh;
end

%BMK-July24-2013: next, threshold the air region(s) away, and take the 
%largest contiguous region juxtaposed to left-image-edge (from
%alignment). The breast is one object, it will be in one piece, 
%along the left-image-edge. In case there are columns of 'zeros' along the
%left edge, like as with the rows above, we wont start on column 1, per-se,
%but on the first column with values

I0mask_cc = bwlabel(I0mask,8);
col_check = max(I0mask_cc,[],1);
init_col = find(col_check>0,1); %first non-zero label on left
if isempty(init_col) % not likely as that would mean there is nothing brighter than air
   init_col=1; %but JUST IN CASE, to avoid crashing out
end
%...now, on the off chance there are multiple objects sticking out of the
%left-image-edge (i.e., like a tag AND the breast) take the largest object
labels_2_check=unique(I0mask_cc(:,init_col)); %always at least 2: 0 and label
labels_2_check(labels_2_check==0)=[]; %drop 0

%if there is more than 1 non-zero label, select by area, otherwise its
%obvious which label we need to keep
if length(labels_2_check)>1 
    area_sizes=zeros(size(labels_2_check));
    for i=1:length(labels_2_check)
        tmp=I0mask_cc==labels_2_check(i);
        area_sizes(i)=sum(tmp(:));
    end
    %biggest object to keep is...
    label_2_keep=labels_2_check(area_sizes==max(area_sizes));
    %and if there just happens to be 2 regions of EXACTLY the same size, just
    %keep the first one, since its on the left an the images are oriented
    %that the breast will be on the left
    label_2_keep=label_2_keep(1);
else
    label_2_keep=labels_2_check;
end

I0mask=double(I0mask_cc==label_2_keep);

for y = 1:dimr
    i=find(I0mask(y,:),1,'last');
    if isempty(i)
        bdrind(1,y) =  1;
    else
        bdrind(1,y) = i;
    end
end

%Cut all points below "narrowest"
bb=min(bdrind(:));
cutpoint=find(bdrind(1,floor(dimr/2):end)==bb,1)+floor(dimr/2)-1;
bdrind(1,cutpoint:end)=1;
