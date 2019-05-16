function [results, breastMask]=libra_run(dcmimg,outroot,mask_fileout,dcminfo,SVMmodel_vars,svm_struct,dcmimg_orig,bw,gw,sig,pecsegflag,max_k,saveIntermediate)
%
%  !!! One should NOT call this function directly !!!
%
%  This function segments breast, classify each group of tissue into
%  dense/non-dense tissue for a given input breast dicom image.
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


%% Segment the breast
if ~exist('outroot','var')
    outroot=fullfile(pwd,'Result_Images','output.jpg');
end
[outpath,outfileroot,ext]=fileparts(outroot);
if ~exist(outpath,'dir')
    mkdir(outpath);
end

vis='off';
results=zeros(1,3);

disp('Masking out Breast')

%% Standardizing orientation
% Michael's modification for image flipping. All image will be flipped to
% the left in the view to continue the process.
flipneed=0;
isflipped=isfield(dcminfo,'FieldOfViewHorizontalFlip') && strcmp(dcminfo.FieldOfViewHorizontalFlip,'YES');
if isflipped
    if strcmp(dcminfo.ImageLaterality,'L')
        flipneed=1;
    end
else
    if strcmp(dcminfo.ImageLaterality,'R')
        flipneed=1;
    end
end % End of Michael's modification

if flipneed==1
    %if we need to flip it
    I0=fliplr(double(dcmimg));
    dcmimg=fliplr(double(dcmimg));
else
    I0=double(dcmimg);
end
% End of Standardizing orientation

%% Breast segmentation
getmask_flg=1;
tmp=I0;
multiplier=0;
while getmask_flg && multiplier<1
    [bdrind, airthresh] = con_breast(I0,tmp);
    
    %%% pectoral segmentation only really needed in MLO view
    if pecsegflag
        % Choose ROI for HT
        [roih] = croi(I0, bdrind,1);
        roih=imdilate(roih,strel('disk',5,0));
        roih=imfilter(roih,fspecial('gaussian',[5 5],1));
        %ignore bottom right quadrent (cut corfcner to corner) other=NaN
        roih(round(size(roih,1).*0.50):end,round(size(roih,2).*0.50):end)=NaN;
        roih(roih<=airthresh)=NaN; %ignore air
        
        % Find pectoral muscle using Hough Transformation
        [peind] = pectoral(I0,roih);
        
    else
        peind=zeros(size(I0,1),1);
    end
    
    % Store breast region information for fuzzy c mean method
    [mask] = breast(I0,peind,bdrind,airthresh);
    
    if sum(mask(floor(size(I0,1)*0.5),:))==0
        %sanity check
        %If the breast mask doesn't include cover a good chunk of the 
        %y-extent of the image (i.e. the breast generally is in the middle of the image), 
        %somehing went wrong and we need to redo; most always this is due to the air
        %thresholding failing, likely because of an extra dense pec region, the
        %solution is to lookk just at the right part of the image to do the
        %thresholding, keep shifting right till it works
        multiplier=multiplier+0.3;
        tmp=I0(:,floor(size(I0,2)*multiplier):end);
    else
        %it worked, w're good to go
        getmask_flg=0;
    end
end

%Need to break out here if the breast mask generation fails completely, to add

%Now add in air-breast boundary correction for GE processed
if isfield(dcminfo,'Manufacturer') && isfield(dcminfo,'PresentationIntentType')
    if strcmp(dcminfo.Manufacturer,'GE MEDICAL SYSTEMS') &&  strcmp(dcminfo.PresentationIntentType,'FOR PRESENTATION')
       peind=floor(peind);
       pecmask=zeros(size(mask));
       for i=1:length(peind) 
           pecmask(i,1:peind(i))=1;
       end
       M2=imerode(mask+pecmask,strel('disk',5,0));
       mask=M2.*mask;
    end
end

%Last step, incase the breast mask is in discontinous parts (e.g, the
%breast and a little of the chest in MLO or the other breast in CC is also
%masked), just keep the big mask
L=bwlabeln(mask);
mask=double(L==mode(L(L>0)));

%erode mask avoid edge issues messing up the histogram
%- but we want the non-eroded area to be our "breast area" for determining
%breast percent-density (and not the area of the downsampled mask, but the
%fullsize mask)
maskarea=sum(sum(imresize(mask,size(dcmimg),'nearest')));
mask2=imerode(mask,strel('disk',15));

mask(mask==0)=NaN;
mask2(mask2==0)=NaN;
% End of Breast segmentation

%% Unsupervised clustering of the mammogram
% I0 and dcmimg are normalized and downsampled image
masked_img=I0.*mask2; %for purposes of fcm-cluster, ignore values of outer 15 pixels

%get values WITHIN mask
tmp=masked_img(:);
non_nan_rows=~isnan(tmp);
%get rid of nan's everything else=foreground aka breast
tvalues=tmp(non_nan_rows);

if size(tvalues,1)==1
    %if for whatever reason, datapoint coords are by column and not
    %row, flip
    tvalues=tvalues';
end
tvalues=sortrows(tvalues,1);

%normalize
stdev=std(tvalues);
mue=mean(tvalues);
tvalues=(tvalues-mue)./stdev;

%...grab range for later windowing
lo_win=tvalues(round(length(tvalues)*0.05));
hi_win=tvalues(round(length(tvalues)*0.95));

%cut at z-score -4 to 4, outside that range=outliers
tvalues=tvalues(find(tvalues>=max(-4,min(tvalues)),1,'first'):find(tvalues<=min(4,max(tvalues)),1,'last'));

%sub-sample the histogram based on its CDF to speed up the FCM
tvalues = subsample(tvalues)';

%clear unneeded variables now
clear t non_nan_rows z I0 masked_img masked_img2 masked_img3 tmp mask2 roih

%Cluster
%find appropriate "k": number of clusters
[k]=zero_crossings(tvalues,bw,gw,sig);
if k<2
    k=2;
elseif k>max_k
    k=max_k;
end
%pre-allocate biggest array now
dense_member=single(zeros([size(dcmimg) k]));

disp(strcat(['Running ' num2str(k) '-class FCM']));
fzy=2;
[cents]=fcm1d(tvalues,k,fzy, 100, 0.0001); %get cluster centers via FCM

%Now determine membership
cents=sortrows(cents,1); %1=lowest (fat) first, -1=highest (dense first)

%determine pixel wise membership - on fullsize img&mask
%apply zscore on full-size input image and mask it
mask=imresize(mask,size(dcmimg),'nearest'); % mask and dcmimg are supposed to be the same size, no?
corrected_img=(double(dcmimg)-mue)./stdev;
masked_img=corrected_img.*mask;

%get distance image
dense_class=uint8(zeros(size(dcmimg)));
%Make distance_2_cluster_center_images
for i=1:length(cents)
    dense_member(:,:,i)=(cents(i)-masked_img).^(-2/(fzy-1));
end
clear masked_img

%Now make membership images
tdist=sum(dense_member,3);
for i=1:length(cents)
    dense_member(:,:,i)=(dense_member(:,:,i)./tdist).*mask;
    %figure(i);imagesc(dense_member(:,:,i));
end
%make everything outside of mask -0.1. Michael: why?
dense_member(isnan(dense_member))=-0.1;

%now get results
%Whos the max
for i=1:length(cents)
    dense_class(:,:)=dense_class+uint8((dense_member(:,:,i)==max(dense_member,[],3)).*i.*mask);
    %figure(i);imagesc(dense_class);
end
dense_class(isnan(dense_class))=0;
clear dense_member tdist
%which cents are closer to dense cluster foo(k) than fat cluster foo(1)

%Sanity check, if upper cluster (which is index k, because it is sorted) is
%<1 cm2 (in original resolution), merge it with next one and check again, 
%then update everything,
flg_uppercluster_size=1;
while flg_uppercluster_size && k>2 %to prevent infinite while loop
    tmp=imresize(logical(dense_class==k),size(dcmimg_orig),'nearest');
    uppercluster_size=sum(tmp(:)).*dcminfo.PixelSpacing(1).*dcminfo.PixelSpacing(2).*0.1.*0.1;
    if uppercluster_size<1 %#ok<BDSCI> %there should almost always be at least 1 cm2 of dense tissue
        %have to supress down so high Intensity pixels of cluster k in
        %corrected_img (i.e., the z-scored img), so as to not screw up some
        %feature calculations down-stream
        Imax_nextcluster=max(corrected_img(dense_class==k-1));
        corrected_img(dense_class==k)=Imax_nextcluster;
        %now to account for going from k to k-1 clusters
        dense_class(dense_class==k)=k-1;%merge cluster k into cluster k-1
        cents(k)=[];%make centroid of cluster k go away
        k=k-1; %merging means 1 less total cluster
    else
        flg_uppercluster_size=0; %sane cluster size, so breakout        
    end
end
clear tmp flg_uppercluster_size flg_uppercluster_size

%tell me all the pds
pdk=zeros(1,k); %the pd% possibilities by cluster
for i=1:k
    seg=logical(dense_class>=i);
    densearea=sum(sum(seg));
    pdk(i)=densearea/maskarea;
end

%Sanity check, if upper cluster is <0.5 cm2, merge it with next one and
%check again, then adjust all 
%% Now lets get the features
% tell me x and y extents of clusters (how wide/tall they are)
xext=zeros(1,k);
yext=zeros(1,k);
for i=1:k
    seg=logical(dense_class>=i);
    xext(i)=max(sum(seg,2))/size(dcmimg,2); % Again, dcmimg is a downsampled image
    yext(i)=max(sum(seg,1))/size(dcmimg,1);
end

max_diff_yext=find(abs(diff(yext))==max(abs(diff(yext))),1,'last')+1;
max_diff_xext=find(abs(diff(xext))==max(abs(diff(xext))),1,'last')+1;

%add GLCM texture features to features_by_aggl
glcm_features=zeros(1,4);
m=dense_class>=1;m_img=m.*dcmimg;%has to be integer
%crop m_img to tight roi: chest to nipple
qqq=max(m,[],1);
clow=find(qqq,1,'first');chigh=find(qqq,1,'last');
qqq=max(m,[],2);
rlow=find(qqq,1,'first');rhigh=find(qqq,1,'last');
m_img=imcrop(m_img,[clow rlow chigh-clow rhigh-rlow]);

glcm=graycomatrix(m_img,'NumLevels',16,'GrayLimits',[lo_win hi_win],'Offset',[0 1; -1 1; -1 0; -1 -1]);
stats=graycoprops(glcm);
glcm_features(1,:)=[mean(stats.Contrast), mean(stats.Correlation),mean(stats.Energy),mean(stats.Homogeneity)];
entropy_value=entropy(m_img);

mask(isnan(mask))=0;
%Shape descriptors of the breast mask
shape_stats=zeros(1,11);
bwstats=regionprops(double(mask),'all');
shape_stats(1)=bwstats.Area;
shape_stats(2)=bwstats.MajorAxisLength;
shape_stats(3)=bwstats.MinorAxisLength;
shape_stats(4)=bwstats.Eccentricity;
shape_stats(5)=bwstats.EquivDiameter;
shape_stats(6)=bwstats.EulerNumber;
shape_stats(7)=bwstats.Orientation;
shape_stats(8)=bwstats.ConvexArea;
shape_stats(9)=bwstats.Solidity;
shape_stats(10)=bwstats.Extent;
shape_stats(11)=bwstats.Perimeter;

histstats=[moment(tvalues,3) moment(tvalues,4) moment(tvalues,5) range(tvalues) skewness(tvalues) kurtosis(tvalues) k str2num(dcminfo.PatientAge(1:3)) dcminfo.BodyPartThickness dcminfo.CompressionForce dcminfo.Exposure dcminfo.XrayTubeCurrent dcminfo.KVP mue stdev lo_win hi_win max_diff_xext max_diff_yext glcm_features entropy_value shape_stats]; %#ok<ST2NM>
features_by_aggl=repmat([0 0 0 histstats],k,1);ccc=3+length(histstats)+1;
% ccc is to keep a record on which column we're at.
%ccc=4 because column 1 will be filled in with case number, 
% column 2 with img_count (1 and 2 are filled by program that calls this) and column 3 with 1:k
features_by_aggl(:,3)=1:k;

%OK Now append features_by_aggl with hierarchical agglomeration features
%ie features on densest cluster, densest 2 clusters (i.e., potential
%segmentations)
disp('...Getting dense-cluster metrics (i.e. hierarchical agglomeration)');
f=zeros(1,length(cents));  % foo is the cluster centers
b=dense_class;a=corrected_img;%a=dcmimg;
%1st order hist features
%...mean
for i=1:length(cents)
    m=b>=i;
    z=a.*double(m);
    z(z==0)=NaN;z=z(1:end); % flatten z
    f(i)=nanmean(z);
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%...variance
for i=1:length(cents)
    m=b>=i;
    z=a.*double(m);
    z(z==0)=NaN;z=z(1:end);
    f(i)=nanvar(z);
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%...skewness
for i=1:length(cents)
    m=b>=i;
    z=a.*double(m);
    z(z==0)=NaN;z=z(1:end);
    f(i)=skewness(z);
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%...Kurtosis
for i=1:length(cents)
    m=b>=i;
    z=a.*double(m);
    z(z==0)=NaN;z=z(1:end);
    f(i)=kurtosis(z);
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%...range
for i=1:length(cents)
    m=b>=i;
    z=a.*double(m);
    z(z==0)=NaN;z=z(1:end);
    f(i)=range(z);
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%...Intensity mean diff
for i=2:length(cents)
    m=b>=i;m2=b<i;
    z=a.*double(m);z2=a.*double(m2);
    z(z==0)=NaN;z=z(1:end);
    z2(z2==0)=NaN;z2=z2(1:end);
    f(i-1)=nanmean(z)-nanmean(z2);
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%Other stuff
%compactness
for i=1:length(cents)
    m=b>=i;
    f(i)=sum(sum(bwperim(m)))/sum(sum(m));
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%Y-extent
for i=1:length(cents)
    m=b>=i;
    f(i)=sum(max(m,[],2))./sum(max(b>=1,[],2));
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%x-extent
for i=1:length(cents)
    m=b>=i;
    f(i)=sum(max(m,[],1))./sum(max(b>=1,[],1));
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%Connected-ness (# of binary regions in agglomerated mass)
for i=1:length(cents)
    m=b>=i;
    f(i)=length(unique(bwlabel(m).*m))-1;
end
features_by_aggl(:,ccc)=f';ccc=ccc+1;

%Shape descriptors
shape_stats=zeros(length(cents),11);
for i=1:length(cents)
    m=b>=i;
    bwstats=regionprops(double(m),'all');
    shape_stats(i,1)=bwstats.Area;
    shape_stats(i,2)=bwstats.MajorAxisLength;
    shape_stats(i,3)=bwstats.MinorAxisLength;
    shape_stats(i,4)=bwstats.Eccentricity;
    shape_stats(i,5)=bwstats.EquivDiameter;
    shape_stats(i,6)=bwstats.EulerNumber;
    shape_stats(i,7)=bwstats.Orientation;
    shape_stats(i,8)=bwstats.ConvexArea;
    shape_stats(i,9)=bwstats.Solidity;
    shape_stats(i,10)=bwstats.Extent;
    shape_stats(i,11)=bwstats.Perimeter;
end
features_by_aggl(:,ccc:ccc+10)=shape_stats;ccc=ccc+11;

%add GLCM and entropy to features_by_aggl
glcm_features=zeros(k,5);
for i=1:length(cents)
    m=b>=i;m_img=m.*dcmimg;%has to be integer
    %crop m_img to tight roi
    qqq=max(m,[],1);
    clow=find(qqq,1,'first');chigh=find(qqq,1,'last');
    qqq=max(m,[],2);
    rlow=find(qqq,1,'first');rhigh=find(qqq,1,'last');
    m_img=imcrop(m_img,[clow rlow chigh-clow rhigh-rlow]);
    
    glcm=graycomatrix(m_img,'NumLevels',16,'GrayLimits',[lo_win hi_win],'Offset',[0 1; -1 1; -1 0; -1 -1]);
    stats=graycoprops(glcm);
    glcm_features(i,:)=[mean(stats.Contrast), mean(stats.Correlation),mean(stats.Energy),mean(stats.Homogeneity),entropy(m_img)];
end
features_by_aggl(:,ccc:ccc+size(glcm_features,2)-1)=glcm_features;

%% Classifying Fatty and Dense, and compute results
fullset=features_by_aggl;

%compute class differences
df=fullset(2:end,4:size(fullset,2))-fullset(1:end-1,4:size(fullset,2));
df=cat(1,zeros(1,size(fullset,2)-3),df);
features_by_aggl=[fullset(:,1:end) df];

%features_by_aggl(:,2)=pdhat;
svm_guess=svmclassify(svm_struct,features_by_aggl(:,SVMmodel_vars));
svm_guess(1)=0;svm_guess(end)=1;%Because because no one is 100% dense, but there is always some density (1 is dense)
answer=find(svm_guess==0, 1,'last')+1;%last contiguous set of dense clusters - this implementation garuantees answer>=2


seg=imresize(logical(dense_class>=answer),size(dcmimg_orig),'nearest');
densearea=sum(sum(seg)).*dcminfo.PixelSpacing(1).*dcminfo.PixelSpacing(2).*0.1.*0.1;

pd=pdk(answer);

mask_up=imresize(mask,size(dcmimg_orig),'nearest');
breastarea=sum(sum(mask_up)).*dcminfo.PixelSpacing(1).*dcminfo.PixelSpacing(2).*0.1.*0.1;

results(1)=breastarea;
results(2)=densearea;
results(3)=pd*100;

%% Save Outpute Images
corrected_img=(dcmimg_orig-mue)./stdev;
%Reverse the segmentations back so they align with original image
if flipneed
    seg=fliplr(seg);
    mask_up=fliplr(mask_up);
    dense_class=fliplr(dense_class);
    corrected_img=fliplr(corrected_img); % Michael: When the image was flipped by the program previously, it has to be flipped again now.
end

mask_toWrite = uint8(imresize(mask_up,size(dcmimg_orig),'nearest'));
disp('Trying to write the result');
maskDir = strcat(outpath,'/totalmask');
mkdir(maskDir);
dicomwrite(mask_toWrite, fullfile(maskDir,'/totalmask.dcm'), dcminfo, 'CreateMode', 'Copy');

% %%%%save only outlined images
figure(1);set(1,'Visible',vis);imagesc(dense_class);colormap jet;colorbar;axis image;
V1=min(tvalues);V2=max(tvalues);N=-3;
[a,c]=hist(tvalues, round(V1*10^(-1*N))*10^N:bw: round(V2*10^(-1*N))*10^N);
w=gw;g=gausswin(w+1,sig);%smooth hist
b=conv(a(1:end),g,'same');
figure(2);set(2,'Visible',vis);bar(c,b./length(tvalues));
yLims = get(gca,'YLim');
hold on
for i=1:length(cents)
    line([cents(i),cents(i)],yLims,'Color','red')
end
hold off

%[imgwind, outlined_image]=wind_and_outline(corrected_img,seg,mask_up,lo_win,2.5);
%imwrite(outlined_image,fullfile(outpath,strcat([outfileroot '_density_segmentation' ext])),ext(2:end));
breastMask=logical(mask_up);

if saveIntermediate
    saveas(1,fullfile(outpath,strcat([outfileroot '_density_imagesc' ext])));
    saveas(2,fullfile(outpath,strcat([outfileroot '_intensity_histogram' ext])));
    imwrite(imgwind,fullfile(outpath,strcat([outfileroot '_Windowed_Original' ext])),ext(2:end));
end    

% saving mat file
res=struct();
res.PD_SVM=pd*100;
res.BreastArea=breastarea;
res.DenseArea=densearea;
res.BreastMask=breastMask;
res.DenseMask=seg;
res.dcm_fname=outfileroot;
save(mask_fileout,'res');

close all;
