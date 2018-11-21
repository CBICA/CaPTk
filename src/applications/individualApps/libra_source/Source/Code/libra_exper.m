function [results, breastMask] = libra_exper(dcm_fname,outdir,outtxt,saveIntermediate)
%
%  This function checks dicom header information and imaging physics,
%  segments, by calling libra_run, breast and dense tissue for a given
%  input breast dicom image.
%
%  Input arguments:
%       dcm_fname          <char>  Absolute path to the dicom images.
%       outdir             <char>  Absolute path to the output folder. If
%                                  not provided, outdir will be
%                                  "Result_Images" under current working
%                                  directory.
%       outtxt             <char>  Absolute path to the output CSV file. If
%                                  not provided, outtxt will be
%                                  "Density.csv" under current working
%                                  directory.
%       saveIntermediate   <bool>  Boolean value to save the intermediate
%                                  in the output (1) or not (0) (Optional.
%                                  Default: 0, not saving.)
%
%  Output arguments:
%       results            breast segmentation values: BreastArea(sqcm),
%                          DenseArea(sqcm), and BreastDensity(%).
%       breastMask         Binarized breast segmentation mask.
%
%  Usage:
%   >> [results, breastMask] =
%   libra_exper(dcm_fname,outdir,outtxt,saveIntermediate) to process a
%   dicom image and save the outputs in outdir and in outtxt. Depending on
%   value in saveIntermediate, the intermediate files could be kept in
%   outDir.
%
%  Version info:
%  $Rev: 450 $:  Revision of last commit
%  $LastChangedBy: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:  Date of last commit
%
%  Developed by Dr. B.M. Keller
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

warning off; %#ok<WNOFF>

%% input checking
if nargin<4
    saveIntermediate=0;
end

if (~exist('outdir','var'))
    outdir=fullfile(pwd,'Result_Images');
end

%Open results file:
if (~exist('outtxt','var'))
    outtxt=fullfile(pwd,'Density.csv');
end
% end of input checking

%% initialization
if (~exist(outtxt,'file'))
    need_header_row_flg=1;
else
    need_header_row_flg=0;
end

fid=fopen(outtxt,'a');
if need_header_row_flg
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n','File Analyzed','Manufacturer','Laterality','ViewPosition','BreastArea(sqcm)','DenseArea(sqcm)','BreastDensity(%)');
end

%disp(strcat(['Analyzing Image: ' dcm_fname]));

%Rootname for output images
[~,outfileroot]=fileparts(dcm_fname);
[outtxtdir,~,~]=fileparts(outtxt);
outroot=fullfile(outdir,strcat([outfileroot '.jpg']));
mask_fileout=fullfile(outdir,strcat(['Masks_' outfileroot]));
skiptxt=fullfile(outtxtdir,'SkippedImg.csv');
% end of initialization

%Now start processing image - read the image and DICOM header into memory
dcmimg=double(dicomread(dcm_fname));
dcminfo=dicominfo(dcm_fname);

%% Series of checks on dicom header/image
checkMammo(outfileroot, dcminfo, skiptxt);
dcminfo=checkPhysics(outfileroot, dcminfo, skiptxt);
% end of Series of checks on dicom header

pecseg_flag=strcmpi(dcminfo.ViewPosition(1),'M');

%% Standardizing intensity
% %Raw needs to be log transformed, inverted and squared for preprocessing,
% in that order
% %Monochrome1 just needs inversion - middle step
% step 1 - raw gets log transform
if strcmp(dcminfo.PresentationIntentType,'FOR PROCESSING')
    disp('...Pre-processing raw mammogram')
    if min(min(dcmimg(:)))<1
        dcmimg=dcmimg+abs(min(min(dcmimg(:))))+1;
    end
    dcmimg=log(dcmimg);
end

%Step 2 - Do pixel inversion if Needed
if strcmp(dcminfo.PhotometricInterpretation,'MONOCHROME1')
    %Need to invert intensities
    %ensure min pixel value is at least "1"
    disp(strcat(['...inverting pixel intensities: ' outfileroot]));
    dcmimg=abs(dcmimg-max(dcmimg(:)));
end

% Step 3 - square transform for density contrast
if strcmp(dcminfo.PresentationIntentType,'FOR PROCESSING')
    dcmimg=dcmimg.^2;
end

dcmimg_orig=dcmimg;
dcmimg=imresize(dcmimg,0.25,'bicubic'); %0.25 or 0.5 ; 'nearest' or 'bicubic'

%If the ImageLaterality field isn't there for some reason, check the
%Laterality field and copy over (because we check .ImageLaterality
if ~isfield(dcminfo,'ImageLaterality') && isfield(dcminfo,'Laterality')
    dcminfo.ImageLaterality=dcminfo.Laterality;
end

% check seriestime field
if ~isfield(dcminfo,'SeriesTime');
    dcminfo.SeriesTime='NA';
else
	% to remove comma as it will go into a csv file.
    dcminfo.SeriesTime=strrep(dcminfo.SeriesTime,',',' ');
end

%% Standardizing orientation
% check whether the image is horizonally flipped
isflipped=isfield(dcminfo,'FieldOfViewHorizontalFlip') && strcmp(dcminfo.FieldOfViewHorizontalFlip,'YES');
if isflipped
    if strcmp(dcminfo.ImageLaterality,'L')
        dcmimg_orig=fliplr(dcmimg_orig);
    end
else
    if strcmp(dcminfo.ImageLaterality,'R')
        dcmimg_orig=fliplr(dcmimg_orig);
    end
end
%Now which side
if strcmp(dcminfo.ImageLaterality,'R')
    side='Right';
elseif strcmp(dcminfo.ImageLaterality,'L')
    side='Left';
else
    side='Special-View';
end
% End of Standardizing orientation

%% Checking image type and vendor
%Lastly - Check Presentation Intent Type, and vendor - for selecting
%training data
type_flg=1;
disp('...identifying appropriate training model');
%first determine  the presentation intent type
if isfield(dcminfo,'PresentationIntentType')
    if strcmp(dcminfo.PresentationIntentType,'FOR PRESENTATION')%Processed
        intent_type='Processed';
    elseif strcmp(dcminfo.PresentationIntentType,'FOR PROCESSING') %RAW
        intent_type='Raw';
    else %unknown field value
        intent_type='Processed'; %most commonly available
        type_flg=0;
    end
else % no field in the header. so choose default
    intent_type='Processed'; %most commonly available
    type_flg=0;
end

%next ID the mammogram Manufacturer
if isfield(dcminfo,'Manufacturer')
    if ~isempty(strfind(dcminfo.Manufacturer,'HOLOGIC'))
        manufact_name='Hologic';
    elseif ~isempty(strfind(dcminfo.Manufacturer,'GE MEDICAL')) || ~isempty(strfind(dcminfo.Manufacturer,'FUJIFILM'))
        manufact_name='GE';
    else %unknown field value
        manufact_name='Hologic'; %default is Hologic
        type_flg=0;
    end
   	% to remove comma as it will go into a csv file.
    dcminfo.Manufacturer=strrep(dcminfo.Manufacturer,',',' ');
else % no field in the header. so choose default
    dcminfo.Manufacturer='NA';
    manufact_name='Hologic'; %default is Hologic
    type_flg=0;
end
% End of Checking image type and vendor

%% Loading SVM model
%now set the training data to be loaded
training_data=strcat(manufact_name,'_',intent_type,'_SVM.mat'); %The trained SVM Model
if type_flg %display chosen model
    disp(strcat('...using model: ',manufact_name,'-',intent_type,'...'))
else %tell user that we had to choose a default because intent and//or manufacturer info is missing from the dicom header
    disp('...*unknown/unavailable presentation intent type and/or manufacturer training available...');
    disp(strcat('...using default model: ',manufact_name,'-',intent_type,'...'));
end

try % try-catch not necessary. Remove
    load(training_data);
catch err
    disp(err.identifier);
    disp(err.message);
    disp(['   > Failed on loading SVM model: ',training_data]);
    error('   > Exiting.')
end
% End of loading SVM model

tic;
%NOW SEGMENT DENSITY!
[results, breastMask]=libra_run(double(dcmimg),outroot,mask_fileout,dcminfo,SVMmodel_vars+1,svm_struct,dcmimg_orig,bw,gw,sig,pecseg_flag,k,saveIntermediate);

fprintf(fid,'%s,%s,%s,%s,%f,%f,%f\n',outfileroot, dcminfo.Manufacturer, side, dcminfo.ViewPosition, results(1:3));
disp(strcat(['Est PD%: ' num2str(results(3)) '%']));
toc;

fclose(fid);

warning on; %#ok<WNON>
end % end of libra_exper.m

function checkMammo(id,dcminfo,txt)
	disp('   > Checking header');
	if isfield(dcminfo,'ImagesInAcquisition')
		TwoDim_check=dcminfo.ImagesInAcquisition==1;
	else
		TwoDim_check=1; % process it anyway
	end

	if isfield(dcminfo,'Modality')
		MG_check=strcmp(dcminfo.Modality,'MG');
	else
		MG_check=1;
	end

	checks=[ TwoDim_check , MG_check ];
	if ~(TwoDim_check && MG_check)
		printoutCSV(id, checks, txt)
		disp('   > This file is not 2D Mammogram, or has multiple slices in acquisition. Skipping');
		error('   > Printing out image info into SkippedImg.csv.');
	end
end % end of checkMammo


function dcminfo_out = checkPhysics(id, dcminfo, txt)
	% New "dicom header fixes" preamble
	% (i.e., because some places blank out age), but there needs to be a valid
	% value. So if something is missing, lets try imputting a value (i.e., a mean value)
	% Michael: moved this part to libra_exper where header checks are, also
	% made it robust to detecting variable containing only spaces.
	if ~isfield(dcminfo,'PatientAge') || isempty(strtrim(dcminfo.PatientAge)) || length(dcminfo.PatientAge)<3
		disp('     > WARNING: Age information is missing. Use Default.');
		dcminfo.PatientAge='055Y'; %default "average value" for age
	end

requiredFields={'PixelSpacing','BodyPartThickness', ...
		'CompressionForce','Exposure','KVP'};

if ~isfield(dcminfo,'PixelSpacing')
    if isfield(dcminfo,'ImagerPixelSpacing') 
        %Use PixelSpacing (calibrated) if available, default to
        %ImagerPixelSpacing (detector resolution) if not. One needs to be
        %there, eitherway
        disp('     > WARNING: PixelSpacing information is missing. Defaulting to detector resolution');
        dcminfo.PixelSpacing=dcminfo.ImagerPixelSpacing;
    else
        disp('     > WARNING: Image Resolution information is missing');
    end
end
        

	if ~isfield(dcminfo,'ViewPosition')
		if ~isfield(dcminfo.ViewCodeSequence.Item_1,'CodeMeaning')
			requiredFields=[requiredFields,'ViewPosition'];
		elseif isfield(dcminfo.ViewCodeSequence.Item_1,'CodeMeaning')
			dcminfo.ViewPosition=dcminfo.ViewCodeSequence.Item_1.CodeMeaning;
		end
	end
	% to remove comma as it will go into a csv file.
	dcminfo.ViewPosition=strrep(dcminfo.ViewPosition,',',' '); 
	
	if ~isfield(dcminfo,'XrayTubeCurrent')
		if ~(isfield(dcminfo,'ExposureInuAs') && isfield(dcminfo,'ExposureTime'))
			requiredFields=[requiredFields,'XrayTubeCurrent'];
		elseif (isfield(dcminfo,'ExposureInuAs') && isfield(dcminfo,'ExposureTime'))
			disp('     > WARNING: XrayTubeCurrent information is missing. Calculate from ExposureInuAs/ExposureTime.');
			dcminfo.XrayTubeCurrent=dcminfo.ExposureInuAs/dcminfo.ExposureTime;
		end
	end
	numFields=size(requiredFields,2);
	isField=zeros(1,numFields);

	for i=1:numFields;
		isField(i)=any(strcmp(requiredFields{i},fieldnames(dcminfo)));
	end

	if(all(isField))
		disp('   > Checking passed.');
		dcminfo_out=dcminfo;
	else
		printoutCSV(id, isField, txt, requiredFields)
		disp('   > Imaging physics information not complete.');
		error('   > Printing out image info into SkippedImg.csv.');
	end
end % end of checkPhysics

function printoutCSV(id, isField, txt, requiredFields)
	if (~exist(txt,'file'))
		s_need_header_row_flg=1;  % To distinguish with the variable in the global scope, s_* prefix means it's defined and only visible in this scope.
	else
		s_need_header_row_flg=0;
	end
	s_fid=fopen(txt,'a');
	if s_need_header_row_flg
		% start from here, fix the header row in a nice way
		fprintf(s_fid,'%s,%s,%s\n','File Skipped','isMammo','MissingFields');
	end

	if length(isField)<3;
		isMammo='No';
		MissingFields='';
	elseif length(isField)==length(requiredFields);
		isMammo='Yes';
		MissingFields(isField==0)=requiredFields(isField==0);
		MissingFields(cellfun('isempty',MissingFields))=[];
	end

	str_missing='';
	for s=1:size(MissingFields,2)
		str_missing=[str_missing MissingFields{s} '; '];
	end

	% Print it out to SkippedImg.csv
	fprintf(s_fid,'%s,%s,%s\n',id,isMammo,str_missing); % not the best way to display the missing fields.
	fclose(s_fid);
end % end of printoutCSV
