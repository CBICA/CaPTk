%  Version info:
%  $Rev: 455 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 16:12:51 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

close all;clc;
filePath = mfilename('fullpath');
[sourceDir,~,~] = fileparts(filePath);
[parentDir,~,~] = fileparts(sourceDir);

%Set output directories, easiest in current folder, not set-up to handle
%directories outside of current folder
testdir=fullfile(sourceDir,'Demo_Test');
outdir=fullfile(testdir,'Result_Images');
mkdir(outdir)
outtxt=fullfile(testdir,'Density.csv'); %Comma-separated file to store Numeric Results
saveIntermediate=0;

% load a ground-truth mat file that has the numbers for each test images.
groundTruth=fullfile(parentDir,'Sample_Data','ground_truth.mat');
load(groundTruth);
failure=0;

idir=fullfile(parentDir,'Sample_Data');

% Running tests on individual images.
for i=1:5
	dcm_fname=fullfile(idir,['Case',num2str(i),'.dcm']);%Full path to data
	[~,NAME,EXT]=fileparts(dcm_fname);
	disp(strcat(['>>>> Performing Density Analysis on Image: ' NAME EXT]));
	
	%Run Density Analysis Code
	try
		[T, result] = evalc('libra(dcm_fname,testdir,saveIntermediate)');
	catch err
        disp(err.identifier);
        disp(err.message);
        disp(['  > Testing of libra on ' NAME ' failed.']);
	end
	
	regression=abs((cell2mat(result(2,2:4))-truth(i,:))./truth(i,:));
	if all(regression < 0.05)  % allows 5% error
        disp(['  > Regression test on Image: ' NAME EXT ' Passed']);
	else
        disp(['  > Regression test on Image: ' NAME EXT ' Failed']);
        failure=failure+1;
	end
end

fprintf('\n');
fprintf('\n');
disp('>>>> Regression test failure report <<<<');
disp(['  > ' num2str(failure) ' has failed during image regression tests.']);
fprintf('\n');
fprintf('\n');

%% Running a test on batch processing
%Set output directories, easiest in current folder, not set-up to handle
%directories outside of current folder
testdir=fullfile(sourceDir,'Demo_Test_Batch');
outdir=fullfile(testdir,'Result_Images');
mkdir(outdir)
outtxt=fullfile(testdir,'Density.csv');
disp(strcat(['>>>> Performing Batch Density Analysis on Sample_Data: ' idir]));
try
	[T,result] = evalc('libra(idir,testdir,saveIntermediate)');
	disp(['  > Batch Density Analysis Finished']);
catch err
    disp(err.identifier);
    disp(err.message);
    disp('  > Testing of libra failed.');
end


