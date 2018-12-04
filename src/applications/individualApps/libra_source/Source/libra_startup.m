%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

filePath = mfilename('fullpath');
[sourceDir,~,~] = fileparts(filePath);

if (~exist(sourceDir,'dir'))
	error('  Error determining the paths. Exiting...');
end

warning off
rmpath(genpath(sourceDir));
warning on

disp('Adding paths...');

addpath(sourceDir);
addpath(fullfile(sourceDir,'Code'));
addpath(fullfile(sourceDir,'Model'));
%addpath(fullfile(sourceDir,'Sample_Data'));

matlabVersion=version('-release');

libraVersion=libra_version;
disp('Welcome to LIBRA!');
disp(['You are using LIBRA-' libraVersion ]);

