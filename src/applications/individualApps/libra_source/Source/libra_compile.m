function libra_compile(install_prefix)
%  Compile LIBRA into executable to a given path
%
%  Argument:
%       install_prefix  <char>  Absolute or relative path to the destination directory
%                               in which "libra" executable will be generated.
%                               Default: an "Executable" directory will be created under
%                               current working directory as the output directory.
%
%  Version info:
%  $Rev: 450 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2015-11-20 11:34:20 -0500 (Fri, 20 Nov 2015) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

libra_startup;
	
filePath = mfilename('fullpath');
[sourceDir,~,~] = fileparts(filePath);

if ~exist('install_prefix','var')
	install_prefix=fullfile(sourceDir,'Executable');
end

modelDir=fullfile(sourceDir,'Model');
mkdir(install_prefix);

disp('Compiling source codes into a binary... It might take a couple of minutes.');
try 
	eval(['mcc -m -R -singleCompThread -R -nosplash -a ''', modelDir, ''' -d ''', install_prefix, ''' -o libra libra.m']);
catch err
	disp(err.identifier);
        disp(err.message);
        error('Compilation failed.');
end

outputBinary=fullfile(install_prefix,'libra');
if isunix
	% Change file permission (chmod 755)
	fileattrib(outputBinary,'+x','a');
	fileattrib(outputBinary,'+w','u');
end
disp(['Compilation finished. The binary is ' outputBinary '.']);
