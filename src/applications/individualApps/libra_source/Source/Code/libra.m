function [output] = libra(input,outDir,saveIntermediate)
%
%  This is a wrapper scipt for processing using libra_exper.m, used for
%  compiling an .exe for distribution. This software is designed to either
%  run interactively (if only libra.exe is called) or non-interactively if
%  arguments corresponding to input and output files/directories are
%  specified in the command line interface.
%
%  The program will perform breast density estimation on the input dicom
%  file or dicom files within a given directory and output three values,
%  breast area(sqcm), dense area(sqcm), and breast density(%). These values
%  are stored in a Density.csv under outDir. Final breast segmentation is
%  saved as a jpg image in subdirectory Result_Images/ in outDir.
%
%  Input arguments -optional: if no arguments are passed, the software runs
%  via a series of prompts for the user to specify the input (a single
%  DICOM File, or directory with mutliple for batch processing) and the
%  target output directory. If arguments are passed, the software will run
%  non-interactively and process the single file specified (or all files
%  under the directory provided) by <input>, saving the results to the
%  specified output folder <outDir>.
%       input              <char>  Absolute path to the dicom file or
%                                  directory containing dicom images
%                                  (Optional)
%       outDir             <char>  Absolute path to the output directory
%                                  (Required if <input> is specified)
%       saveIntermediate   <bool>  Boolean value to save the intermediate files
%                                  and a processing log in the output (1) or not (0) 
%                                  (Optional. Default: 0, not saving.)
%
%  Output arguments -optional:
%       output             Cell array containing dicom filename, BreastArea(sqcm),
%                          DenseArea(sqcm), and BreastDensity(%).
%  Usage:
%
%   >> libra(input, outDir, saveIntermediate) to process all
%	images in iDir and save output in outDir. Depending on value in
%	saveIntermediate, the intermediate files could be kept in outDir.
%
%   while using compiled executable in unix, run
%   $ libra /path/to/input/ /path/to/output/
%
%   in Windows, run
%   > libra C:\path\to\input\ C:\path\to\output\
%
%
%  $Rev: 603 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2016-10-20 16:58:24 -0400 (Thu, 20 Oct 2016) $:    Date of last commit
%
%  Contact:
%     CBIG Group <software at cbica.upenn.edu>

libraVersion = libra_version;
output={}; %preassign here to avoid issues downstream
runlibra0=1;exit_flag=0; %extra flags to allow user to choose to perform additional analyses after a pass
while runlibra0
    if nargin==0
        %%%%%%%%%%%%%%%
        %  Graphical  %
        %%%%%%%%%%%%%%%   
        % If no input arguments given (i.e., double click .exe) use a button driven interface
        % Construct a questdlg with three options - for selecting input files
        choice = questdlg({'Welcome to LIBRA!';'Please choose from the following options:'}, ...
            ['LIBRA-' libraVersion], ...
            'Process a Single DICOM','Batch Processing','Exit','Exit');
        % Handle response
        switch choice
            case 'Process a Single DICOM'
                questdlg({'Please select the DICOM image you wish to analyze.'}, ...
                    ['LIBRA-' libraVersion], ...
                    'OK','OK');
                [FILENAME, PATHNAME]=uigetfile('*','Select a Mammographic Image in DICOM Format');
                input=fullfile(PATHNAME,FILENAME);
                if isnumeric(FILENAME)
                    questdlg({'Input file required.';'Cancelling Current Analysis'}, ...
                        ['LIBRA-' libraVersion], ...
                        'OK','OK');
                    runlibra1 = 0;
                else
                    runlibra1 = 1;
                end
            case 'Batch Processing'
                questdlg({'Please select the folder containing the DICOM images you wish to analyze.'}, ...
                    ['LIBRA-' libraVersion], ...
                    'OK','OK');
                input=uigetdir(pwd,'Select Folder with Mammograms to Analyze');
                if isnumeric(input)
                    questdlg({'Input folder required.';'Cancelling Current Analysis'}, ...
                        ['LIBRA-' libraVersion], ...
                        'OK','OK');
                    runlibra1 = 0;
                else
                    runlibra1 = 1;
                end
            case 'Exit'
                questdlg({'Exiting the LIBRA Software Platform'}, ...
                    ['LIBRA-' libraVersion], ...
                    'OK','OK');
                runlibra1 = 0; exit_flag=1;
            otherwise
                questdlg({'Input selection required.';'Cancelling Current Analysis'}, ...
                    ['LIBRA-' libraVersion], ...
                    'OK','OK');
                runlibra1 = 0;
        end
    %%%%%%%%%%%%%%%
    % COMMANDLINE %
    %%%%%%%%%%%%%%%    
    elseif isdeployed && (any(strfind(input,'version')) || strcmpi(input,'-v'))
        disp(['LIBRA-' libraVersion]);
        disp('Copyright (c) 2014-2016 University of Pennsylvania. All rights reserved.')
        runlibra1=0;
    elseif isdeployed && (any(strfind(input,'help')) || strcmpi(input,'-h') || strcmpi(input,'--h'))

        help_diag={['You are using LIBRA-' libraVersion]; ...
            '   '; ...
            '   '; ...
            'Input parameters -optional: if no parameters are passed, the software runs'; ...
            'via a series of prompts for the user to specify the input (a single'; ...
            'DICOM File, or directory with mutliple for batch processing) and the'; ...
            'target output directory. If arguments are passed, the software will run'; ...
            'non-interactively and process the single file specified (or all files'; ...
            'under the directory provided) by <input>, saving the results to the'; ...
            'specified output folder <outDir>.'; ...
            '     input              <char>  Absolute path to the dicom file or'; ...
            '                                directory containing dicom images'; ...
            '                               (Optional)'; ...
            '     outDir             <char>  Absolute path to the output directory'; ...
            '                                (Required if <input> is specified)'; ...
            '     saveIntermediate   <bool>  Boolean value to save the intermediate'; ...
            '                                in the output (1) or not (0) (Optional;'; ...
            '                                Default: 0, not saving.)'; ...
            '   '; ...
            '  Output arguments -optional:'; ...
            '       output             Cell array containing dicom filename, BreastArea(sqcm),'; ...
            '                          DenseArea(sqcm), and BreastDensity(%).'; ...
            '  Usage:'; ...
            '   In unix, run'; ...
            '   $ libra /path/to/input/ /path/to/output/'; ...
            '   '; ...
            '   in Windows, run'; ...
            '   > libra C:\path\to\input\ C:\path\to\output\'; ...
            'Contact:'; ...
            '   CBIG Group <software at cbica.upenn.edu>' };
        for i=1:size(help_diag,1)
            disp(help_diag{i})
        end
        runlibra1=0;
    elseif nargin<2
        %help libra;
        disp('   You must specify both a valid input file or directory as well as an output directory');
        disp('   when using the command line interface. See "libra.exe -help" for help. Exiting LIBRA...');
        runlibra1=0;
    else
        runlibra1=1;
    end
    
    if runlibra1
        
        % Reads in the file names in each dir, ignore the first two dot and dots
        if isdir(input)
            iDir=input;
            dir_content=dir(iDir); dir_content=dir_content([dir_content.isdir]==0); %get all files in folder, drop . and ..
            raw_images=char(dir_content(:,1).name); % files to check and process
            input_dirflg=1;
        elseif exist(input,'file')
            [iDir,filename,extension]=fileparts(input);
            raw_images=strcat(filename,extension);
            input_dirflg=0;
        else
            %help libra - need to clean for case that user X's out of dialog box
            disp(['   Invalid input choice: ' input]);
            disp('Cancelling Current Analysis')
            return
        end

        
        % Select output directory.
        if ~nargin %no inputs = query the user.
            %%%%%%%%%%%%%%%
            %  Graphical  %
            %%%%%%%%%%%%%%%   
            questdlg({'Please also select the output directory you wish to place your results.'}, ...
                'LIBRA', ...
                'OK','OK');
            outDir = uigetdir(pwd,'Choose Output Directory for Results');
            if isnumeric(outDir)
                questdlg({'Output directory required.';'Cancelling Current Analysis'}, ...
                    'LIBRA', ...
                    'OK','OK');
                runlibra2=0;
            else
                runlibra2=1;
            end
        else
            runlibra2=1; %if there are atleast two nargins, then outDir was provided
        end
        
        if runlibra2
            %dialog box for save intermediate file choice, default 0 for both
            %user and commandline interfaces
            if ~nargin %have to ask user
                %%%%%%%%%%%%%%%
                %  Graphical  %
                %%%%%%%%%%%%%%%   
                choice = questdlg('Would you like to save a processing log and addtional intermediate files as well?', ...
                    'LIBRA', ...
                    'Yes','No','No');
                switch choice
                    case 'Yes'
                        saveIntermediate=1;
                    otherwise
                        saveIntermediate=0;
                end
            elseif nargin==2
                saveIntermediate=0; %set default if not provided
            end
            
            if (~exist(outDir,'dir')) %make output dir if it doesn't exist
                mkdir(outDir);
            end
            
            % Lets also create a log file for the session for all the images to
            % be run when saveIntermediate is true
            if saveIntermediate;
                diary(fullfile(outDir,['LIBRA-logfile_' datestr(now,'mmm-dd-yyyy_HH-MM-SS') '.txt']));
            end
            
            disp(['Running LIBRA-' libraVersion ' on MATLAB R' version('-release')]);
            numFiles=size(raw_images,1);
            if numFiles <= 0;
                disp(['   No files found in ', iDir]);
                disp('Cancelling Current Analysis')
                return
            else
                if input_dirflg
                    disp(['  >> ', num2str(numFiles), ' files found in ', iDir]);
                else
                    disp('  >> Analzying single file...');
                end
            end
        
            % Set output directories and files
            outImgDir=fullfile(outDir, 'Result_Images'); %Directory name in PWD
            outTxt=fullfile(outDir, 'Density.csv'); %comma-separated file to store Numeric Results
            results=zeros(numFiles,3);
            output=cell(numFiles+1,4);
            output(1,:)={'File Analyzed','BreastArea(sqcm)','DenseArea(sqcm)','BreastDensity(%)'};
            
            for i=1:numFiles
                
                %Now Select image to analyze
                dcm_fname=fullfile(iDir,raw_images(i,:)); % Full path to data
                dcm_fname=strtrim(dcm_fname);
                [~,NAME,EXT]=fileparts(dcm_fname);
                display(['  >> Processing File ', num2str(i), '/', num2str(numFiles), ': ',dcm_fname]);
                
                fail={'NA','NA','NA'};
                if isdicom(dcm_fname) % Is it a dicom?
                    try
                        %Run Density Analysis Code
                        [results(i,:)]=libra_exper(dcm_fname,outImgDir,outTxt,saveIntermediate);
                    catch err
                        disp(err.identifier);
                        disp(err.message);
                        disp('   > libra_exper failed.');
                        output(i+1,:)=[NAME,fail(1:3)];
                        continue;
                    end
                else
                    disp(['   > File ', NAME, EXT, ' is not a DICOM images to be processed. Skipping...']);
                    continue
                end
                results_cell=num2cell(results(i,:));
                output(i+1,:)=[NAME,results_cell(1:3)];
                disp('  ');  % added blank line after one image is processed.
            end
            
            if saveIntermediate;
                diary OFF %close log file
            end
            
            %only if user-driven do we want this
            if ~nargin && ispc
                choice=questdlg({'Analysis Complete!';'Open the Results Folder?'}, ...
                    'LIBRA', ...
                    'Yes','No','No');
                switch choice
                    case 'Yes'
                        winopen(outDir)
                end
            end
        end
    end
    runlibra0=0;
    if ~nargin && ~exit_flag
        choice=questdlg({'Would you like to perform additional analyses?'}, ...
            'LIBRA', ...
            'Yes','No','No');
        switch choice
            case 'Yes'
                runlibra0=1;
            otherwise
                questdlg({'Exiting the LIBRA Software Platform'}, ...
                    'LIBRA', ...
                    'OK','OK');
        end
    end
end
