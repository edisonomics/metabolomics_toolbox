function [spectra, outputPath] = load_ft1_files(folderPath)
    % Load all .ft files in the current directory. If any arguments are given,
    % a multi-selection dialog is presented. Outputs to the workspace in a
    % variable called 'spectra'
    % pickFiles: the loading file directory
    % Files will be loaded in numerical order if they start with a digit
    % separated by an underscore.
    %% this is from initial toolbox
    %  YW add some documents 10/10/2018
    if nargin==0
        filenames = dir('*.ft1');
        filenames = {filenames.name};
        %filenames = sortFilenames(filenames); %% this doesn't work, MJ implemented 3JAN2018:
            %filenums = regexp(filenames,'\d+','match');
            %[~,i] = sort(cellfun(@str2num,[filenums{:}]));
            %filenames = filenames(i);
        ftlist = cellfun(@(x) ([cd filesep x]),filenames,'uniformoutput',0);
    elseif nargin==1
        if isnumeric(folderPath)
                error('folderPath should be a string or char array, not a number.');
        end
         % Check if folderPath is an existing directory
        if ~isfolder(folderPath) && ~isfile(folderPath)
            error('The provided folderPath does not exist or is not a directory.');
        end
    
        if isfile(folderPath)
            outputPath = folderPath;
            [~, name, ext] = fileparts(folderPath);
            filename = [name, ext]; 
            spectra=pipe2matlab(folderPath);
            spectra.FileName = filename;
        end
    
        if isfolder(folderPath)
            [filenames,PathName] = uigetfile(fullfile(folderPath,  '*.ft1'),'Select files...','MultiSelect','on');
            outputPath = PathName;
            if isequal(filenames, 0)
                spectra = [];
                return;
            end
            if ischar(filenames)
                filenames = {filenames};
            end
            ftlist = cellfun(@(x) ([PathName x]),filenames,'uniformoutput',0);
    
            for i=1:length(ftlist)
                spectra(i)=pipe2matlab(ftlist{i});
            end
            
            for i = 1:length(ftlist)
                spectra(i).FileName = filenames{i};
            end
            % Display summary
            if length(filenames)>=10
                filenames(1:10)';
            else
                filenames';
            end
            % assignin('caller', 'spectra', spectra)
        end
    end



end