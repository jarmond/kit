function selectedData = kitSelectData(expts,varargin)
% KITSELECTDATA(EXPTS,...) Allows the user to choose specific data from a
% given set of EXPTS for use in downstream analysis tools.
%
% Copyright (c) 2017 C. A. Smith

opts.addMore = [];
opts.channel = 1;
opts.contrast = [0.1 1];
opts.dataType = 'spots'; % can also be 'sisters'
opts.lineProfile = 0; %whether to plot line profiles in x and y across spots
opts.method = 'deselect'; % can also be 'select'
opts.startMovie = 1;
opts.zProject = -1;
opts = processOptions(opts,varargin{:});

% check for previous selections, set up new structure if none given
if isempty(opts.addMore)
    selectedData.dataType = opts.dataType;
    selectedData.selection = {[]};
    preprocExpts = 0;
else
    selectedData = opts.addMore;
    preprocExpts = length(selectedData.selection);
    % ensure that dataTypes match
    if ~strcmp(selectedData.dataType,opts.dataType)
      error('Cannot add more data to the provided selection due to clash in dataType.');
    end
end

% get a filename and save directory if not already provided
[filename,filepath] = uiputfile('*.mat','Save selection file','selectedSpots.mat');

% check whether input is from one or more experiments, change accordingly
if ~iscell(expts{1})
  expts = {expts};
end
nExpts = length(expts);

% get basic information
c = opts.channel;

% get previous results, if any
allSels = selectedData.selection;

% loop over all movies
for iExpt = 1:nExpts
    
    % get the jobs for this experiment
    jobs = expts{iExpt};
    nJobs = length(jobs);
    
    kitLog('Checking experiment %i of %i...',iExpt,nExpts);
    % start progress bar
    prog = kitProgress(0);
    
    % make empty allSels cell
    allSels{iExpt+preprocExpts} = [];
    
    for iJob = opts.startMovie:nJobs
        
        % check whether there is any data contained within this movie
        if ~isfield(jobs{iJob},'dataStruct') || ~isfield(jobs{iJob}.dataStruct{c},'failed') || jobs{iJob}.dataStruct{c}.failed
            continue
        end
        % check whether there are any tracks contained within this movie
        if ~isfield(jobs{iJob}.dataStruct{c},'tracks') || jobs{iJob}.dataStruct{c}.tracks(1).tracksFeatIndxCG == 0
            continue
        end
        % get dataStruct
        dS = jobs{iJob}.dataStruct{c};
        
        switch opts.dataType
          
          case 'sisters'
            % show all sisters
            kitShowAllSisters(jobs{iJob},'channel',c,'contrast',opts.contrast,...
                 'zProject',opts.zProject);
            nData = length(dS.sisterList);
                
          case 'spots'
            % show all spots - defined using tracks
            kitShowAllSpots(jobs{iJob},'channel',c,'contrast',opts.contrast,...
                'zProject',opts.zProject, 'lineProfile',opts.lineProfile);
            nData = length(dS.trackList);
        end
        
        switch opts.method
          % provide a message to request lists of spots/sisters  
          case 'select'
              kitLog('Please list all %s to be selected:',opts.dataType);
              tempList = input('');
              % check for any sisters not in the list, provide warning and
              % remove if so
              incorrect = setdiff(tempList,1:nData);
              if ~isempty(incorrect)
                  warning('The following selected %s do not exist: %s. Will ignore.',opts.dataType,num2str(incorrect));
                  tempList = setdiff(tempList,incorrect);
              end
              
          case 'deselect'
              kitLog('Please list all %s to be ignored:',opts.dataType);
              tempList = input('');
              % check for any sisters not in the list, provide warning and
              % remove if so
              incorrect = setdiff(tempList,1:nData);
              if ~isempty(incorrect)
                  warning('The following selected %s do not exist: %s. Will ignore.',opts.dataType,num2str(incorrect));
              end
              % invert the list to make it selected rather than ignored
              tempList = setdiff(1:nData,tempList);
              
        end
        
        % process the list
        if isempty(tempList)
            continue
        else
            nData = length(tempList);
            tempStruct = ones(nData,2)*iJob;
            tempStruct(:,2) = tempList;
            % collate selections
            allSels{iExpt+preprocExpts} = [allSels{iExpt+preprocExpts}; tempStruct];
        end
        
        % update progress
        prog = kitProgress(iJob/nJobs,prog);
        
    end
    
    % store final list in output, and save
    selectedData.selection = allSels;
    save(fullfile(filepath,filename),'selectedData');
    kitLog('Updated save file up to experiment %i of %i.',iExpt,nExpts);

end

kitLog('Data selection complete.');

end


    
