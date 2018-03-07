function compiledInts = dublIntensityMeasurements(movies,varargin)
% DUBLINTENSITYMEASUREMENTS Produces a structure of population-level
% kinetochore intensity measurements over multiple experiments.
%
%    DUBLINTENSITYMEASUREMENTS(MOVIES,...) Provides all kinetochore
%    intensity measurements for all sisters within all cells across all
%    experiments. The resulting structure provides the tools for deriving
%    population-level intensity analyses for an experiment. Options are 
%    available.
%
%    Options, defaults in {}:-
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. The channels between
%       which to make intensity measurements.
%
%    refMarker: {'self'}, 'inner' or 'outer'. The marker around which
%       intensity measurements are measured. 'self' means intensity
%       measurements are measured about it's own spot centre, while 'inner'
%       and 'outer' mean intensity measurements are measured about the
%       inner or outer kinetochore marker, respectively.
%
%    paired: 0 or {1}. Whether or not to take paired measurements, or raw
%       spot-by-spot measurements.
%
%    prevMeas: {[]} or a structure previously generated. A structure of
%       results from previous experiments to allow new experiment data to
%       be appended.
%
%    spotSelection: {[]} or output from kitSelectData. A structure
%       containing a selection of either sister pair or track IDs per
%       movie, for each experiment. Allows for only specific sisters or
%       spots to be included in the data collection.
%
%
% Copyright (c) 2018 C. A. Smith

% default options
opts.channels = [1 2];
opts.refMarker = 'self';
opts.paired = 1;
opts.prevMeas = [];
opts.spotSelection = [];
% user options
opts = processOptions(opts,varargin{:});

%% Pre-processing input structure

%check structure of movies
if ~iscell(movies{1})
    movies = {movies};
    kitLog('Movie structure provided implies only one experiment. Assuming only one experiment.');
end
%find number of movies
numExpts1 = length(movies);

%process input so that all structs are in cell format
if isempty(opts.spotSelection)
  subset = repmat({[]},numExpts1,1);
  selType = 0;
elseif isstruct(opts.spotSelection) && isfield(opts.spotSelection,'dataType')
  subset = opts.spotSelection.selection;
  switch opts.spotSelection.dataType
    case 'spots' %tracks
      selType = 1;
    case 'sisters'
      selType = 2;
    case 'initCoord'
      selType = 3;
  end
else
  kitLog('Provided spotSelection structure was not derived from kitSelectData. No selection will be imposed.');
  subset = repmat({[]},numExpts1,1);
  selType = 0;
end

%find number of movies and sisters, and ensure they match
numExpts2 = length(subset);
if numExpts1 ~= numExpts2
  error('Have %i spot selections for %i experiments. Please provide spot selection for each experiment.',numExpts2,numExpts1)
end
numExpts = numExpts1;

%% Pre-processing output structure

% find whether any movies have paired spots
paired = 0;
if opts.paired
  for iExpt = 1:numExpts
    for iMov = 1:length(movies{iExpt})
      paired = isfield(movies{iExpt}{iMov}.dataStruct{opts.channels(1)},'sisterList');
      if paired; break; end
    end
    if paired; break; end
  end
end
if selType==3
    if paired
        subset = opts.spotSelection.selection;
        selType = 1;
    else
        subset = opts.spotSelection.rawSelection;
    end
end
    
%check whether or not there is intensity information
ints = 0;
for iExpt = 1:numExpts
  for iMov = 1:length(movies{iExpt})
    ints = isfield(movies{iExpt}{iMov}.dataStruct{opts.channels(1)},'spotInt');
    if ints; break; end
  end
  if ints; break; end  
end
if ints && size(movies{1}{1}.dataStruct{opts.channels(1)}.spotInt.intensity,2)==1 
    opts.intRefMarker = 'self';
    kitLog('Movies do not contain intensity information relative to each other channel. Converting to self-measurement.');
end

if isempty(opts.prevMeas)
    
    % make new intra-measurements structure if no previous measurements
    % provided
    allData = dublMakeIntensityStructure();
    allData = struct2strForm(allData);
    
else
    % get all old data
    allData = struct2strForm(opts.prevMeas);
    
end

% predesignate error arrays
noDS = [];
noSpot = [];
noSis = [];

%% Compiling measurements

for iExpt = 1:numExpts
    
    % get movie and sister list
    theseMovies = movies{iExpt};
    iSubset = subset{iExpt};
    
    % find channel vector
    chanVect = movies{iExpt}{1}.options.neighbourSpots.channelOrientation;
    chanVect = intersect(chanVect,opts.channels,'stable');
    
    % find spot reference channel
    refChan = movies{iExpt}{1}.options.coordSystemChannel;
    
    for iMov = 1:length(theseMovies)
      
      % get the movie index
      movNum = theseMovies{iMov}.index;
      % check whether there is data in this movie
      if ~isfield(theseMovies{iMov},'dataStruct')
        noDS = [noDS; iExpt movNum];
        continue
      end
      
      % get dataStructs
      dSinner = theseMovies{iMov}.dataStruct{chanVect(1)};
      dSouter = theseMovies{iMov}.dataStruct{chanVect(2)};
      refdS   = theseMovies{iMov}.dataStruct{refChan};
      
      % check whether the movie failed
      if ~isfield(refdS,'failed') || refdS.failed || ~isfield(refdS,'failed') || refdS.failed
        noSpot = [noSpot; iExpt movNum];
        continue
      end
      
      % get basic metadata
      nFrames = theseMovies{iMov}.metadata.nFrames;
      
      % get intensities
      sIinner = dSinner.spotInt;
      innerBg = dSinner.cellInt.back;
      sIouter = dSouter.spotInt;
      outerBg = dSouter.cellInt.back;
      
      if paired
        
        % check whether a sisterList is present, and if it contains any sisters
        if ~isfield(refdS,'sisterList') || isempty(refdS.sisterList(1).trackPairs) || ...
                ~isfield(refdS,'sisterList') || isempty(refdS.sisterList(1).trackPairs)
          noSis = [noSis; iExpt movNum];
          continue
        end
        
        % if no sisters given, go through all sisters in movie
        switch selType
            case 0 % no spot selection
                iSubset = 1:length(theseMovies{iMov}.dataStruct{chanVect(1)}.sisterList);
            case 1 % using spots/tracks
                iSubset = 1:length(theseMovies{iMov}.dataStruct{chanVect(1)}.sisterList);
                theseTracks = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
            case 2 % using sisters
                iSubset = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
        end
        
        for iSis = iSubset
            
            % construct sister pair label
            label = sprintf('%02d%02d%03d',iExpt,iMov,iSis);
            
            % start counter for storing data
            c=1;

            % get trackID and spotIDs, make spotIDs nan if deselected
            trackIDs = refdS.sisterList(1).trackPairs(iSis,1:2);
            spotIDs = nan(nFrames,2);
            for iTrack = 1:2
                if selType~=1 %none or sisters
                    spotIDs(:,iTrack) = refdS.trackList(trackIDs(iTrack)).featIndx;
                else %tracks
                    if ismember(trackIDs(iTrack),theseTracks)
                        spotIDs(:,iTrack) = refdS.trackList(trackIDs(iTrack)).featIndx;
                    else
                        trackIDs(iTrack) = NaN;
                    end
                end 
            end
            % if both spots skipped
            if all(isnan(trackIDs))
                continue
            end
            
            % get intensities if required
            switch opts.refMarker
              case 'self'
                if size(sIinner.intensity,2)==1 
                    temp = cat(2,sIinner.intensity);
                    intsInnerMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                    temp = cat(2,sIinner.intensity_max);
                    intsInnerMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                    temp = cat(2,sIouter.intensity);
                    intsOuterMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                    temp = cat(2,sIouter.intensity_max);
                    intsOuterMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                else
                    temp = cat(2,sIinner.intensity(:,chanVect(1)));
                    intsInnerMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                    temp = cat(2,sIinner.intensity_max(:,chanVect(1)));
                    intsInnerMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                    temp = cat(2,sIouter.intensity(:,chanVect(2)));
                    intsOuterMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                    temp = cat(2,sIouter.intensity_max(:,chanVect(2)));
                    intsOuterMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                end
              case 'inner'
                temp = cat(2,sIinner.intensity(:,chanVect(1)));
                intsInnerMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                temp = cat(2,sIinner.intensity_max(:,chanVect(1)));
                intsInnerMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                temp = cat(2,sIinner.intensity(:,chanVect(2)));
                intsOuterMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                temp = cat(2,sIinner.intensity_max(:,chanVect(2)));
                intsOuterMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
              case 'outer'
                temp = cat(2,sIouter.intensity(:,chanVect(1)));
                intsInnerMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                temp = cat(2,sIouter.intensity_max(:,chanVect(1)));
                intsInnerMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                temp = cat(2,sIouter.intensity(:,chanVect(2)));
                intsOuterMean = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
                temp = cat(2,sIouter.intensity_max(:,chanVect(2)));
                intsOuterMax  = temp(spotIDs(~isnan(spotIDs)),:)-innerBg;
            end
            % put data into string format
            newData(c,:) = {'intensity.mean.inner',intsInnerMean}; c=c+1;
            newData(c,:) = {'intensity.mean.outer',intsOuterMean}; c=c+1;
            newData(c,:) = {'intensity.max.inner',intsInnerMax}; c=c+1;
            newData(c,:) = {'intensity.max.outer',intsOuterMax}; c=c+1;
            newData(c,:) = {'intensity.bg.inner',repmat(innerBg,size(intsInnerMean,1),1)}; c=c+1;
            newData(c,:) = {'intensity.bg.outer',repmat(outerBg,size(intsOuterMean,1),1)}; c=c+1;
        
            % compile new data with original
            allData = combineStrForms(allData,newData);
            
            % clear some data to ensure no overlap on next loop
            clear spotIDs newData
        
        end % sisters
        
      else
          
          % start counter for storing data
          c=1;
          
          % if no sisters given, go through all sisters in movie
          switch selType
            case 0 % no spot selection
                spotIDs = 1:size(refdS.initCoord(1).allCoord,1);
            case 1 % using spots/tracks
                trackIDs = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
                spotIDs = cat(2,dSinner.trackList(trackIDs).featIndx);
            case 3 % using initCoord
                spotIDs = subset{iExpt}(subset{iExpt}(:,1)==iMov,2)';
          end
          nSpots = length(spotIDs);
          if nSpots == 0
              continue
          end
          
          % construct spot label
          for iSpot = 1:nSpots
            labels(iSpot,:) = sprintf('%02d%02d%03d',iExpt,iMov,iSpot);
          end
          % put data into string format
          newData(c,:) = {'label',labels}; c=c+1;
          
          % get intensities if required
          switch opts.refMarker
            case 'self'
                if size(sIinner.intensity,2)==1 
                    intsInnerMean = cat(2,sIinner.intensity(spotIDs,:)) - innerBg;
                    intsInnerMax  = cat(2,sIinner.intensity_max(spotIDs,:)) - innerBg;
                    intsOuterMean = cat(2,sIouter.intensity(spotIDs,:)) - outerBg;
                    intsOuterMax  = cat(2,sIouter.intensity_max(spotIDs,:)) - outerBg;
                else
                    intsInnerMean = cat(2,sIinner.intensity(spotIDs,chanVect(1))) - innerBg;
                    intsInnerMax  = cat(2,sIinner.intensity_max(spotIDs,chanVect(1))) - innerBg;
                    intsOuterMean = cat(2,sIouter.intensity(spotIDs,chanVect(2))) - outerBg;
                    intsOuterMax  = cat(2,sIouter.intensity_max(spotIDs,chanVect(2))) - outerBg;
                end
            case 'inner'
                intsInnerMean = cat(2,sIinner.intensity(spotIDs,chanVect(1))) - innerBg;
                intsInnerMax  = cat(2,sIinner.intensity_max(spotIDs,chanVect(1))) - innerBg;
                intsOuterMean = cat(2,sIinner.intensity(spotIDs,chanVect(2))) - outerBg;
                intsOuterMax  = cat(2,sIinner.intensity_max(spotIDs,chanVect(2))) - outerBg;
            case 'outer'
                intsInnerMean = cat(2,sIouter.intensity(spotIDs,chanVect(1))) - innerBg;
                intsInnerMax  = cat(2,sIouter.intensity_max(spotIDs,chanVect(1))) - innerBg;
                intsOuterMean = cat(2,sIouter.intensity(spotIDs,chanVect(2))) - outerBg;
                intsOuterMax  = cat(2,sIouter.intensity_max(spotIDs,chanVect(2))) - outerBg;
          end
          % put data into string format
          newData(c,:) = {'intensity.mean.inner',intsInnerMean}; c=c+1;
          newData(c,:) = {'intensity.mean.outer',intsOuterMean}; c=c+1;
          newData(c,:) = {'intensity.max.inner',intsInnerMax}; c=c+1;
          newData(c,:) = {'intensity.max.outer',intsOuterMax}; c=c+1;
          newData(c,:) = {'intensity.bg.inner',repmat(innerBg,size(intsInnerMean,1),1)}; c=c+1;
          newData(c,:) = {'intensity.bg.outer',repmat(outerBg,size(intsOuterMean,1),1)}; c=c+1;
          
          % compile new data with original
          allData = combineStrForms(allData,newData);  
          clear newData spotIDs
            
      end % paired
    end % movies     
end % expts

%% Save results to structure

compiledInts = strForm2struct(allData);

%% Output any error information

if ~isempty(noDS)
  fprintf('\nThe following cells failed during spot detection:\n');
  for iCell = 1:size(noDS,1)
    fprintf('    Exp %i, Mov %i\n',noDS(iCell,1),noDS(iCell,2));
  end
end
if ~isempty(noSpot)
  fprintf('\nThe following cells found no spots:\n');
  for iCell = 1:size(noSpot,1)
    fprintf('    Exp %i, Mov %i\n',noSpot(iCell,1),noSpot(iCell,2));
  end
end
if ~isempty(noSis)
  fprintf('\nThe following cells contain no sisterList:\n');
  for iCell = 1:size(noSis,1)
    fprintf('    Exp %i, Mov %i\n',noSis(iCell,1),noSis(iCell,2));
  end
end
fprintf('\n');

