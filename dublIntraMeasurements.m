function compiledIntra = dublIntraMeasurements(movies,varargin)
% DUBLINTRAMEASUREMENTS Produces a structure of population-level
% intra-kinetochore measurements over multiple experiments.
%
%    DUBLINTRAMEASUREMENTS(MOVIES,...) Provides all coordinates, inter-
%    and intra-kinetochore measurements for all sisters within all cells
%    across all experiments. The resulting structure provides the tools for
%    deriving population-level analyses for an experiment. Options are
%    available.
%
%    Options, defaults in {}:-
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. The channels between
%       which to make intra-kinetochore measurements. The direction of
%       measurements will be defined by the channel orientation in the
%       neighbourSpots section of options.
%
%    paired: 0 or {1}. Whether or not to take paired measurements, or raw
%       spot-by-spot measurements.
%
%    prevMeas: {[]} or a structure previously generated. A structure of
%       results from previous experiments to allow new experiment data to
%       be appended.
%
%    sisterStructure: {[]} or cell of arrays. A cell containing lists of
%       sister pair IDs, accompanied by its movie number, for each
%       experiment. Allows for only specific sisters to be included in the
%       data collection.
%
%
% Copyright (c) 2016 C. A. Smith

% default options
opts.channels = [1 2];
opts.paired = 1;
opts.prevMeas = [];
opts.sisterStructure = [];
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
if isempty(opts.sisterStructure)
  sisters = repmat({[]},numExpts1,1);
  allocSis = 0;
elseif ~iscell(opts.sisterStructure)
  sisters = {opts.sisterStructure};
  kitLog('Non-cell sister structure provided. Assuming only one experiment.');
  allocSis = 1;
else
  sisters = opts.sisterStructure;
  allocSis = 1;
end

%find number of movies and sisters, and ensure they match
numExpts2 = length(sisters);
if numExpts1 ~= numExpts2
  error('Have %i sister lists for %i experiments. Please provide a sister list for each experiments.',numExpts2,numExpts1)
end
numExpts = numExpts1;

%% Preprocessing output structure

process = movies{1}{1}.options.jobProcess;

% find whether any movies have paired spots
paired = 0;
if opts.paired
  for iExpt = 1:numExpts;
    for iMov = 1:length(movies{numExpts})
      paired = isfield(movies{iExpt}{iMov}.dataStruct{opts.channels(1)},'sisterList');
      if paired; break; end
    end
    if paired; break; end
  end
end

if isempty(opts.prevMeas)
    
    % make new intra-measurements structure if no previous measurements
    % provided
    allData = dublMakeIntraStructure(paired);
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
    theseSisters = sisters{iExpt};
    
    % find channel vector
    chanVect = movies{iExpt}{1}.options.neighbourSpots.channelOrientation;
    remChans = setdiff(chanVect,opts.channels);
    chanVect(chanVect==remChans) = [];
    
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
      
      % check whether the movie failed
      if ~isfield(dSinner,'failed') || dSinner.failed || ~isfield(dSouter,'failed') || dSouter.failed
        noSpot = [noSpot; iExpt movNum];
        continue
      end
        
      if paired
        
        % check whether a sisterList is present, and if it contains any sisters
        if ~isfield(dSinner,'sisterList') || isempty(dSinner.sisterList(1).trackPairs) || ~isfield(dSouter,'sisterList') || isempty(dSouter.sisterList(1).trackPairs)
          noSis = [noSis; iExpt movNum];
          continue
        end
        
        % if no sisters given, go through all sisters in movie
        if ~allocSis
            theseSisters = 1:length(theseMovies{iMov}.dataStruct{chanVect(1)}.sisterList);
        else % otherwise, get the sister list
            theseSisters = sisters{iExpt}(sisters{iExpt}(:,1)==iMov,2)';
        end

        for iSis = theseSisters
            
            % start counter for storing data
            c=1;

            % get sisterLists
            sLinner = dSinner.sisterList(iSis);
            sLouter = dSouter.sisterList(iSis);
            % check that there are sisters
            if isempty(dSinner.sisterList(1).trackPairs)
                continue
            end
            
            % get trackID and spotIDs
            trackIDs = dSinner.sisterList(1).trackPairs(iSis,1:2);
            for iTrack = 1:2
                spotIDs(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).featIndx;
            end
            
            % get basic metadata
            nFrames = size(spotIDs,1);
            pixelSize = theseMovies{iMov}.metadata.pixelSize;

            % get initCoords for these tracks
            for iTime = 1:nFrames
                if ~isnan(spotIDs(iTime,1))
                    iCinner.coords1(iTime,:) = dSinner.initCoord(iTime).allCoord(spotIDs(iTime,1),:);
                    iCouter.coords1(iTime,:) = dSouter.initCoord(iTime).allCoord(spotIDs(iTime,1),:);
                else
                    iCinner.coords1(iTime,:) = NaN(1,6);
                    iCouter.coords1(iTime,:) = NaN(1,6);
                end
                if ~isnan(spotIDs(iTime,2))
                    iCinner.coords2(iTime,:) = dSinner.initCoord(iTime).allCoord(spotIDs(iTime,2),:);
                    iCouter.coords2(iTime,:) = dSouter.initCoord(iTime).allCoord(spotIDs(iTime,2),:);
                else
                    iCinner.coords2(iTime,:) = NaN(1,6);
                    iCouter.coords2(iTime,:) = NaN(1,6);
                end
            end

            %% Kinetochore positions

            % get microscope coordinates of each spot
            micrCoordsInner = [iCinner.coords1(:,1:3) iCinner.coords2(:,1:3)];
            micrCoordsOuter = [iCouter.coords1(:,1:3) iCouter.coords2(:,1:3)];
            micrCoord1_x = [iCinner.coords1(:,1) iCouter.coords1(:,1)];
            micrCoord1_y = [iCinner.coords1(:,2) iCouter.coords1(:,2)];
            micrCoord1_z = [iCinner.coords1(:,3) iCouter.coords1(:,3)];
            micrCoord2_x = [iCinner.coords2(:,1) iCouter.coords2(:,1)];
            micrCoord2_y = [iCinner.coords2(:,2) iCouter.coords2(:,2)];
            micrCoord2_z = [iCinner.coords2(:,3) iCouter.coords2(:,3)];
            % put data into string format
            newData(c,:) = {'microscope.coords.x',[micrCoord1_x micrCoord2_x]}; c=c+1;
            newData(c,:) = {'microscope.coords.y',[micrCoord1_y micrCoord2_y]}; c=c+1;
            newData(c,:) = {'microscope.coords.z',[micrCoord1_z micrCoord2_z]}; c=c+1;
            
            if isfield(dSinner,'planeFit') && ~isempty(dSinner.planeFit)
              plate = 1;
            else
              plate = 0;
            end
            if plate
              % get plate coordinates of each spot
              plateCoordsInner = [sLinner.coords1(:,1:3) sLinner.coords2(:,1:3)];
              plateCoordsOuter = [sLouter.coords1(:,1:3) sLouter.coords2(:,1:3)];
              plateCoord1_x = [sLinner.coords1(:,1) sLouter.coords1(:,1)];
              plateCoord1_y = [sLinner.coords1(:,2) sLouter.coords1(:,2)];
              plateCoord1_z = [sLinner.coords1(:,3) sLouter.coords1(:,3)];
              plateCoord2_x = [sLinner.coords2(:,1) sLouter.coords2(:,1)];
              plateCoord2_y = [sLinner.coords2(:,2) sLouter.coords2(:,2)];
              plateCoord2_z = [sLinner.coords2(:,3) sLouter.coords2(:,3)];
              % put data into string format
              newData(c,:) = {'plate.coords.x',[plateCoord1_x plateCoord2_x]}; c=c+1;
              newData(c,:) = {'plate.coords.y',[plateCoord1_y plateCoord2_y]}; c=c+1;
              newData(c,:) = {'plate.coords.z',[plateCoord1_z plateCoord2_z]}; c=c+1;
            end
            
            %% Inter- and intra-kinetochore measurements
            
            micrData = pairedMeasurements(micrCoordsInner,micrCoordsOuter);
            % put data into string format
            newData(c,:) = {'microscope.sisSep.x',micrData.sisSep_x}; c=c+1;
            newData(c,:) = {'microscope.sisSep.y',micrData.sisSep_y}; c=c+1;
            newData(c,:) = {'microscope.sisSep.z',micrData.sisSep_z}; c=c+1;
            newData(c,:) = {'microscope.sisSep.twoD',micrData.sisSep_2D}; c=c+1;
            newData(c,:) = {'microscope.sisSep.threeD',micrData.sisSep_3D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.x.all',micrData.delta_x}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.y.all',micrData.delta_y}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.z.all',micrData.delta_z}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.oneD',micrData.delta_1D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.twoD.all',micrData.delta_2D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.threeD.all',micrData.delta_3D}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.y.all',micrData.swivel_y}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.z.all',micrData.swivel_z}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.threeD.all',micrData.swivel_3D}; c=c+1;
            if isfield(micrData,'swivel_kMT')
              newData(c,:) = {'microscope.raw.swivel.kMT',micrData.swivel_kMT}; c=c+1;
            else
              newData(c,:) = {'microscope.raw.swivel.kMT',[]}; c=c+1;
            end
            
            if plate
              plateData = pairedMeasurements(plateCoordsInner,plateCoordsOuter,1);
              % put data into string format
              newData(c,:) = {'plate.sisSep.x',plateData.sisSep_x}; c=c+1;
              newData(c,:) = {'plate.sisSep.y',plateData.sisSep_y}; c=c+1;
              newData(c,:) = {'plate.sisSep.z',plateData.sisSep_z}; c=c+1;
              newData(c,:) = {'plate.sisSep.twoD',plateData.sisSep_2D}; c=c+1;
              newData(c,:) = {'plate.sisSep.threeD',plateData.sisSep_3D}; c=c+1;
              newData(c,:) = {'plate.raw.delta.x.all',plateData.delta_x}; c=c+1;
              newData(c,:) = {'plate.raw.delta.y.all',plateData.delta_y}; c=c+1;
              newData(c,:) = {'plate.raw.delta.z.all',plateData.delta_z}; c=c+1;
              newData(c,:) = {'plate.raw.delta.oneD',plateData.delta_1D}; c=c+1;
              newData(c,:) = {'plate.raw.delta.twoD.all',plateData.delta_2D}; c=c+1;
              newData(c,:) = {'plate.raw.delta.threeD.all',plateData.delta_3D}; c=c+1;
              newData(c,:) = {'plate.raw.swivel.y.all',plateData.swivel_y}; c=c+1;
              newData(c,:) = {'plate.raw.swivel.z.all',plateData.swivel_z}; c=c+1;
              newData(c,:) = {'plate.raw.swivel.threeD.all',plateData.swivel_3D}; c=c+1;
              if isfield(plateData,'swivel_kMT')
                newData(c,:) = {'plate.raw.swivel.kMT',plateData.swivel_kMT}; c=c+1;
              else
                newData(c,:) = {'plate.raw.swivel.kMT',[]}; c=c+1;
              end
              newData(c,:) = {'plate.raw.twist.y',plateData.twist_y}; c=c+1;
              newData(c,:) = {'plate.raw.twist.z',plateData.twist_z}; c=c+1;
              newData(c,:) = {'plate.raw.twist.threeD',plateData.twist_3D}; c=c+1;
              newData(c,:) = {'plate.sisterCentreSpeed',plateData.sisCentreSpeed}; c=c+1;

              % get plate thickness measurements
              sisCentre_x = [];
              if length(dSinner.sisterList) == 1
                plateThickness = nan(1,nFrames);
              else
                for jSis = 1:length(dSinner.sisterList)
                  sisCentre_x = [sisCentre_x nanmean([dSinner.sisterList(jSis).coords1(:,1) dSinner.sisterList(jSis).coords2(:,1)],2)];
                end
                plateThickness = nanstd(sisCentre_x,[],2);
              end
              % put data into string format
              newData(c,:) = {'plate.plateThickness',plateThickness}; c=c+1;
            
            end

            %% Quality control
            
            % find which data satisfies the z-depth requirement
            satisfies = +(abs(micrData.delta_z)<0.5*pixelSize(3));
            satisfies(satisfies==0) = NaN;
            % put all data into string format
            newData(c,:) = {'microscope.depthFilter.delta.x.all',micrData.delta_x.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.y.all',micrData.delta_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.z.all',micrData.delta_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.oneD',micrData.delta_1D.*prod(satisfies,2)}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.twoD.all',micrData.delta_2D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.threeD.all',micrData.delta_3D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.y.all',micrData.swivel_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.z.all',micrData.swivel_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.threeD.all',micrData.swivel_3D.*satisfies}; c=c+1;
            if isfield(micrData,'swivel_kMT')
              newData(c,:) = {'microscope.depthFilter.swivel.kMT',micrData.swivel_kMT.*satisfies}; c=c+1;
            else
              newData(c,:) = {'microscope.depthFilter.swivel.kMT',[]}; c=c+1;
            end
            if plate
              newData(c,:) = {'plate.depthFilter.delta.x.all',plateData.delta_x.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.y.all',plateData.delta_y.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.z.all',plateData.delta_z.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.oneD',plateData.delta_1D.*prod(satisfies,2)}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.twoD.all',plateData.delta_2D.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.threeD.all',plateData.delta_3D.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.swivel.y.all',plateData.swivel_y.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.swivel.z.all',plateData.swivel_z.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.swivel.threeD.all',plateData.swivel_3D.*satisfies}; c=c+1;
              if isfield(plateData,'swivel_kMT')
                newData(c,:) = {'plate.depthFilter.swivel.kMT',plateData.swivel_kMT.*satisfies}; c=c+1;
              else
                newData(c,:) = {'plate.depthFilter.swivel.kMT',[]}; c=c+1;
              end
              newData(c,:) = {'plate.depthFilter.twist.y',plateData.twist_y.*nanmax(satisfies,[],2)}; c=c+1;
              newData(c,:) = {'plate.depthFilter.twist.z',plateData.twist_z.*nanmax(satisfies,[],2)}; c=c+1;
              newData(c,:) = {'plate.depthFilter.twist.threeD',plateData.twist_3D.*nanmax(satisfies,[],2)}; c=c+1;
            end
            
            %% Directional information
            
            % get direction of movement
            direc = [];
            for iTrack = 1:2
              if ~isempty(dSinner.trackList(trackIDs(iTrack)).direction)
                direc(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).direction;
              end
            end
            direc(end+1,:) = NaN;
            
            % get potential switch events (i.e. individual timepoints between P and AP)
            switchBuffer = 4;
            switchEvent = zeros(size(direc));
            switchDirec = [diff(direc); NaN NaN];
            for iPoint = 1:2:size(switchEvent,1)-(switchBuffer+1)
                for jSis = 1:2
                    tempDirec = switchDirec(iPoint:iPoint+(switchBuffer-1),jSis);
                    idx = find(tempDirec(2:(switchBuffer-1))==0);
                    if max(abs(tempDirec))==1 && abs(sum(tempDirec))>1 && ~isempty(idx)
                        switchEvent(iPoint+idx(1):iPoint+idx(end),jSis) = 1;
                    end
                end
            end
            % calculate directional information
            direc_P = +(direc==1);               direc_P(direc_P==0) = NaN;
            direc_AP = +(direc==-1);             direc_AP(direc_AP==0) = NaN;
            direc_S = +(switchEvent==1);         direc_S(direc_S==0) = NaN;
            direc_N = +((direc+switchEvent)==0); direc_N(direc_N==0) = NaN;
            % put data into string format
            dirLabel = {'P','AP','S','N'};
            for iDir = 1:4
                eval(['newData(c,:) = {''direction.' dirLabel{iDir} ''',direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.delta.x.' dirLabel{iDir} ''',micrData.delta_x.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.delta.y.' dirLabel{iDir} ''',micrData.delta_y.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.delta.z.' dirLabel{iDir} ''',micrData.delta_z.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.delta.twoD.' dirLabel{iDir} ''',micrData.delta_2D.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.delta.threeD.' dirLabel{iDir} ''',micrData.delta_3D.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.swivel.y.' dirLabel{iDir} ''',micrData.swivel_y.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.swivel.z.' dirLabel{iDir} ''',micrData.swivel_z.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.raw.swivel.threeD.' dirLabel{iDir} ''',micrData.swivel_3D.*direc_' dirLabel{iDir} '}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.delta.x.' dirLabel{iDir} ''',micrData.delta_x.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.delta.y.' dirLabel{iDir} ''',micrData.delta_y.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.delta.z.' dirLabel{iDir} ''',micrData.delta_z.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.delta.twoD.' dirLabel{iDir} ''',micrData.delta_2D.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.delta.threeD.' dirLabel{iDir} ''',micrData.delta_3D.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.swivel.y.' dirLabel{iDir} ''',micrData.swivel_y.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.swivel.z.' dirLabel{iDir} ''',micrData.swivel_z.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                eval(['newData(c,:) = {''microscope.depthFilter.swivel.threeD.' dirLabel{iDir} ''',micrData.swivel_3D.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                if plate
                  eval(['newData(c,:) = {''plate.raw.delta.x.' dirLabel{iDir} ''',plateData.delta_x.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.delta.y.' dirLabel{iDir} ''',plateData.delta_y.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.delta.z.' dirLabel{iDir} ''',plateData.delta_z.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.delta.twoD.' dirLabel{iDir} ''',plateData.delta_2D.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.delta.threeD.' dirLabel{iDir} ''',plateData.delta_3D.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.swivel.y.' dirLabel{iDir} ''',plateData.swivel_y.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.swivel.z.' dirLabel{iDir} ''',plateData.swivel_z.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.raw.swivel.threeD.' dirLabel{iDir} ''',plateData.swivel_3D.*direc_' dirLabel{iDir} '}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.delta.x.' dirLabel{iDir} ''',plateData.delta_x.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.delta.y.' dirLabel{iDir} ''',plateData.delta_y.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.delta.z.' dirLabel{iDir} ''',plateData.delta_z.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.delta.twoD.' dirLabel{iDir} ''',plateData.delta_2D.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.delta.threeD.' dirLabel{iDir} ''',plateData.delta_3D.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.swivel.y.' dirLabel{iDir} ''',plateData.swivel_y.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.swivel.z.' dirLabel{iDir} ''',plateData.swivel_z.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                  eval(['newData(c,:) = {''plate.depthFilter.swivel.threeD.' dirLabel{iDir} ''',plateData.swivel_3D.*direc_' dirLabel{iDir} '.*satisfies}; c=c+1;']);
                end
            end
        
            % compile new data with original
            allData = combineStrForms(allData,newData);
            
            % clear some data to ensure no overlap on next loop
            clear spotIDs iCinner iCouter newData direc
        
        end % sisters
        
      else
        
          % start counter for storing data
          c=1;
          
          % get dataStructs
          dSinner = theseMovies{iMov}.dataStruct{chanVect(1)};
          dSouter = theseMovies{iMov}.dataStruct{chanVect(2)};
          % get initCoords
          iCinner = dSinner.initCoord;
          iCouter = dSouter.initCoord;
          
          % get basic metadata
          nFrames = theseMovies{iMov}.metadata.nFrames;
          pixelSize = theseMovies{iMov}.metadata.pixelSize;
          
          %% Kinetochore positions
          
          % predefine various variables
          micrCoordsInner = []; micrCoordsOuter = [];
          micrCoord_x = []; micrCoord_y = []; micrCoord_z = [];
          
          for iFrame = 1:nFrames
            % get microscope coordinates of each spot
            micrCoordsInner = [micrCoordsInner; iCinner(iFrame).allCoord(:,1:3)];
            micrCoordsOuter = [micrCoordsOuter; iCouter(iFrame).allCoord(:,1:3)];
            micrCoord_x = [micrCoord_x; [iCinner(iFrame).allCoord(:,1) iCouter(iFrame).allCoord(:,1)]];
            micrCoord_y = [micrCoord_y; [iCinner(iFrame).allCoord(:,2) iCouter(iFrame).allCoord(:,2)]];
            micrCoord_z = [micrCoord_z; [iCinner(iFrame).allCoord(:,3) iCouter(iFrame).allCoord(:,3)]];
          end
           
          % put data into string format
          newData(c,:) = {'microscope.coords.x',micrCoord_x}; c=c+1;
          newData(c,:) = {'microscope.coords.y',micrCoord_y}; c=c+1;
          newData(c,:) = {'microscope.coords.z',micrCoord_z}; c=c+1;
            
          if isfield(dSinner,'planeFit') && ~isempty(dSinner.planeFit)
            plate = 1;
            pFinner = dSinner.planeFit; pFouter = dSouter.planeFit;
          else
            plate = 0;
          end
          if plate
            
            % predefine various variables
            plateCoordsInner = []; plateCoordsOuter = [];
            plateCoord_x = []; plateCoord_y = []; plateCoord_z = [];
          
            for iFrame = 1:nFrames
              % get microscope coordinates of each spot
              plateCoordsInner = [plateCoordsInner; pFinner(iFrame).rotatedCoord(:,1:3)];
              plateCoordsOuter = [plateCoordsOuter; pFouter(iFrame).rotatedCoord(:,1:3)];
              plateCoord_x = [plateCoord_x; [pFinner(iFrame).rotatedCoord(:,1) pFouter(iFrame).rotatedCoord(:,1)]];
              plateCoord_y = [plateCoord_y; [pFinner(iFrame).rotatedCoord(:,2) pFouter(iFrame).rotatedCoord(:,2)]];
              plateCoord_z = [plateCoord_z; [pFinner(iFrame).rotatedCoord(:,3) pFouter(iFrame).rotatedCoord(:,3)]];
              % put data into string format
              newData(c,:) = {'plate.coords.x',plateCoord_x}; c=c+1;
              newData(c,:) = {'plate.coords.y',plateCoord_y}; c=c+1;
              newData(c,:) = {'plate.coords.z',plateCoord_z}; c=c+1;
            end
          end
          
          %% Inter- and intra-kinetochore measurements
            
          micrData = soloMeasurements(micrCoordsInner,micrCoordsOuter);
          % put data into string format
          newData(c,:) = {'microscope.raw.delta.x.all',micrData.delta_x}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.y.all',micrData.delta_y}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.z.all',micrData.delta_z}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.twoD.all',micrData.delta_2D}; c=c+1;
          newData(c,:) = {'microscope.raw.delta.threeD.all',micrData.delta_3D}; c=c+1;
          
          if plate
            plateData = soloMeasurements(plateCoordsInner,plateCoordsOuter);
            % put data into string format
            newData(c,:) = {'plate.raw.delta.x.all',plateData.delta_x}; c=c+1;
            newData(c,:) = {'plate.raw.delta.y.all',plateData.delta_y}; c=c+1;
            newData(c,:) = {'plate.raw.delta.z.all',plateData.delta_z}; c=c+1;
            newData(c,:) = {'plate.raw.delta.twoD.all',plateData.delta_2D}; c=c+1;
            newData(c,:) = {'plate.raw.delta.threeD.all',plateData.delta_3D}; c=c+1;
          end
          
          %% Quality control
            
            % find which data satisfies the z-depth requirement
            satisfies = +(abs(micrData.delta_z)<0.5*pixelSize(3));
            satisfies(satisfies==0) = NaN;
            % put all data into string format
            newData(c,:) = {'microscope.depthFilter.delta.x.all',micrData.delta_x.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.y.all',micrData.delta_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.z.all',micrData.delta_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.twoD.all',micrData.delta_2D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.threeD.all',micrData.delta_3D.*satisfies}; c=c+1;
            
            if plate
              newData(c,:) = {'plate.depthFilter.delta.x.all',plateData.delta_x.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.y.all',plateData.delta_y.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.z.all',plateData.delta_z.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.twoD.all',plateData.delta_2D.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.threeD.all',plateData.delta_3D.*satisfies}; c=c+1;
            end
          
          % compile new data with original
          allData = combineStrForms(allData,newData);  
            
      end % paired
    end % movies     
end % expts

%% Save results to structure

compiledIntra = strForm2struct(allData);

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

function measurements = pairedMeasurements(coordsInner,coordsOuter,plate)
    
    if nargin<3 || isempty(plate)
        plate = 0;
    end
        
    %% Inter-kinetochore: separation, twist, and velocities
    
    % coordinate-specific sister separation (using inner-kinetochore)
    sisSep_xyz = diff(reshape(coordsInner,size(coordsInner,1),3,2),[],3);
    measurements.sisSep_x = sisSep_xyz(:,1);
    measurements.sisSep_y = sisSep_xyz(:,2);
    measurements.sisSep_z = sisSep_xyz(:,3);

    % 3D sister separation
    sisSep_3D = sqrt(sum(sisSep_xyz.^2,2));
    measurements.sisSep_3D = sisSep_3D;

    % 2D sister separation
    sisSep_2D = sqrt(sum(sisSep_xyz(:,1:2).^2,2));
    measurements.sisSep_2D = sisSep_2D;

    if plate

        % coordinate-specific twist
        twist_y = sisSep_xyz(:,2)./sisSep_xyz(:,1);
        twist_y = atand(twist_y);
        twist_z = sisSep_xyz(:,3)./sisSep_xyz(:,1);
        twist_z = atand(twist_z);
        measurements.twist_y = twist_y;
        measurements.twist_z = twist_z;

        % 3D twist (dot product with the x-axis with length 1)
        xAxis = repmat([1 0 0],size(sisSep_xyz,1),1);
        twist_3D = dot(sisSep_xyz,xAxis,2);
        twist_3D = twist_3D./(sisSep_3D(:,1));
        twist_3D = acosd(twist_3D);
        twist_3D(twist_3D>90) = 180-twist_3D(twist_3D>90);
        measurements.twist_3D = twist_3D;

        % sister centre velocities in x-coordinate
        sisCentre_x = sum(coordsInner(:,[1 4]),2)/2;
        sisCentreSpeed_x = [diff(sisCentre_x); NaN];
        measurements.sisCentreSpeed = sisCentreSpeed_x;

    end

    %% Intra-kinetochore delta vector
    
    % coordinate-specific delta
    delta1_xyz = coordsOuter(:,1:3) - coordsInner(:,1:3);
    delta2_xyz = coordsOuter(:,4:6) - coordsInner(:,4:6);
    measurements.delta_x = [delta1_xyz(:,1) delta2_xyz(:,1)];
    measurements.delta_y = [delta1_xyz(:,2) delta2_xyz(:,2)];
    measurements.delta_z = [delta1_xyz(:,3) delta2_xyz(:,3)];
    
    % 3D delta
    delta1_3D = sqrt(sum(delta1_xyz.^2,2));
    delta2_3D = sqrt(sum(delta2_xyz.^2,2));
    measurements.delta_3D = [delta1_3D delta2_3D];
    
    % 2D delta
    delta1_2D = sqrt(sum(delta1_xyz(:,1:2).^2,2));
    delta2_2D = sqrt(sum(delta2_xyz(:,1:2).^2,2));
    measurements.delta_2D = [delta1_2D delta2_2D];
    
    % 1D delta
    outerSisSep_xyz = diff(reshape(coordsOuter,size(coordsOuter,1),3,2),[],3);
    outerSisSep_2D = sqrt(sum(outerSisSep_xyz(:,1:2).^2,2));
    delta_1D = (outerSisSep_2D-sisSep_2D)/2;
    measurements.delta_1D = delta_1D;
    
    %% Intra-kinetochore swivel
    
    % y swivel
    dy = sisSep_2D;
    del = delta1_2D;
    eps = sqrt(sum((coordsInner(:,4:5)-coordsOuter(:,1:2)).^2,2));
    swivel1_y = (dy.^2+del.^2-eps.^2)./(2*dy.*del);
    swivel1_y = 180 - acosd(swivel1_y);
    swivel1_y = swivel1_y.*sign(delta1_xyz(:,2));

    del = delta2_2D;
    eps = sqrt(sum((coordsOuter(:,4:5)-coordsInner(:,1:2)).^2,2));
    swivel2_y = (dy.^2+del.^2-eps.^2)./(2*dy.*del);
    swivel2_y = 180 - acosd(swivel2_y);
    swivel2_y = swivel2_y.*sign(delta2_xyz(:,2));

    measurements.swivel_y = [swivel1_y swivel2_y];

    % z swivel
    dz = sqrt(sum(sisSep_xyz(:,[1 3]).^2,2));
    del = sqrt(sum((coordsInner(:,[1 3])-coordsOuter(:,[1 3])).^2,2));
    eps = sqrt(sum((coordsOuter(:,[1 3])-coordsInner(:,[4 6])).^2,2));
    swivel1_z = (dz.^2+del.^2-eps.^2)./(2*dz.*del);
    swivel1_z = 180 - acosd(swivel1_z);
    swivel1_z = swivel1_z.*sign(delta1_xyz(:,3));

    del = sqrt(sum((coordsInner(:,[4 6])-coordsOuter(:,[4 6])).^2,2));
    eps = sqrt(sum((coordsOuter(:,[4 6])-coordsInner(:,[1 3])).^2,2));
    swivel2_z = (dz.^2+del.^2-eps.^2)./(2*dz.*del);
    swivel2_z = 180 - acosd(swivel2_z);
    swivel2_z = swivel2_z.*sign(delta2_xyz(:,3));

    measurements.swivel_z = [swivel1_z swivel2_z];

    % 3D swivel
    swivel1_3D = dot(-sisSep_xyz,delta1_xyz,2);
    swivel1_3D = swivel1_3D./(sisSep_3D(:,1).*delta1_3D(:,1));
    swivel1_3D = acosd(swivel1_3D); 
    swivel2_3D = dot(sisSep_xyz,delta2_xyz,2);
    swivel2_3D = swivel2_3D./(sisSep_3D(:,1).*delta2_3D(:,1));
    swivel2_3D = acosd(swivel2_3D);
    measurements.swivel_3D = [swivel1_3D swivel2_3D];    

end %pairedMeasurements subfunction

function measurements = soloMeasurements(coordsInner,coordsOuter)
    
    % coordinate-specific delta
    delta_xyz = coordsOuter(:,1:3) - coordsInner(:,1:3);
    measurements.delta_x = delta_xyz(:,1);
    measurements.delta_y = delta_xyz(:,2);
    measurements.delta_z = delta_xyz(:,3);
    
    % 3D delta
    delta_3D = sqrt(sum(delta_xyz.^2,2));
    measurements.delta_3D = delta_3D;
    
    % 2D delta
    delta_2D = sqrt(sum(delta_xyz(:,1:2).^2,2));
    measurements.delta_2D = delta_2D;

end %soloMeasurements subfunction

end


