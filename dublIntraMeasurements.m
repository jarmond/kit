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
%    centralise: 0 or {1}. Whether or not to adjust the outer kinetochore
%       components' position so that average delta-x, delta-y and delta-z
%       measurements are zero.
%
%    channels: {[1 2]} or pair of numbers from 1 to 3. The channels between
%       which to make intra-kinetochore measurements. The direction of
%       measurements will be defined by the channel orientation in the
%       neighbourSpots section of options.
%
%    paired: 0 or {1}. Whether or not to take paired measurements, or raw
%       spot-by-spot measurements.
%
%    spotSelection: {[]} or output from kitSelectData. A structure
%       containing a selection of either sister pair or track IDs per
%       movie, for each experiment. Allows for only specific sisters or
%       spots to be included in the data collection.
%
%    verbose: {0} or 1. Whether or not to print progress to command line.
%
%
% Copyright (c) 2017 C. A. Smith

% default options
opts.centralise = 1;
opts.channels = [1 2];
opts.paired = 1;
opts.spotSelection = [];
opts.verbose = 0;
% user options
opts = processOptions(opts,varargin{:});

%% Pre-processing input structure

% check structure of movies
if ~iscell(movies{1})
    movies = {movies};
    if opts.verbose
      kitLog('Movie structure provided implies only one experiment. Assuming only one experiment.');
    end
end

% find number of movies
nExpts = length(movies);

% check whether or not a correct spotSelection has been provided
if isempty(opts.spotSelection)
    spotSel = 0;
else
    if ~iscell(opts.spotSelection)
        opts.spotSelection = {opts.spotSelection};
    end
    if ~isfield(opts.spotSelection{1},'dataType')
        spotSel = 0;
        if opts.verbose
            kitLog('Spot selection provided not produced using kitSelectData. Providing full dataset.')
        end
    else
        spotSel = 1;
    end
    testn = length(opts.spotSelection);
    if testn ~= nExpts
        spotSel = 0;
        if opts.verbose
            kitLog('Spot selection contains only %i experiments, whereas %i are required. Providing full dataset.',nExpts,testn);
        end
    end
end

% find whether any movies have paired spots
paired = 0;
if opts.paired
  for iExpt = 1:nExpts
    for iMov = 1:length(movies{iExpt})
      paired = isfield(movies{iExpt}{iMov}.dataStruct{opts.channels(1)},'sisterList');
      if paired; break; end
    end
    if paired; break; end
  end
end

% make new intra-measurements structure
allData = dublMakeIntraStructure(paired);
allData = struct2strForm(allData);

% predesignate error arrays
noDS = [];
noSpot = [];
noSis = [];

% form experiment label
labprfx = 'expt0';
for iExpt = 1:nExpts
    if iExpt==10
        labprfx(end) = '';
    end
    exptLabs(iExpt,:) = [labprfx num2str(iExpt)];
end 

%% Compiling measurements

% set up progress and counter
prog = kitProgress(0);
totalMovs = 0;
for iExpt = 1:nExpts
    totalMovs = totalMovs+length(movies{iExpt});
end
movCount = 1;

for iExpt = 1:nExpts
    
    % get movie and sister list
    theseMovies = movies{iExpt};
    nMovs = length(theseMovies);
    
    % get the channel vector
    chanVect = movies{iExpt}{1}.options.neighbourSpots.channelOrientation;
    chanVect = intersect(chanVect,opts.channels,'stable');
    
    for iMov = 1:nMovs
      
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
      
      % get pixelSize
      pixelSize = theseMovies{iMov}.metadata.pixelSize;
      
      % get initCoords
      iCinner = dSinner.initCoord;
      iCouter = dSouter.initCoord;
      nFrames = theseMovies{iMov}.metadata.nFrames;
      nSpots = size(iCinner(1).allCoord,1);
      
      % get xyz correction for centralisation
      if opts.centralise
        deltaxyz = [];
        
        for iFrame = 1:nFrames
          deltaxyz = [deltaxyz; iCouter(iFrame).allCoord(:,1:3)-iCinner(iFrame).allCoord(:,1:3)];
        end
        mCentMeds = nanmedian(deltaxyz);
      else
        mCentMeds = [0 0 0];
      end
      
      % check whether or not this movie has a planeFit
      if isfield(dSinner,'planeFit') && ~isempty(dSinner.planeFit)
        plane = 1;
        
        % get planeFits
        pFinner = dSinner.planeFit;
        pFouter = dSouter.planeFit;
        
        % produce centralisation median structures by rotating
        pCentMeds = repmat(mCentMeds,nFrames,1);
        for iFrame = 1:nFrames
          if ~isempty(pFinner(iFrame).planeVectors)
            coordSystem = pFinner(iFrame).planeVectors;
            pCentMeds(iFrame,:) = (coordSystem\mCentMeds')';
          end
        end
        
      else
        plane = 0;
      end
      
      % get experiment label
      exptLab = exptLabs(iExpt,:);
        
      if paired
        
        % check whether a sisterList is present, and if it contains any sisters
        if ~isfield(dSinner,'sisterList') || isempty(dSinner.sisterList(1).trackPairs) || ~isfield(dSouter,'sisterList') || isempty(dSouter.sisterList(1).trackPairs)
          noSis = [noSis; iExpt movNum];
          continue
        end
        
        % get trackIDs
        trackIDs = dSinner.sisterList(1).trackPairs(:,1:2);
        
        % convert spotSelection to spotIdx if required
        if spotSel
            % find the selected indices
            keepIDs = opts.spotSelection{iExpt}(:,1)==iMov;
            keepIDs = opts.spotSelection{iExpt}(keepIDs,2);
            if strcmp(opts.spotSelection{iExpt}.dataType,'sisters')
                % if sisters, get the trackIDs
                trackIDs = trackIDs(keepIDs,:);
            else
                % convert all IDs not included to NaNs
                rmIDs = ~ismember(trackIDs,keepIDs);
                trackIDs(rmIDs) = NaN;
            end
        end
        
        % get number of sisters (two tracks per sister)
        nSisters = size(trackIDs,1);
        
        % get spotInts
        if isfield(dSinner,'spotInt')
            sIinner = dSinner.spotInt;
            innerBg = dSinner.cellInt.back;
            sIouter = dSouter.spotInt;
            outerBg = dSouter.cellInt.back;
        else
            sIinner = repmat(struct('intensity',nan(nFrames,1),'intensity_max',nan(nFrames,1)),1,nSpots);
            innerBg = nan(nFrames,1);
            sIouter = repmat(struct('intensity',nan(nFrames,1),'intensity_max',nan(nFrames,1)),1,nSpots);
            outerBg = nan(nFrames,1);
        end
        
        % Loop over sisters
        for iSis = 1:nSisters
            
            % skip if both trackIDs are NaNs 
            if sum(isnan(trackIDs(iSis,:))==2)
                continue
            end
            
            % start counter for storing data
            c=1;
            % store experiment number
            newData(c,:) = {'exptLabel',exptLab}; c=c+1;
            
            % get spotIDs, unless trackID is NaN
            spotIDs = nan(nFrames,2);
            for iTrk = 1:2
                if ~isnan(trackIDs(iSis,iTrk))
                    spotIDs(:,iTrk) = dSinner.trackList(trackIDs(iSis,iTrk)).featIndx;
                end
            end

            % get microscope coordinates, and plane coordinates if present
            mCoordsInner = nan(nFrames,6);
            mCoordsOuter = nan(nFrames,6);
            if plane
                pCoordsInner = nan(nFrames,6);
                pCoordsOuter = nan(nFrames,6);
            end
            for iFrame = 1:nFrames
                for iTrk = 1:2
                    % track one stored in (:,1:3), two in (:,4:6)
                    rng = (3*(iTrk-1)+1):3*iTrk;
                    if ~isnan(spotIDs(iFrame,1))
                        mCoordsInner(iFrame,rng) = iCinner(iFrame).allCoord(spotIDs(iFrame,iTrk),1:3);
                        mCoordsOuter(iFrame,rng) = iCouter(iFrame).allCoord(spotIDs(iFrame,iTrk),1:3) - mCentMeds;
                        if plane
                            pCoordsInner(iFrame,rng) = pFinner(iFrame).rotatedCoord(spotIDs(iFrame,iTrk),1:3);
                            pCoordsOuter(iFrame,rng) = pFouter(iFrame).rotatedCoord(spotIDs(iFrame,iTrk),1:3) - pCentMeds;
                        end
                    end
                end
            end
            
            % get coordinates in sister-sister axis
            sCoordsInner = mCoordsInner; sCoordsOuter = mCoordsOuter;
            sisSep = sCoordsInner(:,1:3) - sCoordsInner(:,4:6); % points towards first ID
            % calculate angles for rotation, in XY and XZ
            theta = sisSep(:,2)./sisSep(:,1); theta = atan(theta);
            rho = sign(sisSep(:,1)).*sisSep(:,3)./sqrt(sum(sisSep(:,1:2).^2,2)); rho = atan(rho);
            % loop over spots
            for iRot = 1:length(theta)
                % get rotation angles
                iThe = theta(iRot); iRho = rho(iRot);
                xyRotMat = [cos(iThe), -sin(iThe), 0;
                            sin(iThe),  cos(iThe), 0;
                            0,          0,         1];
                xzRotMat = [cos(iRho),  0, -sin(iRho);
                            0,          1,  0;
                            sin(iRho),  0,  cos(iRho)];
                rotMat = xyRotMat*xzRotMat;
                % rotate inner coordinates
                sCoordsInner(iRot,1:3) = sCoordsInner(iRot,1:3)*rotMat;
                sCoordsInner(iRot,4:6) = sCoordsInner(iRot,4:6)*rotMat;
                % rotate outer coordinates
                sCoordsOuter(iRot,1:3) = sCoordsOuter(iRot,1:3)*rotMat;
                sCoordsOuter(iRot,4:6) = sCoordsOuter(iRot,4:6)*rotMat;
                                           
            end

            %% Kinetochore positions

            % get microscope coordinates of each spot
            mCoords_x = [mCoordsInner(:,[1 4]) mCoordsOuter(:,[1 4])];
            mCoords_y = [mCoordsInner(:,[2 5]) mCoordsOuter(:,[2 5])];
            mCoords_z = [mCoordsInner(:,[3 6]) mCoordsOuter(:,[3 6])];
            % put data into string format
            newData(c,:) = {'microscope.coords.x',mCoords_x(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'microscope.coords.y',mCoords_y(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'microscope.coords.z',mCoords_z(:,[1 3 2 4])}; c=c+1;
            
            if plane
              % get plate coordinates of each spot
              pCoords_x = [pCoordsInner(:,[1 4]) pCoordsOuter(:,[1 4])];
              pCoords_y = [pCoordsInner(:,[2 5]) pCoordsOuter(:,[2 5])];
              pCoords_z = [pCoordsInner(:,[3 6]) pCoordsOuter(:,[3 6])];
              % put data into string format
              newData(c,:) = {'plate.coords.x',pCoords_x(:,[1 3 2 4])}; c=c+1;
              newData(c,:) = {'plate.coords.y',pCoords_y(:,[1 3 2 4])}; c=c+1;
              newData(c,:) = {'plate.coords.z',pCoords_z(:,[1 3 2 4])}; c=c+1;
            end
            
            % get microscope coordinates of each spot
            sCoords_x = [sCoordsInner(:,[1 4]) sCoordsOuter(:,[1 4])];
            sCoords_y = [sCoordsInner(:,[2 5]) sCoordsOuter(:,[2 5])];
            sCoords_z = [sCoordsInner(:,[3 6]) sCoordsOuter(:,[3 6])];
            % put data into string format
            newData(c,:) = {'sisters.coords.x',sCoords_x(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'sisters.coords.y',sCoords_y(:,[1 3 2 4])}; c=c+1;
            newData(c,:) = {'sisters.coords.z',sCoords_z(:,[1 3 2 4])}; c=c+1;
            
            %% Inter- and intra-kinetochore measurements
            
            mData = pairedMeasurements(mCoordsInner,mCoordsOuter,0);
            % put data into string format
            newData(c,:) = {'microscope.sisSep.x',mData.sisSep_x}; c=c+1;
            newData(c,:) = {'microscope.sisSep.y',mData.sisSep_y}; c=c+1;
            newData(c,:) = {'microscope.sisSep.z',mData.sisSep_z}; c=c+1;
            newData(c,:) = {'microscope.sisSep.twoD',mData.sisSep_2D}; c=c+1;
            newData(c,:) = {'microscope.sisSep.threeD',mData.sisSep_3D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.x',mData.delta_x}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.y',mData.delta_y}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.z',mData.delta_z}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.oneD',mData.delta_1D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.twoD',mData.delta_2D}; c=c+1;
            newData(c,:) = {'microscope.raw.delta.threeD',mData.delta_3D}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.y',mData.swivel_y}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.z',mData.swivel_z}; c=c+1;
            newData(c,:) = {'microscope.raw.swivel.threeD',mData.swivel_3D}; c=c+1;
            if isfield(mData,'swivel_kMT')
              newData(c,:) = {'microscope.raw.swivel.kMT',mData.swivel_kMT}; c=c+1;
            else
              newData(c,:) = {'microscope.raw.swivel.kMT',[]}; c=c+1;
            end
            
            if plane
              pData = pairedMeasurements(pCoordsInner,pCoordsOuter,1);
              % put data into string format
              newData(c,:) = {'plate.sisSep.x',pData.sisSep_x}; c=c+1;
              newData(c,:) = {'plate.sisSep.y',pData.sisSep_y}; c=c+1;
              newData(c,:) = {'plate.sisSep.z',pData.sisSep_z}; c=c+1;
              newData(c,:) = {'plate.sisSep.twoD',pData.sisSep_2D}; c=c+1;
              newData(c,:) = {'plate.sisSep.threeD',pData.sisSep_3D}; c=c+1;
              newData(c,:) = {'plate.raw.delta.x',pData.delta_x}; c=c+1;
              newData(c,:) = {'plate.raw.delta.y',pData.delta_y}; c=c+1;
              newData(c,:) = {'plate.raw.delta.z',pData.delta_z}; c=c+1;
              newData(c,:) = {'plate.raw.delta.oneD',pData.delta_1D}; c=c+1;
              newData(c,:) = {'plate.raw.delta.twoD',pData.delta_2D}; c=c+1;
              newData(c,:) = {'plate.raw.delta.threeD',pData.delta_3D}; c=c+1;
              newData(c,:) = {'plate.raw.swivel.y',pData.swivel_y}; c=c+1;
              newData(c,:) = {'plate.raw.swivel.z',pData.swivel_z}; c=c+1;
              newData(c,:) = {'plate.raw.swivel.threeD',pData.swivel_3D}; c=c+1;
              if isfield(pData,'swivel_kMT')
                newData(c,:) = {'plate.raw.swivel.kMT',pData.swivel_kMT}; c=c+1;
              else
                newData(c,:) = {'plate.raw.swivel.kMT',[]}; c=c+1;
              end
              newData(c,:) = {'plate.raw.twist.y',pData.twist_y}; c=c+1;
              newData(c,:) = {'plate.raw.twist.z',pData.twist_z}; c=c+1;
              newData(c,:) = {'plate.raw.twist.threeD',pData.twist_3D}; c=c+1;
              newData(c,:) = {'plate.sisterCentreSpeed',pData.sisCentreSpeed}; c=c+1;

              % get plate thickness measurements
%               sisCentre_x = [];
%               if length(dSinner.sisterList) == 1
%                 plateThickness = nan(nFrames,1);
%               else
%                 for jSis = 1:length(dSinner.sisterList)
%                   sisCentre_x = [sisCentre_x; nanmean([dSinner.sisterList(jSis).coords1(:,1) dSinner.sisterList(jSis).coords2(:,1)],2)];
%                 end
%                 plateThickness = nanstd(sisCentre_x,[],2);
%               end
%               % put data into string format
%               newData(c,:) = {'plate.plateThickness',plateThickness}; c=c+1;
            
            end
            
            sData = pairedMeasurements(sCoordsInner,sCoordsOuter,0);
            % put data into string format
            newData(c,:) = {'sisters.sisSep.x',sData.sisSep_x}; c=c+1;
            newData(c,:) = {'sisters.sisSep.y',sData.sisSep_y}; c=c+1;
            newData(c,:) = {'sisters.sisSep.z',sData.sisSep_z}; c=c+1;
            newData(c,:) = {'sisters.sisSep.twoD',sData.sisSep_2D}; c=c+1;
            newData(c,:) = {'sisters.sisSep.threeD',sData.sisSep_3D}; c=c+1;
            newData(c,:) = {'sisters.raw.delta.x',sData.delta_x}; c=c+1;
            newData(c,:) = {'sisters.raw.delta.y',sData.delta_y}; c=c+1;
            newData(c,:) = {'sisters.raw.delta.z',sData.delta_z}; c=c+1;
            newData(c,:) = {'sisters.raw.delta.oneD',sData.delta_1D}; c=c+1;
            newData(c,:) = {'sisters.raw.delta.twoD',sData.delta_2D}; c=c+1;
            newData(c,:) = {'sisters.raw.delta.threeD',sData.delta_3D}; c=c+1;
            newData(c,:) = {'sisters.raw.swivel.y',sData.swivel_y}; c=c+1;
            newData(c,:) = {'sisters.raw.swivel.z',sData.swivel_z}; c=c+1;
            newData(c,:) = {'sisters.raw.swivel.threeD',sData.swivel_3D}; c=c+1;

            %% Quality control
           
            % find which data satisfies the z-depth requirement
            satisfies = +(abs(mData.delta_z)<0.5*pixelSize(3));
            satisfies(satisfies==0) = NaN;
            
            % put all data into string format
            newData(c,:) = {'microscope.depthFilter.delta.x',mData.delta_x.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.y',mData.delta_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.z',mData.delta_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.oneD',mData.delta_1D.*prod(satisfies,2)}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.twoD',mData.delta_2D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.delta.threeD',mData.delta_3D.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.y',mData.swivel_y.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.z',mData.swivel_z.*satisfies}; c=c+1;
            newData(c,:) = {'microscope.depthFilter.swivel.threeD',mData.swivel_3D.*satisfies}; c=c+1;
            if isfield(mData,'swivel_kMT')
              newData(c,:) = {'microscope.depthFilter.swivel.kMT',mData.swivel_kMT.*satisfies}; c=c+1;
            else
              newData(c,:) = {'microscope.depthFilter.swivel.kMT',[]}; c=c+1;
            end
            if plane
              newData(c,:) = {'plate.depthFilter.delta.x',pData.delta_x.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.y',pData.delta_y.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.z',pData.delta_z.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.oneD',pData.delta_1D.*prod(satisfies,2)}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.twoD',pData.delta_2D.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.delta.threeD',pData.delta_3D.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.swivel.y',pData.swivel_y.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.swivel.z',pData.swivel_z.*satisfies}; c=c+1;
              newData(c,:) = {'plate.depthFilter.swivel.threeD',pData.swivel_3D.*satisfies}; c=c+1;
              if isfield(pData,'swivel_kMT')
                newData(c,:) = {'plate.depthFilter.swivel.kMT',pData.swivel_kMT.*satisfies}; c=c+1;
              else
                newData(c,:) = {'plate.depthFilter.swivel.kMT',[]}; c=c+1;
              end
              newData(c,:) = {'plate.depthFilter.twist.y',pData.twist_y.*nanmax(satisfies,[],2)}; c=c+1;
              newData(c,:) = {'plate.depthFilter.twist.z',pData.twist_z.*nanmax(satisfies,[],2)}; c=c+1;
              newData(c,:) = {'plate.depthFilter.twist.threeD',pData.twist_3D.*nanmax(satisfies,[],2)}; c=c+1;
            end
            % put all data into string format
            newData(c,:) = {'sisters.depthFilter.delta.x',sData.delta_x.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.delta.y',sData.delta_y.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.delta.z',sData.delta_z.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.delta.oneD',sData.delta_1D.*prod(satisfies,2)}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.delta.twoD',sData.delta_2D.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.delta.threeD',sData.delta_3D.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.swivel.y',sData.swivel_y.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.swivel.z',sData.swivel_z.*satisfies}; c=c+1;
            newData(c,:) = {'sisters.depthFilter.swivel.threeD',sData.swivel_3D.*satisfies}; c=c+1;
            
            % get intensities if required
            intsInnerMean = [sIinner(trackIDs(iSis,1)).intensity(:)-innerBg     sIinner(trackIDs(iSis,2)).intensity(:)-innerBg];
            intsInnerMax  = [sIinner(trackIDs(iSis,1)).intensity_max(:)-innerBg sIinner(trackIDs(iSis,2)).intensity_max(:)-innerBg];
            intsOuterMean = [sIouter(trackIDs(iSis,1)).intensity(:)-outerBg     sIouter(trackIDs(iSis,2)).intensity(:)-outerBg];
            intsOuterMax  = [sIouter(trackIDs(iSis,1)).intensity_max(:)-outerBg sIouter(trackIDs(iSis,2)).intensity_max(:)-outerBg];
            % put data into string format
            newData(c,:) = {'intensity.mean.inner',intsInnerMean}; c=c+1;
            newData(c,:) = {'intensity.mean.outer',intsOuterMean}; c=c+1;
            newData(c,:) = {'intensity.max.inner',intsInnerMax}; c=c+1;
            newData(c,:) = {'intensity.max.outer',intsOuterMax}; c=c+1;
            newData(c,:) = {'intensity.bg.inner',innerBg}; c=c+1;
            newData(c,:) = {'intensity.bg.outer',outerBg}; c=c+1;
            
            %% Directional information
            
            if nFrames > 1
                % get direction of movement
                direc = [];
                for iTrack = 1:2
                  if ~isempty(dSinner.trackList(trackIDs(iTrack)).direction)
                      if selType==1 && ~ismember(trackIDs(iTrack),theseTracks)
                        direc(:,iTrack) = nan(nFrames-1,1);
                      else
                        direc(:,iTrack) = dSinner.trackList(trackIDs(iTrack)).direction;
                      end
                  end
                end
                direc(end+1,:) = NaN;

                % get potential switch events (i.e. individual timepoints between P and AP)
                switchBuffer = 4;
                switchEvent = zeros(size(direc));
                switchDirec = [diff(direc); NaN NaN];
                for iPoint = 1:2:size(switchEvent,1)-(switchBuffer+1)
                    for iTrack = 1:2
                        tempDirec = switchDirec(iPoint:iPoint+(switchBuffer-1),iTrack);
                        idx = find(tempDirec(2:(switchBuffer-1))==0);
                        if max(abs(tempDirec))==1 && abs(sum(tempDirec))>1 && ~isempty(idx)
                            switchEvent(iPoint+idx(1):iPoint+idx(end),iTrack) = 1;
                        end
                    end
                end
                % calculate directional information
                direc_P = +(direc==1);               direc_P(direc_P==0) = NaN;
                direc_AP = +(direc==-1);             direc_AP(direc_AP==0) = NaN;
                direc_S = +(switchEvent==1);         direc_S(direc_S==0) = NaN;
                direc_N = +((direc+switchEvent)==0); direc_N(direc_N==0) = NaN;
                % put data into string format
                newData(c,:) = {'direction.P',direc_P}; c=c+1;
                newData(c,:) = {'direction.AP',direc_AP}; c=c+1;
                newData(c,:) = {'direction.S',direc_S}; c=c+1;
                newData(c,:) = {'direction.N',direc_N};
                
            end
        
            % compile new data with original
            allData = combineStrForms(allData,newData);
            
            % clear data to ensure no overlap on next loop
            clear newData
        
        end % sisters
        
      else % if paired
        
          % start counter for storing data
          c=1;
          
          % get dataStructs
          dSinner = theseMovies{iMov}.dataStruct{chanVect(1)};
          dSouter = theseMovies{iMov}.dataStruct{chanVect(2)};
          % get initCoords
          iCinner = dSinner.initCoord;
          iCouter = dSouter.initCoord;
          nSpots = size(iCinner(1).allCoord,1);
          
          % get basic metadata
          nFrames = theseMovies{iMov}.metadata.nFrames;
          pixelSize = theseMovies{iMov}.metadata.pixelSize;
          
          % get spotInts
          if isfield(dSinner,'spotInt')
            sIinner = dSinner.spotInt;
            innerBg = dSinner.cellInt.back;
            sIouter = dSouter.spotInt;
            outerBg = dSouter.cellInt.back;
          else
            sIinner = repmat(struct('intensity',nan(nFrames,1),'intensity_max',nan(nFrames,1)),1,nSpots);
            innerBg = nan(nFrames,1);
            sIouter = repmat(struct('intensity',nan(nFrames,1),'intensity_max',nan(nFrames,1)),1,nSpots);
            outerBg = nan(nFrames,1);
          end
          
          % get basic metadata
          nFrames = theseMovies{iMov}.metadata.nFrames;
          pixelSize = theseMovies{iMov}.metadata.pixelSize;
          
          % store experiment numbers
          newData(c,:) = {'exptLabel',repmat(exptLab,nSpots,1)}; c=c+1;
          
          %% Kinetochore positions
          
          % predefine various variables
          micrCoordsInner = []; micrCoordsOuter = [];
          micrCoord_x = []; micrCoord_y = []; micrCoord_z = [];
          
          for iFrame = 1:nFrames
            % get microscope coordinates of each spot
            micrCoordsInner = [micrCoordsInner; iCinner(iFrame).allCoord(:,1:3)];
            micrCoordsOuter = [micrCoordsOuter; iCouter(iFrame).allCoord(:,1:3)-mCentMeds];
            micrCoord_x = [micrCoord_x; [iCinner(iFrame).allCoord(:,1) iCouter(iFrame).allCoord(:,1)-mCentMeds(1)]];
            micrCoord_y = [micrCoord_y; [iCinner(iFrame).allCoord(:,2) iCouter(iFrame).allCoord(:,2)-mCentMeds(2)]];
            micrCoord_z = [micrCoord_z; [iCinner(iFrame).allCoord(:,3) iCouter(iFrame).allCoord(:,3)-mCentMeds(3)]];
          end
           
          % put data into string format
          newData(c,:) = {'microscope.coords.x',micrCoord_x}; c=c+1;
          newData(c,:) = {'microscope.coords.y',micrCoord_y}; c=c+1;
          newData(c,:) = {'microscope.coords.z',micrCoord_z}; c=c+1;
            
          if isfield(dSinner,'planeFit') && ~isempty(dSinner.planeFit)
            plane = 1;
            pFinner = dSinner.planeFit; pFouter = dSouter.planeFit;
          else
            plane = 0;
          end
          if plane
            
            % predefine various variables
            plateCoordsInner = []; plateCoordsOuter = [];
            plateCoord_x = []; plateCoord_y = []; plateCoord_z = [];
          
            for iFrame = 1:nFrames
              
              if ~isempty(pFinner(iFrame).planeVectors)
                %get the coordinate system of each frame with a plane
                coordSystem = pFinner(iFrame).planeVectors;
                pCentMeds = (coordSystem\mCentMeds')';
              end
              
              % get microscope coordinates of each spot
              plateCoordsInner = [plateCoordsInner; pFinner(iFrame).rotatedCoord(:,1:3)];
              plateCoordsOuter = [plateCoordsOuter; pFouter(iFrame).rotatedCoord(:,1:3)-pCentMeds];
              plateCoord_x = [plateCoord_x; [pFinner(iFrame).rotatedCoord(:,1) pFouter(iFrame).rotatedCoord(:,1)-pCentMeds(1)]];
              plateCoord_y = [plateCoord_y; [pFinner(iFrame).rotatedCoord(:,2) pFouter(iFrame).rotatedCoord(:,2)-pCentMeds(2)]];
              plateCoord_z = [plateCoord_z; [pFinner(iFrame).rotatedCoord(:,3) pFouter(iFrame).rotatedCoord(:,3)-pCentMeds(3)]];
              % put data into string format
              newData(c,:) = {'plate.coords.x',plateCoord_x}; c=c+1;
              newData(c,:) = {'plate.coords.y',plateCoord_y}; c=c+1;
              newData(c,:) = {'plate.coords.z',plateCoord_z}; c=c+1;
              
            end
          end
          
          % get intensities
          intsInnerMean=[]; intsInnerMax=[]; intsOuterMean=[]; intsOuterMax=[];
          for iSpot = 1:nSpots
            intsInnerMean = [intsInnerMean; sIinner(iSpot).intensity(:)-innerBg];
            intsInnerMax = [intsInnerMax;   sIinner(iSpot).intensity_max(:)-innerBg];
            intsOuterMean = [intsOuterMean; sIouter(iSpot).intensity(:)-outerBg];
            intsOuterMax = [intsOuterMax;   sIouter(iSpot).intensity_max(:)-outerBg];
          end
          % put data into string format
          newData(c,:) = {'intensity.mean.inner',intsInnerMean}; c=c+1;
          newData(c,:) = {'intensity.mean.outer',intsOuterMean}; c=c+1;
          newData(c,:) = {'intensity.max.inner',intsInnerMax}; c=c+1;
          newData(c,:) = {'intensity.max.outer',intsOuterMax}; c=c+1;
          
      %% Inter- and intra-kinetochore measurements

      mData = soloMeasurements(micrCoordsInner,micrCoordsOuter);
      % put data into string format
      newData(c,:) = {'microscope.raw.delta.x',mData.delta_x}; c=c+1;
      newData(c,:) = {'microscope.raw.delta.y',mData.delta_y}; c=c+1;
      newData(c,:) = {'microscope.raw.delta.z',mData.delta_z}; c=c+1;
      newData(c,:) = {'microscope.raw.delta.twoD',mData.delta_2D}; c=c+1;
      newData(c,:) = {'microscope.raw.delta.threeD',mData.delta_3D}; c=c+1;

      if plane
        pData = soloMeasurements(plateCoordsInner,plateCoordsOuter);
        % put data into string format
        newData(c,:) = {'plate.raw.delta.x',pData.delta_x}; c=c+1;
        newData(c,:) = {'plate.raw.delta.y',pData.delta_y}; c=c+1;
        newData(c,:) = {'plate.raw.delta.z',pData.delta_z}; c=c+1;
        newData(c,:) = {'plate.raw.delta.twoD',pData.delta_2D}; c=c+1;
        newData(c,:) = {'plate.raw.delta.threeD',pData.delta_3D}; c=c+1;
      end

      %% Quality control

        % find which data satisfies the z-depth requirement
        satisfies = +(abs(mData.delta_z)<0.5*pixelSize(3));
        satisfies(satisfies==0) = NaN;
        % put all data into string format
        newData(c,:) = {'microscope.depthFilter.delta.x',mData.delta_x.*satisfies}; c=c+1;
        newData(c,:) = {'microscope.depthFilter.delta.y',mData.delta_y.*satisfies}; c=c+1;
        newData(c,:) = {'microscope.depthFilter.delta.z',mData.delta_z.*satisfies}; c=c+1;
        newData(c,:) = {'microscope.depthFilter.delta.twoD',mData.delta_2D.*satisfies}; c=c+1;
        newData(c,:) = {'microscope.depthFilter.delta.threeD',mData.delta_3D.*satisfies}; c=c+1;

        if plane
          newData(c,:) = {'plate.depthFilter.delta.x',pData.delta_x.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.y',pData.delta_y.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.z',pData.delta_z.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.twoD',pData.delta_2D.*satisfies}; c=c+1;
          newData(c,:) = {'plate.depthFilter.delta.threeD',pData.delta_3D.*satisfies}; c=c+1;
        end

      % compile new data with original
      allData = combineStrForms(allData,newData);  
            
    end % paired
    
    % update progress
    prog = kitProgress(movCount/totalMovs,prog);
    
  end % movies     
end % expts

%% Save results to structure

compiledIntra = strForm2struct(allData);
compiledIntra.exptLabel = num2cell(compiledIntra.exptLabel,2);
compiledIntra.options.channels = opts.channels;
compiledIntra.options.paired = paired;
compiledIntra.options.plate = plane;

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
  fprintf('\n');
end
if ~isempty(noSis)
  fprintf('\nThe following cells contain no sisterList:\n');
  for iCell = 1:size(noSis,1)
    fprintf('    Exp %i, Mov %i\n',noSis(iCell,1),noSis(iCell,2));
  end
  fprintf('\n');
end

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


