function jobs = dublManualPairSisters(jobs,varargin)
% KITMANUALPAIRSISTERS Manually select pairs of kinetochores.
%
%    JOB = DUBLMANUALPAIRSISTERS(JOB,...) Images of each kinetochore found
%    in a JOB are shown and their centre plotted (in white). Centres of
%    nearby kinetochores (potential sisters, in green) are also shown, so
%    that the user can double-click the spot they believe to be the
%    original spot's sister.
%    Kinetochores with no clear sister can be omitted from the sister list
%    by pressing any key.
%    Nearby kinetochores already allocated (red), and kinetochores
%    previously decided to have no clear pair (yellow), are also plotted.
%
%    Options, defaults in {}:-
%
%    contrast: {{[0.1 1],[0.1 1],[0.1 1]}} or similar structure. A 1x3 cell
%                   of pairs of numbers for use in image contrast.
%
%    coordChans: {[1 2]} or any subset of [1,2,3]. The channels for which
%                   sisterList structures will be produced.
%
%    imageChans: {[1 2]} or any subset of [1,2,3]. The channels which will
%                   be plotted when inspecting spots for pairing.
%
%    maxSisSep: {2} or distance in µm. Maximum distance between the
%                   original spot and potential sisters.
%
%    metaphasePlate: {0} or 1. Whether or not to plot a metaphase plate.
%                   N.B. This doesn't work yet. Future plan.
%
%    mode: {'dual'}, 'zoom' or 'full'. Whether to show images as:
%                   - 'full'
%                     full images.
%                   - 'zoom'
%                     zooming in on kinetochores to an area defined by the
%                     'maxSisSep' option.
%                   - 'dual'
%                     a combination of the two.
%
%    plotChan: {1}, 2 or 3. The channel which will be used for plotting
%                   spot centres when inspecting spots for pairing.
%                   N.B. This only works for one channel so far, but will
%                   be increased in future.
%
%    redo: {0} or 1. Whether to completely redo sister pairing, i.e. don't
%                   plot spots as being previously allocated etc.
%                   N.B. This doesn't work yet. Future plan.
%
%    subset: {[]} or list of movie numbers. Sub-list of movies to process.
%
%    verbose: {0} or 1. Whether or not to print progress to the command
%                   line.
%
%    zProjRange: {2} or distance in pixels. Total distance about the
%                   original spot over which to project in the
%                   z-coordinate.
%
%
%    Future directions:-
%
%    - Incorporate plotting of metaphase plate
%    - Allow processing of movies as well as single z-stacks
%       - Looking for previously discovered sisters before manual selection
%       - Take three points per track to ensure faithful pairing
%
% Copyright (c) 2016 C. A. Smith

% default options
opts.chanOrder = [2 1 3]; % order of the channels in the movie
opts.contrast = {[0.1 1] [0.1 1] [0.1 1]};
opts.coordChans = [1 2]; % channels to process
opts.imageChans = [1 2]; % imaging more channels may provide more information
opts.maxSisSep = 2; % maximum distance in µm between potential sisters being plotted
opts.mode = 'dual';
opts.plotChan = 1; % channel coordinates to be plotted for allocation
opts.redo = 0;
opts.subset = [];
opts.verbose = 0;
opts.zProjRange = 2; % number of pixels over which to project in z

% process user-defined options
opts = processOptions(opts,varargin{:});

% turn off all warnings
warning('off','all')

%% FRONT END
% if no subset given, use full length of jobs
if isempty(opts.subset)
    opts.subset = 1:length(jobs);
end

% start counter
counter = 1;
nMovies = length(opts.subset);

% loop over movies to run kitManualPairSisters
for iMov = opts.subset
    
    % output log information
    kitLog('Pairing kinetochores in movie %i (%i of %i):',iMov,counter,nMovies);
    counter = counter+1;
    
    % perform the sister pairing
    [jobs{iMov},userStatus] = pairSisters(jobs{iMov},opts);
    
    
    if strcmp(userStatus,'userPaused')
        [~,lastMovie] = find(opts.subset == iMov);
        if lastMovie == 1
            lastMovie = 0;
            break
        else
            lastMovie = opts.subset(lastMovie-1);
            break
        end
    end
end

% check user input
switch userStatus
    case 'userPaused'
        if lastMovie == 0
            kitLog('User paused sister pairing before any movies saved.')
            return
        else
            kitLog('User paused sister pairing. Sister pairing saved up to movie %i.',lastMovie)
        end
    case 'completed'
        kitLog('Sister pairing completed for all movies.')
end

end

function [job,userStatus] = pairSisters(job,opts)
%% GET REQUIRED IMAGE AND METADATA
[md, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie),job.metadata);
movieIdx = job.index;

% get crop information, if any
crop = job.ROI(movieIdx).crop;
cropSize = job.ROI(movieIdx).cropSize;
if isempty(crop)
    cropSize = md.frameSize;
end
% specify RGB channel order [R G B]
chanOrder = opts.chanOrder;

% produce image file
fullImg = zeros([cropSize(1:2),3]);
for iChan = opts.imageChans
    img(:,:,:,iChan) = kitReadImageStack(reader, md, 1, iChan, crop, 0);
    fullImg(:,:,chanOrder(iChan)) = max(img(:,:,:,iChan),[],3); % full z-project
    irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
    fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
end

% produce figure environment
figure(1)
clf

plotTitle = sprintf('Check cell orientation\nPress: y to continue, n to skip, q to quit.');
if length(opts.imageChans) == 1
    imshow(fullImg(:,:,chanOrder(opts.imageChans)));
    title(plotTitle,'FontSize',14)
else
    imshow(fullImg)
    title(plotTitle,'FontSize',14)
end
[~,~,buttonPress] = ginput(1);

% calculate size of cropped region in pixels based on maxSisSep
pixelSize = job.metadata.pixelSize(1:3);
cropRange = 1.25*repmat(opts.maxSisSep,1,3)./pixelSize;
cropRange = round(cropRange);
chrShift = job.options.chrShift.result;
% find chrShift rounded to nearest pixel
pixChrShift = cellfun(@times,chrShift,repmat({[pixelSize pixelSize]},3),'UniformOutput',0);


%% GET ALL COORDINATES

% get coordinates in both µm and pixels
coords = job.dataStruct{opts.plotChan}.initCoord(1).allCoord(:,[2 1 3]);
coordsPix = job.dataStruct{opts.plotChan}.initCoord(1).allCoordPix(:,[2 1 3]);

nCoords = size(coords,1);


%% PLOT IMAGE, REQUEST INPUT

% pre-allocate index vectors
pairedIdx = [];
unpairedIdx = [];
if buttonPress == 110 %110 corresponds to user pressing n key
    unallocatedIdx = [];
elseif buttonPress == 121 %121 corresponds to user pressing y key
    unallocatedIdx = 1:nCoords;
    if opts.verbose
        fprintf('Cell accepted. Continuing with pairing of sisters.\n');
    end
elseif buttonPress == 113 %113 corresponds to user pressing q key
    userStatus = 'userPaused';
    return
else
    unallocatedIdx = 1:nCoords;
    fprintf('Another key other than y, n or q pressed. Continuing with kinetochore pairing.\n');
end

doublePairingIdx = [];
sisterIdxArray = [];

plotTitle = sprintf('Locate the white cross'' sister.\nClick on an: unallocated (g), ignored (y) or pre-paired (r) cross.\nPress backspace if none applicable.');

% preconfigure progress information
prog = kitProgress(0);

while ~isempty(unallocatedIdx)
    
    % give progress information
    nRemaining = size(unallocatedIdx,2);
    prog = kitProgress((nCoords-nRemaining)/nCoords,prog);
    
    % take the earliest spotIdx
    iCoord = unallocatedIdx(1);

    % get nnDist information
    unallNNidx    = getNNdistIdx(coords,iCoord,unallocatedIdx,opts.maxSisSep);
    pairedNNidx   = getNNdistIdx(coords,iCoord,pairedIdx,opts.maxSisSep);
    unpairedNNidx = getNNdistIdx(coords,iCoord,unpairedIdx,opts.maxSisSep);
    
    % if spot has no nearest neighbour, remove and continue
    if isempty([unallNNidx pairedNNidx unpairedNNidx])
        unpairedIdx = [unpairedIdx iCoord];
        unallocatedIdx(1) = [];
        continue
    end
    
    % get origin pixel-coordinates
    iCoordsPix = coordsPix(iCoord,:);
    centreCoords = round(iCoordsPix);
    % get pixel-coordinates for each nnDistIdx
    unallCoordsPix = coordsPix(unallNNidx,:);
    pairedCoordsPix = coordsPix(pairedNNidx,:);
    unpairedCoordsPix = coordsPix(unpairedNNidx,:);
    % compile all coords for later
    frameCoordsPix = [iCoordsPix; unallCoordsPix; pairedCoordsPix; unpairedCoordsPix];
    
    % IMAGE PROCESSING
    % predesignate cropImg structure
    coordRange = [max(1,centreCoords(1)-cropRange(1)) min(centreCoords(1)+cropRange(1),cropSize(1));...
                  max(1,centreCoords(2)-cropRange(2)) min(centreCoords(2)+cropRange(2),cropSize(2));...
                  max(1,centreCoords(3)-opts.zProjRange) min(centreCoords(3)+opts.zProjRange,cropSize(3))];
    cropImg = zeros(coordRange(1,2)-coordRange(1,1)+1,...
                    coordRange(2,2)-coordRange(2,1)+1,...
                    3);
    
	% plot image in figure 1
    figure(1)
    clf
                
    switch opts.mode
        
        case 'zoom'
    
            % get zoomed image centred at origin coordinate
            tempImg = img(coordRange(1,1):coordRange(1,2), coordRange(2,1):coordRange(2,2), :, :);
            for iChan = opts.imageChans
                if iChan ~= opts.plotChan
                    iCoordRange = coordRange - pixChrShift{opts.plotChan,iChan}(1:2);
                else
                    iCoordRange = coordRange;
                end
                cropImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, iCoordRange(3,1):iCoordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(cropImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                cropImg(:,:,chanOrder(iChan)) = imadjust(cropImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end

            if length(opts.imageChans) == 1
                imshow(cropImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(cropImg)
            end
            title(plotTitle,'FontSize',14)

            % PLOTTING ZOOMED COORDINATES
            hold on
            % origin coordinates in white
            scatter(iCoordsPix(2) - coordRange(2,1)+1, iCoordsPix(1) - coordRange(1,1)+1,'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,2) - coordRange(2,1)+1, unallCoordsPix(:,1) - coordRange(1,1)+1,'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,2) - coordRange(2,1)+1, unpairedCoordsPix(:,1) - coordRange(1,1)+1,'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,2) - coordRange(2,1)+1, pairedCoordsPix(:,1) - coordRange(1,1)+1,'xr','sizeData',200,'LineWidth',1.25);
        
        case 'full'
        
            % get image
            tempImg = img;
            for iChan = opts.imageChans
                fullImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, coordRange(3,1):coordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end

            if length(opts.imageChans) == 1
                imshow(fullImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(fullImg)
            end
            title(plotTitle,'FontSize',14)
        
            % PLOTTING COORDINATES
            hold on
            % origin coordinates in white
            scatter(iCoordsPix(2),iCoordsPix(1),'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,2),unallCoordsPix(:,1),'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,2),unpairedCoordsPix(:,1),'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,2),pairedCoordsPix(:,1),'xr','sizeData',200,'LineWidth',1.25);
        
        case 'dual'
            
            subplot(2,4,[3,4,7,8])
            % get full image
            tempImg = img;
            for iChan = opts.imageChans
                fullImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, coordRange(3,1):coordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(fullImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                fullImg(:,:,chanOrder(iChan)) = imadjust(fullImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end
            
            if length(opts.imageChans) == 1
                imshow(fullImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(fullImg)
            end
            
            % PLOTTING COORDINATES
            hold on
            % origin coordinates in white
            scatter(iCoordsPix(2),iCoordsPix(1),'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,2),unallCoordsPix(:,1),'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,2),unpairedCoordsPix(:,1),'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,2),pairedCoordsPix(:,1),'xr','sizeData',200,'LineWidth',1.25);
            
            subplot(2,4,[1,2,5,6])
            % get zoomed image
            tempImg = img(coordRange(1,1):coordRange(1,2), coordRange(2,1):coordRange(2,2), :, :);
            for iChan = opts.imageChans
                cropImg(:,:,chanOrder(iChan)) = max(tempImg(:,:, coordRange(3,1):coordRange(3,2), iChan),[],3);
                irange(iChan,:) = stretchlim(cropImg(:,:,chanOrder(iChan)),opts.contrast{iChan});
                cropImg(:,:,chanOrder(iChan)) = imadjust(cropImg(:,:,chanOrder(iChan)),irange(iChan,:), []);
            end
            
            if length(opts.imageChans) == 1
                imshow(cropImg(:,:,chanOrder(opts.imageChans)));
            else
                imshow(cropImg)
            end
            title(plotTitle,'FontSize',14)
            
            % PLOTTING ZOOMED COORDINATES
            hold on
            % origin coordinates in white
            scatter(iCoordsPix(2) - coordRange(2,1)+1, iCoordsPix(1) - coordRange(1,1)+1,'xw','sizeData',200,'LineWidth',1.25);
            % unallocated in green
            scatter(unallCoordsPix(:,2) - coordRange(2,1)+1, unallCoordsPix(:,1) - coordRange(1,1)+1,'xg','sizeData',200,'LineWidth',1.25);
            % unpaired in yellow
            scatter(unpairedCoordsPix(:,2) - coordRange(2,1)+1, unpairedCoordsPix(:,1) - coordRange(1,1)+1,'xy','sizeData',200,'LineWidth',1.25);
            % paired in red
            scatter(pairedCoordsPix(:,2) - coordRange(2,1)+1, pairedCoordsPix(:,1) - coordRange(1,1)+1,'xr','sizeData',200,'LineWidth',1);
            
    end
    
    
    % USER INPUT AND SISTER-SISTER PAIRING
    
    % uses a crosshair to allow user to provide a pair of coordinates
    [userY,userX,key] = ginput(1);
    
    if ismember(key,1:3)
        
        % correct user-provided information if zoomed
        if any(strcmp(opts.mode,{'dual','zoom'}))
            userX = userX + coordRange(1,1)-1; userY = userY + coordRange(2,1)-1;
        end
        
        % find the user-selected cross
        userNNdist = createDistanceMatrix([userX userY],frameCoordsPix(:,1:2));
        [~,userIdx] = min(userNNdist);
        [userIdx,~] = find((coordsPix - repmat(frameCoordsPix(userIdx,:),size(coordsPix,1),1)) == 0);
        userIdx = userIdx(1);
        
        % check user hasn't given the original cross
        if userIdx == iCoord
            fprintf('White cross cannot be selected. Please try again.\n')
            continue
        end
            
        % check from which group the new spot is, and process accordingly
        if any(userIdx == unallocatedIdx)
            
            % remove new coordinates from unallocated lists, add to paired
            unallocatedIdx = setdiff(unallocatedIdx,[iCoord userIdx]);
            pairedIdx = [pairedIdx iCoord userIdx];
            
        elseif any(userIdx == unpairedIdx)
            
            % remove new coordinates from unpaired and unallocated lists,
            % add to paired
            unallocatedIdx = setdiff(unallocatedIdx,iCoord);
            unpairedIdx = setdiff(unpairedIdx,userIdx);
            pairedIdx = [pairedIdx iCoord userIdx];
            
        elseif any(userIdx == pairedIdx)
            
            % remove new coordinate from unallocated, add to paired
            unallocatedIdx = setdiff(unallocatedIdx,iCoord);
            pairedIdx = [pairedIdx iCoord];
            
            % locate the original pair for the new paired index
            [iRow,iCol] = find(sisterIdxArray == userIdx);
            oldPairIdx = sisterIdxArray(iRow,setdiff([1 2],iCol));
            
            % add originally-paired index to unallocated, remove old pair
            % from sister array
            unallocatedIdx = [unallocatedIdx oldPairIdx];
            pairedIdx = setdiff(pairedIdx,oldPairIdx);
            sisterIdxArray(iRow,:) = [];
            
            % save that the common sister has been chosen twice - this may
            % not be necessary, but keep for the time being
            doublePairingIdx = [doublePairingIdx; oldPairIdx userIdx];
            
        end
        % save sister pair indices
        if opts.verbose
            fprintf('Sisters paired: %i and %i.\n',iCoord,userIdx)
        end
        sisterIdxArray = [sisterIdxArray; iCoord userIdx];
    
    else % coordinate not paired with anything
        unpairedIdx = [unpairedIdx iCoord];
        unallocatedIdx = setdiff(unallocatedIdx,iCoord);
        if opts.verbose
            fprintf('No sister paired with %i.\n',iCoord)
        end
    end
    
end

prog = kitProgress(1,prog);

userStatus = 'completed';

%% STORE INFORMATION

% set up empty sisterLists
emptySisterList = struct('trackPairs',[],...
                        'coords1',[],...
                        'coords2',[],...
                        'sisterVectors',[],...
                        'distances',[]);
emptyTracks = struct('tracksFeatIndxCG',0,...
               'tracksCoordAmpCG',[],...
               'seqOfEvents', [1 1 1; 1 2 1],...
               'coordAmp4Tracking',[]);
                    
for iChan = opts.coordChans
    
    % construct new sisterList
    sisterList = emptySisterList;
    sisterList(1).trackPairs(:,1) = 1:size(sisterIdxArray,1);
    sisterList(1).trackPairs(:,2) = sisterList(1).trackPairs(:,1)+size(sisterIdxArray,1);
    
    % get microscope coordinates and standard deviations in both µm
    coords = job.dataStruct{iChan}.initCoord(1).allCoord;
    % get plane coordinates if applicable
    if isfield(job.dataStruct{iChan},'planeFit') && ~isempty(job.dataStruct{iChan}.planeFit(1).plane)
      rotCoords = job.dataStruct{iChan}.planeFit(1).rotatedCoord;
    else
      rotCoords = coords;
    end
    
    for iSis = 1:size(sisterIdxArray,1)
        
        sisterList(iSis).coords1 = rotCoords(sisterIdxArray(iSis,1),:);
        sisterList(iSis).coords2 = rotCoords(sisterIdxArray(iSis,2),:);
        sisterList(iSis).distances = sqrt(sum((sisterList(iSis).coords1 - sisterList(iSis).coords2).^2,2));
        
    end
    
    % ensure tracks are re-allocated for new cell
    tracks = emptyTracks;
    % construct new tracks
    for iTrack = 1:length(sisterIdxArray(:))
        
        tracks(iTrack) = emptyTracks;
        tracks(iTrack).tracksFeatIndxCG = sisterIdxArray(iTrack);
        tracks(iTrack).tracksCoordAmpCG = coords(sisterIdxArray(iTrack),:);
        tracks(iTrack).coordAmp4Tracking = rotCoords(sisterIdxArray(iTrack),:);
        
    end
                 
    job.dataStruct{iChan}.tracks = tracks;
    job.dataStruct{iChan}.sisterList = sisterList;
    job = kitExtractTracks(job,iChan);
    job = kitSaveJob(job);
    
end
end
   
function distIdx = getNNdistIdx(coords,origIdx,otherIdx,maxSisSep)

    nnDists = createDistanceMatrix(coords(origIdx,:),coords(otherIdx,:));
    nnCoords = (nnDists < maxSisSep & nnDists > 0);
    distIdx = otherIdx(nnCoords);
    
end