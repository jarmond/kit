function dataStruct=kitFindCoordsMaki(job,reader,dataStruct,channel)
%KITFINDCOORDSMAKI Find coordinates using maki method
%
% SYNOPSIS: job=kitFindCoordsMaki(job,reader)
%
% INPUT job: Struct containing tracking job setup options. Following
%            fields are required:
%
%       reader: BioFormats reader.
%
%       dataStruct: Maki compatible dataStruct.
%           .initCoord: Coordinates of identified spots in each frame.
%
%       channel: One-based channel index.
%
% OUTPUT job: as input but with updated values.
%
% Copyright (c) 2007 Jonas Dorn
% Copyright (c) 2012 Ed Harry
% Copyright (c) 2015 Jonathan W. Armond

dataProperties = dataStruct.dataProperties;
metadata = job.metadata;
opts = job.options;

% Majority of following code comes from makiInitCoord.m. Modifications are
% for using BioFormats directly to read movie and stripping out cruft.

% 2D movies?
if metadata.is3D
    nDims = 3;
else
    nDims = 2;
end

filters = createFilters(nDims,dataProperties);

% for conversion to microns in the end: get pixelSize
pixelSize = metadata.pixelSize;

% get halfPatchSize to adjust centroid result. The center pixel of a 5x5x5
% array is (3,3,3), thus, we have to subtract that from centroid coordinates
halfPatchSize = dataProperties.FILTERPRM(4:6)/2+0.5;

% setup lists
nTimepoints = metadata.nFrames;
initCoordRaw = cell(nTimepoints, 1); % raw initial coords (before testing)

cropSize = job.cropSize;

% set up initCoord-Structure
tmp(1:nTimepoints,1) = struct('allCoord',[],...
    'allCoordPix',[],'correctionMu',[],'nSpots',[],...
    'initAmp',[],'amp',[],'data4MMF',[]);
dataStruct.initCoord = tmp;

% Minimum spots per frame.
if metadata.is3D
    minSpotsPerFrame = opts.minSpotsPerFrame;
else
    minSpotsPerFrame = opts.minSpotsPerFrame/4;
end

% Loop over frames.
prog = kitProgress(0);
for t=1:nTimepoints
  raw = kitReadImageStack(reader,metadata,t,channel,job.crop);

  % Remove offset.
  offset = min(raw(:));
  raw = raw - offset;

  if opts.betterBackground
    filtered = fastGauss3D(raw,[],filters.signalP,filters.border,filters.signal,1);
    % filter signal first, then run smaller filter with background to
    % save some time

    background = fastGauss3D(filtered,[],filters.backgroundP,filters.border,filters.background,1);

    % mask signal pixels, then recalculate background
    sigMask = raw>background;

    % remove some of the spurious hits. Use 2d mask for speed and
    % memory (3d strel won't work on binary image with imopen)
    sigMask = imopen(sigMask, strel('disk',1));

    % fill in binary mask holes (eharry)
    sigMask = imfill(sigMask, 'holes');

    rawMsk = raw.*(~sigMask);
    sigMask = [];

    rawMsk(rawMsk==0) = NaN;
    background = fastGauss3D(rawMsk,[],filters.signalP,1,filters.signal,1);
    background(isnan(background)) = 0;
  else
    filtered = fastGauss3D(raw,[],filters.signalP,filters.border,filters.signal);
    % filtering with sigma=15 is equal to filtering with 1 and then
    % with 14
    background = fastGauss3D(filtered,[],filters.backgroundP,filters.border,filters.background);
  end


  % amplitude is filtered image - background. This underestimates the
  % true signal. The underestimation becomes stronger if betterBackground
  % is not used.
  amplitude = filtered - background;


  % noise is local average (averaged over filter support) of squared
  % residuals of raw-filtered image
  noise = (raw-filtered).^2;
  noise = fastGauss3D(noise,[],filters.signalP,filters.border,filters.noise);

  % find local maxima
  if metadata.is3D
    locMax = loc_max3Df(filtered,[3 3 3]);
    locMaxIdx = sub2ind(cropSize,locMax(:,1),locMax(:,2),locMax(:,3));
  else
    locMaxImg = locmax2d(filtered,[3 3]);
    locMaxIdx = find(locMaxImg);
    locMax = zeros(length(locMaxIdx),3);
    [locMax(:,1),locMax(:,2)] = ind2sub(cropSize(1:2),locMaxIdx);
  end

  % read amplitude and noise at local maxima
  initCoordTmp = [locMax, amplitude(locMaxIdx), noise(locMaxIdx), ...
                  zeros(length(locMaxIdx),1)];
  initCoordTmp = sortrows(initCoordTmp,-4);
  % only take MAXSPOTS highest amplitudes. This will leave plenty of noise
  % spots, but it will avoid huge arrays
  initCoordTmp = initCoordTmp(1:min(dataProperties.MAXSPOTS,size(initCoordTmp,1)),:);

  % loop through all to get sub-pixel positions.
  raw = raw - background; %overwrite raw to save memory
  for iSpot = 1:size(initCoordTmp,1)

    % read volume around coordinate
    patch = stamp3d(raw,filters.signalP,initCoordTmp(iSpot,1:nDims),0);

    %HLE,KJ - calculate low-index edge patch adjustment if relevant
    edgeAdjustTmp = initCoordTmp(iSpot,1:3) - halfPatchSize;
    edgeAdjustTmp = abs(edgeAdjustTmp) .* edgeAdjustTmp<0;

    % subpixel coord is integer coord plus centroid (subtract
    % halfPatchSize so that center coordinate of patch is (0,0,0))
    centroid = centroid3D(patch);
    if any(centroid>2*halfPatchSize | centroid < 0)
        % Try again with exponent 2.
        centroid = centroid3D(patch,2);
        if any(centroid>2*halfPatchSize | centroid < 0)
            warning('Centroid fitting error: initCoordTmp(%d) t=%d',iSpot,t);
        end
    end
    initCoordTmp(iSpot,1:3) = ...
        initCoordTmp(iSpot,1:3) + ...
        centroid - halfPatchSize;

    %HLE,KJ - correct position due to patch truncation due to proximity
    %to a low-index edge
    initCoordTmp(iSpot,1:3) = initCoordTmp(iSpot,1:3) + edgeAdjustTmp;

    % amplitude guess is integral.
    initCoordTmp(iSpot,6) = nanmean(patch(:));

  end

  % use maxPix-amplitude to calculate cutoff - meanInt is not consistent
  % with the rest of the measures!
  initCoordTmp = [initCoordTmp,...
                  initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)),...
                  initCoordTmp(:,4)./sqrt(initCoordTmp(:,5)./max(initCoordTmp(:,4),eps))];

  initCoordRaw{t} = initCoordTmp;

  % Report progress.
  prog = kitProgress(t/nTimepoints, prog);
end

clear initCoordTmp
allCoord = cat(1,initCoordRaw{:});

%% Find cutoff

% find cutoff based on amplitude/sqrt(noise/amp), though the others are
% very similar. Allow fallback if less than 25 spots per frame are found
% (this indicates that cutFirstHistmode of splitModes failed)
cutoff = zeros(3,1);
cutoff1 = splitModes(allCoord(:,4)); % amplitude
cutoff2 = splitModes(allCoord(:,7)); % amplitude/sqrt(nse) - dark noise
cutoff3 = splitModes(allCoord(:,8)); % amplitude/sqrt(nse/amp) - poisson
if ~isempty(cutoff1)
    cutoff(1) = cutoff1;
else
    cutoff(1) = NaN;
end
if ~isempty(cutoff2)
    cutoff(2) = cutoff2;
else
    cutoff(2) = NaN;
end
if ~isempty(cutoff3)
    cutoff(3) = cutoff3;
else
    cutoff(3) = NaN;
end

% note: ask only for spots in good frames
minGood = minSpotsPerFrame*nTimepoints;
nn(3) = sum(allCoord(:,8)>cutoff(3))/minGood;
nn(2) = sum(allCoord(:,7)>cutoff(2))/minGood;
nn(1) = sum(allCoord(:,4)>cutoff(1))/minGood;
if nn(3) > 1
    cutoffIdx = 3;
    cutoffCol = 8;
elseif nn(2) > 1
    cutoffIdx = 2;
    cutoffCol = 7;
elseif nn(1) > 1
    cutoffIdx = 1;
    cutoffCol = 4;
else
    warning('less than %i spots per frame found. kitFindCoordsMaki failed',minSpotsPerFrame)
    dataStruct.failed = 1;
    return
end

dataStruct.failed = 0;
% remember the cutoff criterion used
dataStruct.statusHelp{3,3} = [cutoffIdx,cutoffCol];

% loop and store only good locMax.

for t=1:nTimepoints
    goodIdxL = initCoordRaw{t}(:,cutoffCol) > cutoff(cutoffIdx);

    % count good spots
    dataStruct.initCoord(t).nSpots = sum(goodIdxL);

    % Store pixel coords in image coordinate system (i.e. x == cols, y ==
    % rows). Uncertainty is 0.25 pix.
    dataStruct.initCoord(t).allCoordPix = ...
        [initCoordRaw{t}(goodIdxL,[2 1 3]),...
        0.25*ones(dataStruct.initCoord(t).nSpots,3)];

    % store estimated amplitude and noise
    dataStruct.initCoord(t).initAmp = initCoordRaw{t}(goodIdxL,4:5);

    % store integral amplitude
    dataStruct.initCoord(t).amp = [initCoordRaw{t}(goodIdxL,6),...
        zeros(dataStruct.initCoord(t).nSpots,1)];

    % store coords in microns and correct
    dataStruct.initCoord(t).allCoord = ...
        dataStruct.initCoord(t).allCoordPix.*...
        repmat(pixelSize,dataStruct.initCoord(t).nSpots,2);


    % Visualize final result.
    if opts.debug.showCentroidFinal ~= 0
      % If 3D image, max project.
      raw = kitReadImageStack(reader,metadata,t,channel,job.crop);
      img = max(raw,[],3);
      figure(1);
      imshow(img,[]);

      % Plot image and overlay spots.
      hold on;
      if ~isempty(initCoordRaw{t})
        plot(initCoordRaw{t}(:,2),initCoordRaw{t}(:,1),'b+');
      end
      if any(goodIdxL)
        plot(dataStruct.initCoord(t).allCoordPix(:,1), ...
             dataStruct.initCoord(t).allCoordPix(:,2),'rx');
      end

      title('Centroid fitted spots (r), cands (b)');
      hold off;
      drawnow;
      switch opts.debug.showCentroidFinal
        case -1
          pause;
        case -2
          keyboard;
      end
    end

end
