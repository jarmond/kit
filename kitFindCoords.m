function job=kitFindCoords(job, reader, channel)
%KITFINDCOORD Find kinetochore coordinates
%
% SYNOPSIS: job=kitFindCoords(job, raw, channel)
%
% INPUT job: Struct containing tracking job setup options.
%            Requires at least the following fields:
%
%       reader: BioFormats reader.
%
%       channel: Channel to find coords in.
%
% OUTPUT job: as input but with updated values.
%
% Copyright (c) 2012 Jonathan W. Armond

% Method of fixing spot locations: centroid or Gaussian MMF, or none.
method = job.options.coordMode{channel};
% Handle old jobset versions.
if ~isfield(job.options,'spotMode')
  spotMode = 'histcut';
else
  % Method of identify first spot candidates: Histogram cut 'histcut' or
  % multiscale wavelet product 'wavelet'.
  spotMode = job.options.spotMode{channel};
end

% Set up data struct.
dataStruct = kitMakeMakiDatastruct(job, channel);
job.dataStruct{channel} = dataStruct;
options = job.options;

nFrames = job.metadata.nFrames;
is3D = job.metadata.is3D;
ndims = 2 + is3D;
filters = createFilters(ndims,dataStruct.dataProperties);

% Read image
movie = kitReadWholeMovie(reader,job.metadata,channel,job.crop,0,1);
[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);

% Initialize output structure
localMaxima = repmat(struct('cands',[]),nFrames,1);

% Find candidate spots.
switch spotMode
  case 'histcut'
    modeName = 'histogram mode cut';
  case 'adaptive'
    modeName = 'adaptive thresholding';
  otherwise
    error('Unknown spot detector: %s',spotMode);
end
kitLog(['Detecting spot candidates using ' modeName]);

switch spotMode
  case 'histcut'
    spots = cell(nFrames,1);
    for i=1:nFrames
      img = movie(:,:,:,i);
      spots{i} = histcutSpots(img,options,dataStruct.dataProperties);
    end

  case 'adaptive'
    spots = adaptiveSpots(movie,options.adaptiveLambda,options.debug.showAdaptive);
end

nSpots = zeros(nFrames,1);
for i=1:nFrames
  nSpots(i) = size(spots{i},1);

  % Round spots to nearest pixel and limit to image bounds.
  spots{i} = bsxfun(@min,bsxfun(@max,round(spots{i}),1),[imageSizeX,imageSizeY,imageSizeZ]);

  % Store the cands of the current image
  % TODO this is computed in both spot detectors, just return it.
  img = movie(:,:,:,i);
  background = imgaussfilt3(img,filters.backgroundP(1:3),'FilterSize',filters.backgroundP(4:6));
  localMaxima(i).cands = spots{i};
  spots1D = sub2ind(size(img),spots{i}(:,1),spots{i}(:,2),spots{i}(:,3));
  localMaxima(i).candsAmp = img(spots1D);
  localMaxima(i).candsBg = background(spots1D);

  % Visualize candidates.
  if options.debug.showMmfCands ~= 0
    showSpots(img,spots{i});
    title(['Local maxima cands n=' num2str(size(spots{i},1))]);
    drawnow;
    switch options.debug.showMmfCands
      case -1
        pause;
      case -2
        keyboard;
    end
  end
end
kitLog('Average spots per frame: %.1f +/- %.1f',mean(nSpots),std(nSpots));

% Refine spot candidates.
switch method
  case 'centroid'
    job = kitCentroid(job,movie,localMaxima,channel);
  case 'gaussian'
    job = kitMixtureModel(job,movie,localMaxima,channel);
  case 'norefine'
    % No refinement. Copy localMaxima to initCoords.
    initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0, ...
                                  'amp',[],'bg',[]);
    initCoord(1).localMaxima = localMaxima;
    for i=1:nFrames
      initCoord(i).nSpots = size(localMaxima(i).cands,1);
      initCoord(i).allCoordPix = [localMaxima(i).cands(:,[2 1 3]) ...
                          0.25*ones(initCoord(i).nSpots,3)];
      initCoord(i).allCoord = bsxfun(@times, initCoord(i).allCoordPix,...
        repmat(job.metadata.pixelSize,[1 2]));
      initCoord(i).amp = [localMaxima(i).candsAmp zeros(initCoord(i).nSpots,1)];
    end
    % Store data.
    job.dataStruct{channel}.initCoord = initCoord;
    job.dataStruct{channel}.failed = 0;
  otherwise
    error(['Unknown coordinate finding mode: ' job.coordMode]);
end

