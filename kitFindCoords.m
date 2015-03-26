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
  case 'wavelet'
    modeName = 'multiscale wavelet product';
end
kitLog(['Detecting spot candidates using ' modeName]);
prog = kitProgress(0);
for iImage = 1 : nFrames
  % get frame
  img = movie(:,:,:,iImage);
  background = fastGauss3D(img,[],filters.backgroundP,filters.border,filters.background);

  switch spotMode
    case 'histcut'
      spots = histcutSpots(img,options,dataStruct.dataProperties);
    case 'wavelet'
      spots = waveletSpots(img); % TODO allow configuration by options
  end

  %store the cands of the current image
  localMaxima(iImage).cands = spots;
  spots1D = sub2ind(size(img),round(spots(:,1)),round(spots(:,2)),round(spots(:,3)));
  localMaxima(iImage).candsAmp = img(spots1D);
  localMaxima(iImage).candsBg = background(spots1D);

  % Visualize candidates.
  if options.debug.showMmfCands ~= 0
    % If 3D image, max project.
    img = max(img,[],3);

    % Make 3 layers out of original image (normalized).
    img = img/max(img(:));
    img = repmat(img,[1 1 3]);

    % Set candidate pixels below p-value threshold red, and others blue.
    % Two passes to ensure accepted on top.
    spotsPix = round(spots);
    for i=1:length(locMax)
      if passIdx(i)==0
        img(spotsPix(i,1),spotsPix(i,2),:) = [0 0 1];
      end
    end
    for i=1:length(spotsPix)
      if passIdx(i)==1
        img(spotsPix(i,1),spotsPix(i,2),:) = [1 0 0];
      end
    end

    % Plot image.
    figure(1);
    imshow(img);
    title('Local maxima cands');
    drawnow;
    if options.debug.showMmfCands < 0
      pause;
    end
  end

  %display progress
  prog = kitProgress(iImage/nFrames,prog);
end

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

