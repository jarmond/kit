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

if strcmp(spotMode,'wavelet') && options.waveletLevelAdapt
  % Compare successive frames on nSpots and covaraince of distance matrix difference
  kitLog('Determining adaptive wavelet threshold');
  tk = 1;
  tkinc = 0.1;
  tkincfac = 1.1;
  tkmax = 50;
  %f = [1 floor(nFrames/2) nFrames-1];
  f = [1 nFrames-1];
  first = 1;
  i = 1;
  while first || (tk < tkmax && nSpots > 0 && ~isnan(ld))
    options.waveletLevelThresh=tk;
    for k=1:length(f)
      [A,~,ld] = waveletSpots(movie(:,:,:,f(k)),options,dataStruct.dataProperties);
      B = waveletSpots(movie(:,:,:,f(k)+1),options,dataStruct.dataProperties);
      % Compute minimum difference between each point in A and set B.
      meanMinDiffA = 0;
      for j=1:size(A,1)
        meanMinDiffA = meanMinDiffA + min(createDistanceMatrix(A(j,:),B));
      end
      % Compute minimum difference between each point in B and set A.
      meanMinDiffB = 0;
      for j=1:size(B,1)
        meanMinDiffB = meanMinDiffB + min(createDistanceMatrix(B(j,:),A));
      end
      % Combine metrics to estimate difference in point clouds.
      frameDiff(i,k) = (meanMinDiffA + meanMinDiffB)/(size(A,1)+size(B,1));
    end

    % Increment tk.
    first = 0;
    tkvec(i) = tk;
    tk = tk + tkinc;
    tkinc = tkinc*tkincfac;
    nSpots = size(A,1);
    i = i+1;
  end
  % Pick tk which minimises frameDiff metric, with small penalty for increasing tk.
  lambda = 0.01; % penalty factor.
  frameDiff = sum(frameDiff,2) + lambda*tkvec';
  pp = pchip(tkvec,frameDiff);
  tk = fminbnd(@(x) ppval(pp,x),tkvec(1),tkvec(end));
  kitLog('Using wavelet threshold: %g',tk);
  options.waveletLevelThresh = tk;
  job.options = options;
  kitSaveJob(job); % Record used value.
end


prog = kitProgress(0);
nSpots = 0;
for iImage = 1 : nFrames
  % get frame
  img = movie(:,:,:,iImage);

  switch spotMode
    case 'histcut'
      spots = histcutSpots(img,options,dataStruct.dataProperties);
    case 'wavelet'
      spots = waveletSpots(img,options,dataStruct.dataProperties);
  end
  nSpots = nSpots + size(spots,1);

  % Round spots to nearest pixel and limit to image bounds.
  spots = bsxfun(@min,bsxfun(@max,round(spots),1),[imageSizeX,imageSizeY,imageSizeZ]);

  % Store the cands of the current image
  background = imgaussfilt3(img,filters.backgroundP(1:3),'FilterSize',filters.backgroundP(4:6));
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

    % Show candidate pixels in red.
    for i=1:size(spots,1)
      img(spots(i,1),spots(i,2),:) = [0 0 1];
    end

    % Plot image.
    figure(1);
    imshow(img);
    title(['Local maxima cands n=' num2str(size(spots,1))]);
    drawnow;
    if options.debug.showMmfCands < 0
      pause;
    end
  end

  %display progress
  prog = kitProgress(iImage/nFrames,prog);
end
kitLog('Average spots per frame: %g',nSpots/nFrames);

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

