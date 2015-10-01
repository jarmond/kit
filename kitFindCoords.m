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
  tk = 0.1;
  tkinc = 0.1;
  tkincfac = 1.1;
  tkmax = 50;
  %f = [1 floor(nFrames/4) floor(nFrames/2) floor(3*nFrames/4) nFrames-1];
  f = [1  floor(nFrames/2) nFrames-1];
  %f = [1 nFrames-1];
  first = 1;
  i = 1;
  while first || (tk < tkmax && nSpots > 0 && ~isnan(ld))
    options.waveletLevelThresh=tk;
    for k=1:length(f)
      [A,~,ld] = waveletSpots(movie(:,:,:,f(k)),options);
      B = waveletSpots(movie(:,:,:,f(k)+1),options);
      frameDiff(i,k) = meanMinDiff(A,B);
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
  lambda = 0.1; % penalty factor.
  frameDiff = mean(frameDiff,2) + lambda*tkvec';
  [~,minIdx] = min(frameDiff);
  pp = pchip(tkvec,frameDiff);
  tk = fminbnd(@(x) ppval(pp,x),tkvec(max(minIdx-1,1)),tkvec(min(minIdx+1,length(tkvec))));
  kitLog('Using wavelet threshold: %g',tk);

  if options.debug.showWaveletAdapt
    figure;
    x = linspace(tkvec(1),tkvec(end));
    plot(x,ppval(pp,x));
    hold on
    plot([tk tk],ylim,'r--');
    hold off
    ylabel('Frame-frame point cloud difference');
    xlabel('Wavelet threshold');
    drawnow;
  end

  options.waveletLevelThresh = tk;
  job.options = options;
  kitSaveJob(job); % Record used value.
end


prog = kitProgress(0);
nSpots = zeros(nFrames,1);
for iImage = 1 : nFrames
  % get frame
  img = movie(:,:,:,iImage);

  switch spotMode
    case 'histcut'
      spots = histcutSpots(img,options,dataStruct.dataProperties);
    case 'wavelet'
      spots = waveletSpots(img,options);
  end
  nSpots(iImage) = size(spots,1);

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
    showSpots(img,spots);
    title(['Local maxima cands n=' num2str(size(spots,1))]);
    drawnow;
    switch options.debug.showMmfCands
      case -1
        pause;
      case -2
        keyboard;
    end
  end

  %display progress
  prog = kitProgress(iImage/nFrames,prog);
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

