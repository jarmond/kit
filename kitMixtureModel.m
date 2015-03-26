function job = kitMixtureModel(job,movie,localMaxima,channel)
% Refine spots in 3D using Gaussian mixture model fits
%
% Created by: Jacques Boisvert and K. Jaqaman
% Modified by: J. W. Armond
% Copyright (c) 2014 Jonathan W. Armond

%% Input + initialization

dataStruct = job.dataStruct{channel};
options = job.options;

% get number of frames
nFrames = job.metadata.nFrames;
is3D = job.metadata.is3D;
ndims = 2 + is3D;
%get initial guess of PSF sigma
filterPrm = dataStruct.dataProperties.FILTERPRM;
%initialize some variables
emptyFrames = [];

[imageSizeX,imageSizeY,imageSizeZ,~] = size(movie);

% Go over all frames and register empty frames.
for iImage = 1 : nFrames
  %if there are no cands, register that this is an empty frame
  if isempty(localMaxima(iImage).cands)
    emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
  end
end

%make a list of images that have local maxima
goodImages = setxor(1:nFrames,emptyFrames,'legacy');

% get psf sigma from filterPrm
if is3D
  psfSigma = [filterPrm(1) filterPrm(3)];
else
  psfSigma = filterPrm(1);
end

% Minimum distance in multiples of initial PSF sigma between features
% for estimating.
featSep = [1 1]*options.clusterSeparation;

%% PSF sigma estimation
if options.maxPsfSigmaIters

  %specify which parameters to fit for sigma estimation
  if is3D
    fitParameters = [{'X1'} {'X2'} {'X3'} {'A'} {'Sxy'} {'S3'} {'B'}]; % 3d  approx
  else
    fitParameters = [{'X1'} {'X2'} {'A'} {'Sxy'} {'B'}];
  end

  %store original input sigma
  psfSigmaIn = psfSigma;

  %give a dummy value for psfSigma0 and acceptCalc to start while loop
  psfSigma0 = [0 0];
  acceptCalc = 1;

  %initialize variable counting number of iterations
  numIter = 0;

  %iterate as long as estimated sigma is larger than initial sigma
  while numIter <= options.maxPsfSigmaIters && acceptCalc && any((abs(psfSigma - psfSigma0) ./ psfSigma0) > [0.05 0.05]) % fit both sigma-xy and sigma-z
                                                                                                                         %add one to number of iterations
    numIter = numIter + 1;

    %save input PSF sigma in new variable and empty psfSigma for estimation
    psfSigma0 = psfSigma;
    psfSigma = [];

    %calculate some numbers that get repeated many times
    psfSigma5 = ceil(5*psfSigma0);

    %initialize progress display
    switch numIter
      case 1
        currEst = sprintf('[%1.3f %1.3f]',psfSigma0(1),psfSigma0(2));
        kitLog('Estimating PSF sigma = %s',currEst);
      otherwise
        kitLog('Repeating PSF sigma estimation');
    end
    prog = kitProgress(0);

    %go over images
    for iImage = goodImages(1:min(50,length(goodImages)))

      %read image
      image = movie(:,:,:,iImage);

      %get feature positions
      featPos = localMaxima(iImage).cands;
      featBg = localMaxima(iImage).candsBg;
      featAmp = localMaxima(iImage).candsAmp;

      %retain only features that are more than 5*psfSigma0 away from boundaries
      if is3D
        feat2use = find(featPos(:,1) > psfSigma5(1) & ...
                        featPos(:,1) < imageSizeX - psfSigma5(1) & ...
                        featPos(:,2) > psfSigma5(1) & ...
                        featPos(:,2) < imageSizeY - psfSigma5(1) & ...
                        featPos(:,3) > psfSigma5(2) & ...
                        featPos(:,3) < imageSizeZ - psfSigma5(2));
      else
        feat2use = featPos(:,1) > psfSigma5(1)...
            & featPos(:,1) < imageSizeX - psfSigma5(1)...
            & featPos(:,2) > psfSigma5(1)...
            & featPos(:,2) < imageSizeY - psfSigma5(1);
      end
      featPos = featPos(feat2use,:);

      %if there is more than one feature ...
      if length(feat2use) > 1

        % XY distance matrix
        distXY = createDistanceMatrix(featPos(:,1:2),featPos(:,1:2));
        % find xy and z distance above min distances
        distXY = distXY > ceil(featSep(1)*psfSigma0(1));
        % make the diagonals for each matrix true
        distXY(logical(eye(size(distXY)))) = true;

        if is3D
          % Z distance matrix
          distZ = abs(createDistanceMatrix(featPos(:,3),featPos(:,3)));
          distZ = distZ > ceil(featSep(2)*psfSigma0(2));
          % make the diagonals for each matrix true
          distZ(logical(eye(size(distZ)))) = true;
        end

        % now find rows that have all true in XY OR Z
        if is3D
          feat2use = all(distXY | distZ, 2);
        else
          feat2use = all(distXY, 2);
        end

        featPos = featPos(feat2use,:);
        featBg = featBg(feat2use);
        featAmp = featAmp(feat2use);
      end

      %go over the selected features and estimate psfSigma
      numFeats = size(featPos,1);

      parameters = zeros(numFeats,length(fitParameters));
      if numFeats >= 1

        for iFeat = 1 : numFeats

          %crop image around selected feature
          lowerBoundXY = featPos(iFeat,1:2) - psfSigma5(1);
          upperBoundXY = featPos(iFeat,1:2) + psfSigma5(1);
          if is3D
            lowerBoundZ = featPos(iFeat,3) - psfSigma5(2);
            upperBoundZ = featPos(iFeat,3) + psfSigma5(2);
            imageCropped = image(lowerBoundXY(1):upperBoundXY(1),...
                                    lowerBoundXY(2):upperBoundXY(2),...
                                    lowerBoundZ:upperBoundZ);
          else
            imageCropped = image(lowerBoundXY(1):upperBoundXY(1),...
                                    lowerBoundXY(2):upperBoundXY(2));
          end

          %estimate sigma if image region contains no NaNs
          %NaNs appear due to cropping
          if all(~isnan(imageCropped(:)))

            % initial guess of parameters
            if is3D
              initGuess = [psfSigma5(1)+1 psfSigma5(1)+1 psfSigma5(2)+1 ...
                           featAmp(iFeat) psfSigma0(1) psfSigma0(2) featBg(iFeat)];
            else
              initGuess = [psfSigma5(1)+1 psfSigma5(1)+1 ...
                           featAmp(iFeat) psfSigma0(1) featBg(iFeat)];
            end

            %fit image and estimate sigma of Gaussian
            parameters(iFeat,:) = GaussFitND(imageCropped,[],fitParameters,initGuess);

          else %otherwise assign NaN
            parameters(iFeat,:) = NaN;
          end

        end

        %add to array of sigmas
        if is3D
          psfSigma = [psfSigma; parameters(:,5:6)]; %#ok<AGROW>
        else
          psfSigma = [psfSigma; parameters(:,4)]; %#ok<AGROW>
        end

      end %(if numFeats >= 1)

      %display progress
      prog = kitProgress(iImage/min(50,length(goodImages)),prog);

    end %(for iImage = images2use)

    %estimate psfSigma as the robust mean of all the sigmas from the fits
    %get rid of NaNs from cropped regions
    if ~isempty(psfSigma)
      if is3D
        psfSigma = psfSigma((~isnan(psfSigma(:,1)) & ~isnan(psfSigma(:,2))),:);
      else
        psfSigma = psfSigma(~isnan(psfSigma));
      end
    end
    numCalcs = size(psfSigma,1);
    if numCalcs > 0
      if is3D
        [psfSigmaXY_temp,~,inlierIndxXY_temp] = robustMean(psfSigma(:,1));
        [psfSigmaZ_temp,~,inlierIndxZ_temp] = robustMean(psfSigma(:,2));

        psfSigma = [psfSigmaXY_temp psfSigmaZ_temp];

        acceptCalc = (numCalcs >= 100 && length(inlierIndxXY_temp) >= 0.5*numCalcs && length(inlierIndxZ_temp) >= 0.5*numCalcs) || ...
            (numCalcs >= 50 && length(inlierIndxXY_temp) >= 0.7*numCalcs && length(inlierIndxZ_temp) >= 0.7*numCalcs) || ...
            (numCalcs >= 10 && length(inlierIndxXY_temp) >= 0.9*numCalcs ...
             && length(inlierIndxZ_temp) >= 0.9*numCalcs);
      else
        [psfSigma,~,inlierIndx] = robustMean(psfSigma);
        %accept new sigma if there are enough observations and inliers
        acceptCalc = (numCalcs >= 100 && length(inlierIndx) >= 0.7*numCalcs) || ...
            (numCalcs >= 50 && length(inlierIndx) >= 0.9*numCalcs) || ...
            (numCalcs >= 10 && length(inlierIndx) == numCalcs);
      end
    else

      acceptCalc = 0;

    end

    %show new sigma if estimation is accepted
    if acceptCalc
      if is3D
        currEst = sprintf('[%1.3f %1.3f]',psfSigma(1),psfSigma(2));
      else
        currEst = sprintf('%1.3f',psfSigma(1));
      end
      kitLog('Estimated PSF sigma = %s from %d features',currEst,numCalcs);
    else %otherwise alert user that input sigma was retained
      psfSigma = psfSigmaIn;
      kitLog('Not enough observations to change PSF sigma, using input PSF sigma');
    end

  end %(while numIter <= numSigmaIter && acceptCalc && ((psfSigma-psfSigma0)/psfSigma0 > 0.05))

  %if maximum number of iterations has been performed but sigma value is not converging
  if numIter == options.maxPsfSigmaIters+1 && acceptCalc && any((abs(psfSigma - psfSigma0) ./ psfSigma0) > [0.05 0.05])
    psfSigma = psfSigmaIn;
    kitLog('Estimation terminated (no convergence), using input PSF sigma');
  end

end %(if numSigmaIter)


%% Mixture-model fitting

% save the pixelSize
pixelSize = job.metadata.pixelSize;

%initialize initCoord
initCoord(1:nFrames) = struct('allCoord',[],'allCoordPix',[],'nSpots',0,'amp',[],'bg',[]);

kitLog('Refining spots using mixture-model fitting');
prog = kitProgress(0);
%go over all non-empty images ...
for iImage = goodImages

    % Get frame.
    imageRaw = movie(:,:,:,iImage);

    % Get candidate maxima.
    cands = localMaxima(iImage).cands;

    % Fit with mixture-models.
    [coordList,ampList,bgList,rejects] = mixtureModelFit(cands,imageRaw,psfSigma,options);
    nSpots = size(coordList,1);
    if ~is3D
      coordList = [coordList(:,1:2) zeros(nSpots,1) coordList(:,3:4) zeros(nSpots,1)];
    end

    % Visualize final result.
    if options.debug.showMmfFinal ~= 0
      % If 3D image, max project.
      img = max(imageRaw,[],3);
      figure(1);
      imshow(img,[]);

      % Plot image and overlay spots.
      hold on;
      plot(cands(:,2),cands(:,1),'b+');
      if ~isempty(rejects.amp)
        plot(rejects.amp(:,1),rejects.amp(:,2),'gx');
      end
      if ~isempty(rejects.dist)
        plot(rejects.dist(:,1),rejects.dist(:,2),'yx');
      end
      if ~isempty(coordList)
        plot(coordList(:,1),coordList(:,2),'rx');
      end

      title('MMF fitted spots (r), cands (b), amp rej (g), dist rej (y)');
      hold off;
      drawnow;
      switch options.debug.showMmfFinal
        case -1
          pause;
        case -2
          keyboard;
      end
    end

    if options.debug.showMmfPvals ~= 0
      figure(2);
      subplot(2,1,1);
      if ~isempty(rejects.amp)
        histogram(rejects.amp(:,ndims+1));
      end
      title('Amplitude reject p-vals');
      subplot(2,1,2);
      if ~isempty(rejects.dist)
        histogram(rejects.dist(:,ndims));
      end
      title('Distance reject p-vals');
      drawnow;
      switch options.debug.showMmfPvals
        case -1
          pause;
        case -2
          keyboard;
      end
    end

    %save results
    initCoord(iImage).nSpots = nSpots;
    initCoord(iImage).allCoordPix = coordList;
    initCoord(iImage).amp = ampList;
    initCoord(iImage).bg = bgList;


    %check whether frame is empty
    if initCoord(iImage).nSpots == 0
      emptyFrames = [emptyFrames; iImage]; %#ok<AGROW>
      initCoord(iImage).allCoord = initCoord(iImage).allCoordPix;
    else
      % calc real space coordinates
      initCoord(iImage).allCoord = initCoord(iImage).allCoordPix .* repmat(pixelSize,initCoord(iImage).nSpots,2);
    end

    %display progress
    prog = kitProgress(iImage/length(goodImages),prog);
end

%% Post-processing

%sort list of empty frames, keep only unique frames
emptyFrames = unique(emptyFrames);

%store empty frames and frames where detection failed in structure
%exceptions
exceptions = struct('emptyFrames',emptyFrames);

% save results
initCoord(1).exceptions = exceptions;
initCoord(1).localMaxima = localMaxima;
dataStruct.dataProperties.psfSigma = psfSigma;
dataStruct.initCoord = initCoord;

dataStruct.failed = 0;

job.dataStruct{channel} = dataStruct;
