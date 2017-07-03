function job=kitLocalIntensity(job, reader, metadata, channel, opts)
% KITLOCALINTENSITY Measure intensity local to kinetochores
%
%

%% Preparation: get metadata, produce structures

% get some metadata
nFrames = metadata.nFrames;

% check if this channel contains coordinate or track information, or neither
if length(job.dataStruct)<channel || isempty(job.dataStruct{channel})
  refChan = job.options.coordSystemChannel;
else
  refChan = channel;
end
chrShift = job.options.chrShift.result{job.options.coordSystemChannel,channel};
chrShift = chrShift(1:3)./metadata.pixelSize(1:3);
initCoord = job.dataStruct{refChan}.initCoord;
planeFit = job.dataStruct{refChan}.planeFit;
if strcmp(job.options.jobProcess,'zonly') && ~isfield(job.dataStruct{refChan},'trackList')
  useTracks = 0;
  nKTs = initCoord(1).nSpots;
else
  useTracks = 1;
  trackList = job.dataStruct{refChan}.trackList;
  nKTs = length(trackList);
end

% predesignate intensity structure
intStruct(1:nKTs) = struct(...
    'intensity',nan(nFrames,1),...
    'intensity_median',nan(nFrames,1),...
    'intensity_min',nan(nFrames,1),...
    'intensity_max',nan(nFrames,1),...
    'intensity_ratio',nan(nFrames,1),...
    'maskCoord',nan(nFrames,3),...
    'maxCoord',nan(nFrames,3),...
    'angleToMax',nan(nFrames,1),...
    'distToMax',nan(nFrames,1));

% convert pole shift to pixels
poleShiftPixels = round(opts.poleShift / metadata.pixelSize(1));

% Gaussian filter
hgauss = fspecial('gauss');


%% Create mask
r = round(opts.maskRadius / metadata.pixelSize(1));
if ~strcmp(opts.maskShape, 'cone')
  se = strel('disk', r, 0);
  asymMask = 0;
  mask = double(se.getnhood());
  mask(mask == 0) = nan;
  if strcmp(opts.maskShape, 'semicirc')
    mask(r+2:end,:) = nan;
    asymMask = 1;
  end
else
  mask = ones(r);
  % Make conical mask.
  ycut = tand(opts.maskConeAngle)*(0:r-1);
  for i=1:r
    ycutFlr = floor(ycut(i));
    if ycutFlr>=1
      mask(1:ycutFlr,i) = nan;
    end
  end
  % Add asymmetry and orient.
  mask = [flipud(mask); nan(r-1,r)];
  % Mirror
  mask = [fliplr(mask) mask(:,2:r)];
  asymMask = 1;
end

if asymMask == 1
  lMask = mask; % For -ve x pole attached.
  rMask = flipud(mask); % For +ve x pole attached.
end

maskWarning=0;


%% Get intensities

prog = kitProgress(0);
for t=1:nFrames
    stack = kitReadImageStack(reader, metadata, t, channel, job.crop, 0);
    if opts.gaussFilterSpots
      % Gaussian filter each z-plane.
      for z=1:size(stack,3)
        stack(:,:,z) = imfilter(stack(:,:,z), hgauss, 'replicate');
      end
    end

    if t==1
      % Take overall fluorescence estimate.
      back = mean(stack(:));

      % Estimate background and signal.
      [backMode,sigVal,imgHist] = estBackgroundAndSignal(stack);
      intensityDistF(:) = imgHist(:,2);
      intensityDistX(:) = imgHist(:,1);

      % Estimate background as difference between signal mode and background mode.
      backDiff = sigVal - backMode;
    end
    
    for j=1:nKTs
      if useTracks
        % Map track back to coords. FIXME use trackList instead
        pixIdx = trackList(j).featIndx(t);
      else
        pixIdx = j;
      end
        
      if pixIdx > 0 && pixIdx < size(initCoord(t).allCoordPix,1)
        % Read off intensity in radius r around spot.
        pixCoords = initCoord(t).allCoordPix(pixIdx,1:3);

        if sum(isnan(pixCoords)) == 0

          if useTracks && trackList(j).attach ~= 0 && ...
                size(pixCoords,2)==size(planeFit(t).planeVectors,2) && ...
                refChan ~= channel
            % Shift mask toward pole if measuring intensity of
            % non-localised channel, if known sister.
            maskCoords = pixCoords + [sign(trackList(j).attach),0,0] * ...
                poleShiftPixels*planeFit(t).planeVectors;
          else
            maskCoords = pixCoords;
          end
          % chromatic shift from coordSysChan to the channel being measured
          maskCoords = maskCoords + chrShift;

          % Extract intensity.
          x = max(1,round(maskCoords(1)));
          y = max(1,round(maskCoords(2)));
          z = max(1,round(maskCoords(3)));
          if z > size(stack,3)
            warning('Spot outside frame boundaries');
            continue
          end
          if size(stack,3)>1
            imgPlane = stack(:,:,z);
          else
            imgPlane = stack;
          end
          [my,mx] = size(imgPlane);
          % If too close to edge, skip.
          if x <= r || y <= r || x >= mx-r || y >= my-r
            continue
          end

          % If asymmetric mask, choose left/right mask based on pole
          % attachment.
          if useTracks && asymMask
            if trackList(j).attach > 0
              mask = rMask;
            else
              mask = lMask;
            end
            if ~isempty(planeFit(t).planeVectors)
              % Rotate mask into coordinate system.
              angle = acos(planeFit(t).planeVectors(1,1))*180/pi;
              mask = imrotate(mask, angle, 'nearest', 'crop');
            else
              if ~maskWarning
                warning('No rotation vector for mask');
                maskWarning=1;
              end
            end
          end
          % Coordinates are in image system, so plot(x,y) should draw the
          % spots in the correct place over the image. However, the image
          % matrix is indexed by (row,col) => (y,x).
          maskImg = mask .* imgPlane(y-r:y+r, x-r:x+r);
          nonNanPix = maskImg(~isnan(maskImg));
          intStruct(j).intensity(t) = mean(nonNanPix);
          intStruct(j).intensity_median(t) = median(nonNanPix);
          intStruct(j).intensity_max(t) = max(nonNanPix);
          intStruct(j).intensity_min(t) = min(nonNanPix);
          intStruct(j).intensity_ratio(t) = ...
              intStruct(j).intensity_max(t)/intStruct(j).intensity_min(t);
          intStruct(j).maskCoord(t,:) = maskCoords;
          %[maxX,maxY] = ind2sub(size(maskImg),maxIdx);
          [~,maxIdx] = max(maskImg(:));
          [maxY,maxX] = ind2sub(size(maskImg),maxIdx);
          intStruct(j).maxCoord(t,:) = [maxX+x-r-1,maxY+y-r-1,z];

          % Calculate angle between spot and max point.
          vector = pixCoords-intStruct(j).maxCoord(t,:);
          intStruct(j).distToMax(t) = norm(vector(1:2),2)*metadata.pixelSize(1);
          % Rotate angle into coordinate system.
          if size(pixCoords,2)==size(planeFit(t).planeVectors,2)
            vector = vector*planeFit(t).planeVectors;
          end
          angle = -atan2(vector(2),-vector(1)); % See doc atan2 for diagram.
          intStruct(j).angleToMax(t) = angle;
          intStruct(j).referenceChannel = refChan;
          if useTracks
            intStruct(j).referenceStruct = 'trackList';
          else
            intStruct(j).referenceStruct = 'initCoord';
          end
        end
      end
    end

  % Report progress.
  prog = kitProgress(t/nFrames, prog);
end

if ~isempty(opts.photobleachCorrect)
  % Compute photobleach from entire image.
  pbProfile=kitIntensityDistn(job,reader,metadata,channel,[],[],1,job.crop);
  t=((1:size(pbProfile,1))-1)';

  % Normalize PB if necessary.
  if pbProfile(1) ~= 1 || ~all(pbProfile(:) < 1)
    pbProfile(:) = pbProfile(:)/pbProfile(1);
  end

  if numel(t) > 4 && license('test','Curve_Fitting_Toolbox')
    % Fit double exp.
    pbF = fit(t,pbProfile(:,channel),'exp2');
    % Correct.
    for j=1:nKTs
      intStruct(j).intensity(:) = intStruct(j).intensity(:)./pbF(t);
      intStruct(j).intensity_median(:) = intStruct(j).intensity_median(:)./pbF(t);
      intStruct(j).intensity_max(:) = intStruct(j).intensity_max(:)./pbF(t);
      intStruct(j).intensity_min(:) = intStruct(j).intensity_min(:)./pbF(t);
    end
  else
    kitLog('Warning: Not correcting for photobleach');
  end

end

% record background intensity
cellInt.back = back;
cellInt.backMode = backMode;
cellInt.backDiff = backDiff;
cellInt.maskPixelRadius = r;
cellInt.intensityDistF = intensityDistF;
cellInt.intensityDistX = intensityDistX;

% return data
job.dataStruct{channel}.cellInt = cellInt;
job.dataStruct{channel}.spotInt = intStruct;

end