function job=kitIterative(job,reader,channel)
% KITITERATIVE Iterates on particle detection using existing tracks as guide.
%
% Assumes no merging/splitting.
%
% Copyright (c) 2015 Jonathan W. Armond

nFrames = job.metadata.nFrames;
options = job.options;
options.mmfAddSpots = 1; % Force adding new spots.
options.debug.showIterative=1;

ds = job.dataStruct{channel};

%movieInfo = genMovieInfo(ds);
ndims = 2 + job.metadata.is3D;
filterPrm = ds.dataProperties.FILTERPRM;
psfSigma = filterPrm([1 1 3]);
roiWnd = psfSigma*4;

% Sort tracks by length. Lengthen longest first.
stend = horzcat(ds.tracks.seqOfEvents);
stend = stend(:,1:4:end);
len = stend(2,:)-stend(1,:)+1;
[~,idx] = sort(len);
ds.tracks = ds.tracks(idx);
stend = stend(:,idx);
nTracks = length(ds.tracks);

warning('off','HISTOGRAM:notEnoughDataPoints');

% Forward iteration.
nNewCount = 0;
for t=3:nFrames
  % All tracks which end in previous frame and have at least one frame before that.
  idx = find(stend(2,:)==t-1 & stend(1,:)<=t-2);

  newSpots = [];
  newAmps = [];
  newBg = [];
  if ~isempty(idx)
    img = kitReadImageStack(reader,job.metadata,t,channel,job.crop,0);
    [sx,sy,sz] = size(img);

    for i=1:length(idx)
      % TODO: Kalman filter to predict position?

      % Previous two positions.
      featIndx = ds.tracks(i).tracksFeatIndxCG(end-1:end);
      pos = [ds.initCoord(t-2).allCoordPix(featIndx(1),1:3); ds.initCoord(t-1).allCoordPix(featIndx(2),1:3)];

      % Extrapolate to current position, assuming continuation of previous motion.
      vec = pos(2,:) - pos(1,:);
      newPos = pos(2,:) + vec;
      if any(newPos<1) || any(newPos>[sx sy sz])
        continue
      end

      % Try to detect particle in this frame.
      roiBnds = [newPos-roiWnd; newPos+roiWnd];
      roiBnds = [max(roiBnds(1,:),[1 1 1]); min(roiBnds(2,:),[sx sy sz])];
      r = round(roiBnds);
      roi = img(r(1,1):r(2,1),r(1,2):r(2,2),r(1,3):r(2,3));
      if numel(roi)>=9
        cands = histcutSpots(roi,options,ds.dataProperties);
        if ~isempty(cands)
          [spots,amps,bg] = mixtureModelFit(cands,roi,psfSigma([1 3]),options);
          if ~isempty(spots)
            % Correct particles for ROI offset.
            spots(:,1:3) = bsxfun(@plus,spots(:,1:3),roiBnds(1,:));
            newSpots = [newSpots; spots];
            newAmps = [newAmps; amps];
            newBg = [newBg; bg];
          end
        end
      end
    end

    if ~isempty(newSpots)
      % Check particles are new.
      prevSpots = ds.initCoord(t).allCoordPix;
      dists = createDistanceMatrix(newSpots(:,1:3),prevSpots(:,1:3));
      dists(logical(eye(size(dists)))) = inf;
      valid = all(dists>psfSigma(1),2);

      if any(valid)
        % Store new particles.
        nValid = sum(valid);
        nNewCount = nNewCount + nValid;
        ds.initCoord(t).allCoordPix = [prevSpots; newSpots(valid,:)];
        newSpotsUm = bsxfun(@times, newSpots(valid,:), repmat(job.metadata.pixelSize,[1 2]));
        ds.initCoord(t).allCoord = [ds.initCoord(t).allCoord; newSpotsUm];
        ds.initCoord(t).nSpots = ds.initCoord(t).nSpots + sum(valid);
        ds.initCoord(t).amp = [ds.initCoord(t).amp; newAmps(valid,:)];
        ds.initCoord(t).bg = [ds.initCoord(t).bg; newBg(valid,:)];

        if options.debug.showIterative
          showSpots(img,{prevSpots(:,1:3),newSpots(valid,1:3)});
          title(sprintf('t = %d, %d new spots',t,nValid));
          drawnow;
        end
      end
    end
  end
end

% Reverse extrapolation for beginning of tracks.


warning('on','HISTOGRAM:notEnoughDataPoints');

% Re-generate tracks.
kitLog('Refitting plane');
ds = kitFitPlane(job,reader,ds,channel,0);
kitLog('Regenerating tracks');
ds = kitGenerateTracks(ds);

% Store results.
job.dataStruct{channel} = ds;
