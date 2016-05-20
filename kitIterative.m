function job=kitIterative(job,reader,channel)
% KITITERATIVE Iterates on particle detection using existing tracks as guide.
%
% Assumes no merging/splitting.
%
% Copyright (c) 2015 Jonathan W. Armond

% Extract trackLists.
jobTmp = kitExtractTracks(job, channel);
trackList = jobTmp.dataStruct{channel}.trackList;
clear jobTmp

nFrames = job.metadata.nFrames;
options = job.options;
options.mmfAddSpots = 0; % Force adding new spots.
% Default minimum, 50% tracked.
minLength = round(0.5*nFrames);

ds = job.dataStruct{channel};

filterPrm = ds.dataProperties.FILTERPRM;
psfSigma = filterPrm([1 1 3])./job.metadata.pixelSize
roiWnd = psfSigma/2
minFitPixels = max(round(2*roiWnd-1))^3

% Sort tracks by length. Lengthen longest first.
stend = horzcat(ds.tracks.seqOfEvents);
stend = stend(:,1:4:end);
len = stend(2,:)-stend(1,:)+1;
[~,idx] = sort(len,'descend');
nTracks = length(ds.tracks);

warning('off','HISTOGRAM:notEnoughDataPoints');

% Storage for new particles.
new(1:nFrames) = struct('spots',[],'amps',[],'bg',[]);

% Algorithm:
% Find tracks with missing timepoints.
% Search through unassigned spots. If there is an unambiguous spot in predicted position, assign it.
% Else fit spot and check for one in place.

verbose=4
% Go over all tracks in sorted order.
for i=1:nTracks
  tidx = idx(i);
  tl = trackList(tidx);
  gaps = find(isnan(tl.featIndx));
  % If tracked fraction is below threshold, skip.
  if nFrames-length(gaps)<minLength
    continue
  end

  % Attempt to fill in each gap, going forward.
  for j=1:length(gaps)
    t=gaps(j);
    if t==1
      continue
    end
    prevIndx = tl.featIndx(t-1);
    % If not tracked in previous frame, can't fill in this frame.
    if isnan(prevIndx)
      continue
    end

    % Frame to fill in.
    img = kitReadImageStack(reader,job.metadata,t,channel,job.crop,0);
    [sx,sy,sz] = size(img);

    prevPos = ds.initCoord(t-1).allCoordPix(prevIndx,1:3);
    newPos = prevPos; % Assume no movement for now.
    if any(newPos<1) || any(newPos>[sy sx sz])
      continue
    end

    if verbose>3
      figure(1);
      % unassigned particles.
      assigned = false(ds.initCoord(t).nSpots,1);
      for k=1:length(trackList)
        if ~isnan(trackList(k).featIndx(t))
          assigned(trackList(k).featIndx(t)) = true;
        end
      end
      unassigned = ~assigned;
      subplot(1,2,1);
      prevImg = kitReadImageStack(reader,job.metadata,t-1,channel,job.crop,0);
      showSpots(prevImg,prevPos);
      subplot(1,2,2);
      showSpots(img,{prevPos,ds.initCoord(t).allCoordPix(unassigned,1:3),ds.initCoord(t).allCoordPix(assigned,1:3)});
      drawnow
    end


    % Try to detect particle in this frame.
    roiBnds = [newPos-roiWnd; newPos+roiWnd];
    roiBnds = [max(roiBnds(1,:),[1 1 1]); min(roiBnds(2,:),[sy sx sz])];
    r = round(roiBnds);
    roi = img(r(1,2):r(2,2),r(1,1):r(2,1),r(1,3):r(2,3));
    if numel(roi)>=minFitPixels
      if options.debug.showIterativeFrames
        figure(1);
        prevImg = kitReadImageStack(reader,job.metadata,t-1,channel,job.crop,0);
        % Trim img.
        proiBnds = [prevPos-roiWnd; prevPos+roiWnd];
        proiBnds = [max(proiBnds(1,:),[1 1 1]); min(proiBnds(2,:),[sy sx sz])];
        r = round(proiBnds);
        prevRoiMax = max(prevImg(r(1,2):r(2,2),r(1,1):r(2,1),r(1,3):r(2,3)),[],3);
        imshowpair(imresize(prevRoiMax,25),imresize(max(roi,[],3),25),'montage');
        drawnow;
      end

      cands = histcutSpots(roi,options,ds.dataProperties);
      if ~isempty(cands)
        [spots,amps,bg] = mixtureModelFit(cands,roi,psfSigma([1 3]),options);
        if ~isempty(spots)
          % Correct particles for ROI offset.
          spots(:,1:3) = bsxfun(@plus,spots(:,1:3),roiBnds(1,:));
          new(t).spots = [new(t).spots; spots];
          new(t).amps = [new(t).amps; amps];
          new(t).bg = [new(t).bg; bg];
          if options.debug.showIterativeFrames
            figure(1);
            clf;
            subplot(1,2,1);
            showSpots(prevImg,prevPos);
            subplot(1,2,2);
            showSpots(img,{spots,prevPos});
          end
        end
      end
    end
  end % t
end % i


% Check particles are new.
nNewCount = 0;
for t=1:nFrames
  if ~isempty(new(t).spots)
    prevSpots = ds.initCoord(t).allCoordPix;
%     dists = createDistanceMatrix(new(t).spots(:,1:3),prevSpots(:,1:3));
%     dists(logical(eye(size(dists)))) = inf;
%     valid = all(dists>psfSigma(1),2);
    distsXY = createDistanceMatrix(new(t).spots(:,1:2),prevSpots(:,1:2));
    distsXY(logical(eye(size(distsXY)))) = inf;
    distsZ = createDistanceMatrix(new(t).spots(:,3),prevSpots(:,3));
    distsZ(logical(eye(size(distsZ)))) = inf;
    valid = all(distsXY>psfSigma(1) & distsZ>psfSigma(3),2);

    if any(valid)
      % Store new particles.
      nValid = sum(valid);
      nNewCount = nNewCount + nValid;
      ds.initCoord(t).allCoordPix = [prevSpots; new(t).spots(valid,:)];
      newSpotsUm = bsxfun(@times, new(t).spots(valid,:), repmat(job.metadata.pixelSize,[1 2]));
      ds.initCoord(t).allCoord = [ds.initCoord(t).allCoord; newSpotsUm];
      ds.initCoord(t).nSpots = ds.initCoord(t).nSpots + sum(valid);
      ds.initCoord(t).amp = [ds.initCoord(t).amp; new(t).amps(valid,:)];
      ds.initCoord(t).bg = [ds.initCoord(t).bg; new(t).bg(valid,:)];

      if options.debug.showIterative
        showSpots(img,{prevSpots(:,1:3),new(t).spots(valid,1:3)});
        title(sprintf('t = %d, %d new particles',t,nValid));
        drawnow;
      end
    end
  end
end

% TODO Reverse extrapolation for beginning of tracks.


warning('on','HISTOGRAM:notEnoughDataPoints');
kitLog('%d particles added',nNewCount);

% Re-generate tracks.
kitLog('Refitting plane');
ds = kitFitPlane(job,reader,ds,channel,0);
kitLog('Regenerating tracks');
ds = kitGenerateTracks(ds);

% Store results.
job.dataStruct{channel} = ds;
