function dataStruct=kitDiagnostics(dataStruct,elapsed)
% KITDIAGNOSTICS Compute some basic tracking diagnostics
%
% Copyright (c) 2013 Jonathan W. Armond

diag.elapsedTime = elapsed;

% Spots.
diag.nSpotsPerFrame = mean(vertcat(dataStruct.initCoord.nSpots));

% Track counts.
diag.nSisters = length(dataStruct.sisterList);
if diag.nSisters == 1 && isempty(dataStruct.sisterList(1).coords1)
  diag.nSisters = 0;
end
diag.nTracks = length(dataStruct.trackList);
if diag.nTracks == 1 && isempty(dataStruct.trackList(1).coords)
  diag.nTracks = 0;
end

% Track lengths.
if diag.nSisters > 0
  coords = horzcat(dataStruct.sisterList.coords1);
  diag.avgSisterTrackLength = nanmean(sum(~isnan(coords(:,1:6:end))));
else
  diag.avgSisterTrackLength = 0;
end
if diag.nTracks > 0
  coords = horzcat(dataStruct.trackList.coords);
  diag.avgTrackLength = nanmean(sum(~isnan(coords(:,1:6:end))));
else
  diag.avgTrackLength = 0;
end

% Plane fit.
nFrames = length(dataStruct.planeFit);
nFramesWithPlane = size(vertcat(dataStruct.planeFit.plane),1);
diag.percentWithPlane = 100*nFramesWithPlane/nFrames;

% Sister track displacement mean and variance.
if diag.nSisters > 0
  coords1 = horzcat(dataStruct.sisterList.coords1);
  coords2 = horzcat(dataStruct.sisterList.coords2);
  % Displacements.
  coords = diff([coords1(:,1:6:end) coords2(:,1:6:end)]);
  diag.sisterDisp  = nanmean(nanmean(abs(coords)));
  diag.sisterVar = nanmean(nanvar(coords));
  diag.sisterPoints = sum(sum(~isnan(coords1(:,1:6:end))));
else
  diag.sisterDisp = 0;
  diag.sisterVar = 0;
  diag.sisterPoints = 0;
end

% Number of sister tracks with 75% length.
if diag.nSisters > 0
  nFrames = size(coords,1);
  cutOff = floor(0.75*nFrames);
  coords = horzcat(dataStruct.sisterList.coords1);
  diag.nLongSisters = sum(sum(~isnan(coords(:,1:6:end))) > cutOff);
else
  diag.nLongSisters = 0;
end



% Store result.
dataStruct.diagnostics = diag;
