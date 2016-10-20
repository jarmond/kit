function dataStruct = kitAssembleNeighbourStructs(fullDataStruct,channel,options)
% Uses spots found in a neighbour channel to construct tracks and
% sisterList structures.
%
%
% Created by: C. A. Smith
% Copyright (c) 2016 C. A. Smith

%% Input + initialization

% create output structure
dataStruct = fullDataStruct{channel};

% get initCoord and planeFit from which new coordinates will be organised
initCoord = dataStruct.initCoord;
planeFit = dataStruct.planeFit;
% get reference tracks, trackList and sisterList
refTracks = fullDataStruct{options.coordSystemChannel}.tracks;
refTrackList = fullDataStruct{options.coordSystemChannel}.trackList;
refSisterList = fullDataStruct{options.coordSystemChannel}.sisterList;
nSisters = length(refSisterList);
nFrames = size(refSisterList(1).coords1,1);

% get trackIDs and spotIDs
trackIDs = refSisterList(1).trackPairs(:,1:2);
spotIDs = repmat({[]},nSisters,1);
for iSisPair = 1:nSisters
  spotIDs{iSisPair} = [refTrackList(trackIDs(iSisPair,1)).featIndx refTrackList(trackIDs(iSisPair,2)).featIndx];
end

% create tracks and sisterList structures for compiling neighbour data
tracks = refTracks;
sisterList = refSisterList;


%% Move data from initCoord to tracks and sisterList

% loop through sister pairs, then between each sister, then over time
for iSisPair = 1:nSisters
  for iSis = 1:2
    for iFrame = 1:nFrames
      
      % get spotID
      iSpotID = spotIDs{iSisPair}(iFrame,iSis);
      % get trackID
      iTrackID = trackIDs(iSisPair,iSis);
      
      % find time period of track sequence
      firstTime = tracks(iTrackID).seqOfEvents(1,1);
      lastTime = tracks(iTrackID).seqOfEvents(2,1);
      % if this frame is contained within the time period, transfer the
      % information
      if ismember(iFrame,firstTime:lastTime) && ~isnan(iSpotID)
        
        % adjust the origin of frames to the start of the track, where the
        % 8-factor comes from the 1x(40x8) structure of tracksCoordAmpCG
        adjFrame = 8*(iFrame-firstTime);
        
        % move data from initCoord to tracks
        try
        tracks(iTrackID).tracksCoordAmpCG(adjFrame+1:adjFrame+3) = initCoord(iFrame).allCoord(iSpotID,1:3);
        tracks(iTrackID).tracksCoordAmpCG(adjFrame+5:adjFrame+7) = initCoord(iFrame).allCoord(iSpotID,4:6);
        tracks(iTrackID).tracksCoordAmpCG([adjFrame+4 adjFrame+8]) = initCoord(iFrame).amp(iSpotID,1:2);
        
        tracks(iTrackID).coordAmp4Tracking(adjFrame+1:adjFrame+3) = planeFit(iFrame).rotatedCoord(iSpotID,1:3);
        tracks(iTrackID).coordAmp4Tracking(adjFrame+5:adjFrame+7) = planeFit(iFrame).rotatedCoord(iSpotID,4:6);
        tracks(iTrackID).coordAmp4Tracking([adjFrame+4 adjFrame+8]) = initCoord(iFrame).amp(iSpotID,1:2);
        catch
           qq=1; 
        end
        
      end
      
      % move data from initCoord to sisterList, checking a spotID exists
      if iSis==1
        if isnan(iSpotID)
          sisterList(iSisPair).coords1(iFrame,:) = NaN;
        else
          sisterList(iSisPair).coords1(iFrame,:) = planeFit(iFrame).rotatedCoord(iSpotID,:);
        end
      else
        if isnan(iSpotID)
          sisterList(iSisPair).coords2(iFrame,:) = NaN;
        else
          sisterList(iSisPair).coords2(iFrame,:) = planeFit(iFrame).rotatedCoord(iSpotID,:);
        end
      end
      
    end
  end
  
  % measure sister pair-specific vectors and distances
  sisterList(iSisPair).sisterVectors = sisterList(iSisPair).coords2(:,1:3) - sisterList(iSisPair).coords1(:,1:3);
  sisterList(iSisPair).distances(:,1) = sqrt(sum(sisterList(iSisPair).sisterVectors.^2,2));
  sisterList(iSisPair).distances(:,2) = sqrt(sum((sisterList(iSisPair).coords2(:,4:6) - sisterList(iSisPair).coords1(:,4:6)).^2,2));
  
end

% save results to final structure
dataStruct.tracks = tracks;
dataStruct.sisterList = sisterList;
