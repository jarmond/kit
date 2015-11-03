function dataStruct = kitGenerateTracks(dataStruct,verbose)
% KITGENERATETRACKS Tracks kinetochores throughout a movie
%
%    DATASTRUCT = KITGENERATETRACKS(DATASTRUCT) Connect spots to generate
%    kinetochore tracks. DATASTRUCT should be output from KITFITPLANE. Adds
%    field .tracks.
%
% Copyright (c) 2007 K. Jaqaman
% Copyright (c) 2012 Jonathan W. Armond

if nargin<2
  verbose=0;
end

movieInfo = genMovieInfo(dataStruct);

%get tracking parameters
gapCloseParam = dataStruct.dataProperties.tracksParam.gapCloseParam;
costMatrices = dataStruct.dataProperties.tracksParam.costMatrices;
kalmanFunctions = dataStruct.dataProperties.tracksParam.kalmanFunctions;

%track the kinetochores
tracks = trackCloseGapsKalman(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,3,0);

%replace the coordinate used for tracking (whether the rotated
%coordinates or the original coordinates shifted by center of mass) by
%the original coordinates
for iTrack = 1 : length(tracks)

    %store coordinates used for tracking in another field
    tracks(iTrack).coordAmp4Tracking = tracks(iTrack).tracksCoordAmpCG;

    %fetch the start and end time of this track
    startTime = tracks(iTrack).seqOfEvents(1,1);
    endTime = tracks(iTrack).seqOfEvents(2,1);

    %go over all frames where this track exists
    for iFrame =  startTime : endTime

        %get the feature making up this track in this frame
        iFeature = tracks(iTrack).tracksFeatIndxCG(iFrame-startTime+1);

        %if there is a feature (not a gap)
        if iFeature ~= 0

            %replace coordiantes and their stds
            tracks(iTrack).tracksCoordAmpCG(1,(iFrame-startTime)*8+1:...
                (iFrame-startTime)*8+3) = dataStruct.initCoord(iFrame).allCoord(iFeature,1:3);
            tracks(iTrack).tracksCoordAmpCG(1,(iFrame-startTime)*8+5:...
                (iFrame-startTime)*8+7) = dataStruct.initCoord(iFrame).allCoord(iFeature,4:6);

        end

    end

end %(for iTrack = 1 : numTracks)


%store tracks in dataStruct
dataStruct.tracks = tracks;

if verbose
  % Plot all tracks.
  trackStats = catStruct(3,'dataStruct.tracks.seqOfEvents');

  figure;
  hold on;
  trackLengths = zeros(length(tracks),1);
  for j = 1:length(tracks) % loop cols
    % read track coordinates. coordinates are unrotated/translated
    colCoords = trackData(j,dataStruct,trackStats,1);

    % plot individual tracks
    plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),...
          'Color',extendedColors(j))

    startTime = dataStruct.tracks(j).seqOfEvents(1,1);
    endTime = dataStruct.tracks(j).seqOfEvents(2,1);
    trackLengths(j) = endTime-startTime;
    if verbose >= 2
      % Label points by sequence
      labels = cellstr(num2str((startTime:endTime)'));
      text(colCoords(:,1),colCoords(:,2),colCoords(:,3),labels,'verticalalignment','bottom','horizontalalignment','right');
    end
  end
  hold off;
  fprintf('Mean track length: %g\n',mean(trackLengths));
end
