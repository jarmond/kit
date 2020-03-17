function sigma_through_time = jonathanAnalyseSpotShape(job,varargin)
%
%Estimate spot shape for tracked pairs at all frames through movie
%Uses jonathanSingleSpotEstimateSigma
%
%Jonathan U Harrison 2020-01-30
%%%%%%%%%%%

% set default options
opts.contrast = [0.1 1];
opts.channel = 1;
opts.newFig = 0;
opts.sigmaScale = 0.8;
opts.title = [];
opts.transpose = 0;
opts.withinFig = 0;
opts.zoomScale = 1;
opts.zoom = 1;
opts.zoomRangeMicrons = 0.5;
opts.verbose = 0;
opts.everyKFrames = 1;

% process options
opts = processOptions(opts, varargin{:});

[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),...
    job.metadata,0);

chan = opts.channel;
nFrames = job.metadata.nFrames;
sisterList = job.dataStruct{chan}.sisterList;

%%%%%%%%
% get pixel resolution
pixelSize = job.metadata.pixelSize;
nPairs = size(sisterList(1).trackPairs,1);
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,chan}(1:3);

% all_coords = nan(2,3,nPairs); %coords x sister
% for sisPair = 1:nPairs
%     trackIDs = sisterList(1).trackPairs(sisPair,1:2);
%     for iSis = 1:2
%         tk = trackIDs(iSis);
%         track = job.dataStruct{chan}.tracks(tk);
%
%         startTime = track.seqOfEvents(1,1);
%         endTime   = track.seqOfEvents(2,1);
%         if timePoint < startTime || timePoint > endTime
%             all_coords(iSis,:,sisPair) = nan(1,3);
%         else
%             all_coords(iSis,:,sisPair) = ...
%                 track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
%             all_coords(iSis,:,sisPair) = all_coords(iSis,:,sisPair) + chrShift;
%             all_coords(iSis,:,sisPair) = all_coords(iSis,:,sisPair)./pixelSize;
%         end
%     end
% end
%%%%%%%%

%note that all time points are stored, but only a subset are computed
sigma_through_time = nan(nFrames,3,2,nPairs); %frames by XYZ by sister
all_coords = nan(nFrames,3,2,nPairs); %coords x sister
amplitude = nan(nFrames,1,2,nPairs);
% accumulate track information by channel and sister
for timePoint = 1:opts.everyKFrames:nFrames
% for a given frame, need to first get coordinates of all spots to use to fit each individual spot    
    for sisPair = 1:nPairs
        % get track information
        trackIDs = sisterList(1).trackPairs(sisPair,1:2);
        for iSis = 1:2
            tk = trackIDs(iSis);
            track = job.dataStruct{chan}.tracks(tk);
            
            startTime = track.seqOfEvents(1,1);
            endTime   = track.seqOfEvents(2,1);
            if timePoint < startTime || timePoint > endTime
                all_coords(timePoint,:,iSis,sisPair) = nan(1,3);
            else
                all_coords(timePoint,:,iSis,sisPair) = ...
                    track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
                all_coords(timePoint,:,iSis,sisPair) = all_coords(timePoint,:,iSis,sisPair) + chrShift;
                all_coords(timePoint,:,iSis,sisPair) = all_coords(timePoint,:,iSis,sisPair)./pixelSize;
                amplitude(timePoint,:,iSis,sisPair) = track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-4);
            end
        end
    end
    fprintf('Time point %d\n',timePoint);
    for sisPair = 1:nPairs
        trackIDs = sisterList(1).trackPairs(sisPair,1:2);
        for iSis = 1:2
            tk = trackIDs(iSis);
            track = job.dataStruct{chan}.tracks(tk);
            
            startTime = track.seqOfEvents(1,1);
            endTime   = track.seqOfEvents(2,1);
            if timePoint < startTime || timePoint > endTime
                sigma_through_time(timePoint,:,iSis,sisPair) = nan(1,3);
            else
                sigma_through_time(timePoint,:,iSis,sisPair) = ...
                    diag(jonathanSingleSpotEstimateSigma(job,md,reader,...
                    all_coords, amplitude, ...
                    'sisterPair',sisPair,'timePoint',timePoint,...
                    'sigmaScale',1.5,'verbose',opts.verbose, ...
                    'contrast',opts.contrast, ...
                    'zoomRangeMicrons',opts.zoomRangeMicrons, ...
                    'sisID',iSis));
                %error('enough so stop now');
            end
        end
    end
end
