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

% process options
opts = processOptions(opts, varargin{:});

[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),...
    job.metadata,0);

chan = opts.channel;
nFrames = job.metadata.nFrames;
sisterList = job.dataStruct{chan}.sisterList;

nPairs = size(sisterList(1).trackPairs,1);
sigma_through_time = zeros(nFrames,3,2,nPairs); %frames by XYZ by sister

for sisPair =1:nPairs
    % get track information
    trackIDs = sisterList(1).trackPairs(sisPair,1:2);
    
    % accumulate track information by channel and sister
    for timePoint = 1:nFrames
        for iSis = 1:2
            tk = trackIDs(iSis);
            track = job.dataStruct{chan}.tracks(tk);
            
            startTime = track.seqOfEvents(1,1);
            endTime   = track.seqOfEvents(2,1);
            if timePoint < startTime || timePoint > endTime
                sigma_through_time(timePoint,:,iSis) = nan(1,3);
            else
                sigma_through_time(timePoint,:,iSis,sisPair) = ...
                    diag(jonathanSingleSpotEstimateSigma(job,md,reader,...
                    'sisterPair',sisPair,'timePoint',timePoint,...
                    'sigmaScale',1.5));
            end
        end
    end
end
