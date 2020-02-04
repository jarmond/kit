function sigma1 = jonathanSingleSpotEstimateSigma(job,md,reader,varargin)
%
%Compare also JonathanEstimateSigma and kitShowSisterPair
%
%Jonathan U Harrison 2020-01-22
%%%%%%%%%%%


if nargin<1
    error('Must supply a job.');
end

% set default options
opts.contrast = [0.1 1];
opts.channel = 1;
opts.newFig = 0;
opts.sigmaScale = 0.8;
opts.sisterPair = 1;
opts.timePoint = 10;
opts.title = [];
opts.transpose = 0;
opts.withinFig = 0;
opts.zoomScale = 1;
opts.zoom = 1;
opts.zoomRangeMicrons = 0.5;
opts.verbose = 0;

% process options
opts = processOptions(opts, varargin{:});

%% Image and coordinate acquisition

if nargin < 3
    % open movie
    [md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata,0);
end

% get coordinate system and plot channels
chan = opts.channel;
refChan = job.options.coordSystemChannel;
% get chromatic shift
chrShift = job.options.chrShift.result{refChan,chan}(1:3);

% get how close in to zoom to
zoomRangeMicrons = opts.zoomRangeMicrons;

% get sister information
sisPair = opts.sisterPair;
sisterList = job.dataStruct{chan}.sisterList;

if sisPair > length(sisterList)
    kitLog('Sister pair ID provided too large. Only have %i sisters in movie. Quitting.',length(sisterList))
    return
end

% get track information
trackIDs = sisterList(1).trackPairs(sisPair,1:2);
timePoint = opts.timePoint;

% get pixel resolution
pixelSize = job.metadata.pixelSize;

coords = nan(2,3); %coords x sister

% accumulate track information by channel and sister
for iSis = 1:2
    tk = trackIDs(iSis);
    track = job.dataStruct{chan}.tracks(tk);
    
    startTime = track.seqOfEvents(1,1);
    endTime   = track.seqOfEvents(2,1);
    if timePoint < startTime || timePoint > endTime
        coords(iSis,:) = nan(1,3);
    else
        coords(iSis,:) = ...
            track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
        coords(iSis,:) = coords(iSis,:) + chrShift;
        coords(iSis,:) = coords(iSis,:)./pixelSize;
    end
end

if sum(isnan(coords(:))) == 6
    kitLog('No coordinates found for sister pair %i at time point %i. Ignoring.',sisPair,timePoint);
    sigma1 = nan(3,3);
    return
end


% calculate pair centre and convert to pixels
% centrePxl = nanmean(coords);
% centrePxl = round(centrePxl);
if ~any(isnan(coords(1,:)))
    centrePxl = round(coords(1,:)); %choose sister 1 (check if is nan)
else 
    centrePxl = round(coords(2,:)); %already checked if both are nan so sister 2 should not be nan
end

%% Image production

% read stack
img = kitReadImageStack(reader,md,timePoint,chan,job.ROI.crop,0);

xReg = [centrePxl(1)-ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(1)))+1 ...
    centrePxl(1)+ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(1)))+1];
yReg = [centrePxl(2)-ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(2)))+1 ...
    centrePxl(2)+ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(2)))+1];
zReg = [centrePxl(3)-ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(3)))+1 ...
    centrePxl(3)+ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(3)))+1];
if opts.transpose
    imgCrpd = img(xReg(1):xReg(2),yReg(1):yReg(2),zReg(1):zReg(2));
else
    imgCrpd = img(yReg(1):yReg(2),xReg(1):xReg(2),zReg(1):zReg(2));
end
if opts.verbose
    figure; subplot(1,2,1);
    imshow(max(imgCrpd,[],3),[],'InitialMagnification',1000);
    title('Max projection')
    subplot(1,2,2);
    imshow(imgCrpd(:,:,round(size(imgCrpd,3)/2)),[],'InitialMagnification',1000);
    title('Median section');
    
    if size(imgCrpd,3)<10
        figure;
        for i =1:size(imgCrpd,3)
            subplot(1,size(imgCrpd,3),i);
            imshow(imgCrpd(:,:,i),opts.contrast,'InitialMagnification',1000);
            title(sprintf('z=%d',i));
        end
    end
end
[mu,sigma1] = jonathanEstimateSigma(imgCrpd,opts.verbose,0,opts.sigmaScale);

end
