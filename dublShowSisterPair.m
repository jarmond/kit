function dublShowSisterPair(job,varargin)
% DUBLSHOWSISTERPAIR Plots image of dual-channel movie with coordinates
% for a given sister pair.
%
%    DUBLSHOWSISTERPAIR(JOB,...) Plots coordinates in two channels over
%    the movie image for a given sister pair at a given timepoint.
%
%    Options, defaults in {}:-
%
%    channelMap: {[2 1 3]} or some perturbation. Order in which the titled
%    channels are presented in the figure, where typically:
%    1=red, 2=green, 3=blue.
%
%    plotChannels: {[1 2]} or some subset of [1 2 3]. Vector of channels
%    for plotting, where typically:
%    1=red, 2=green, 3=blue.
%
%    plotTitles: {'eGFP','tagRFP','mCate'} or similar. Titles of each plot
%    denoting which channel is which.
%
%    sisterPair: {1} or number. Sister pair within JOB being plotted.
%
%    timePoint: {1} or number. Timepoint at which to plot the sister pair.
%
%    transpose: {0} or 1. Whether to tranpose the image.
%
%    zoomScale: {1} or number. Magnification of a zoomed image as a
%    proportion of the default, so that 0.5 will zoom out, and 1.5 will
%    zoom in.
%
%    zoomTrack: 0 or {1}. Whether or not to zoom into the specific sister
%    pair. A value of 0 plots the whole cell, with just the coordinates of
%    the sister pair plotted.
%
%    zProject: {0}, 1 or -1. Whether or not to project in the z-direction.
%    -1 will project in the 3 z-slices surrounding the sister pair.
%
% Copyright (c) 2014 C. A. Smith

if nargin<1
  error('Must supply a job.');
end

% set default options
opts.channelMap = [2 1 3]; % green, red, blue
opts.subpixelate = 9;
opts.contrast = repmat({[0.1 1]},1,2);
opts.newFig = 0;
opts.plotChannels = 1:2;
opts.plotTitles = {'eGFP','tagRFP','mCate'};
opts.sisterPair = 1;
opts.timePoint = 1;
opts.transpose = 0;
opts.zoomScale = 1;
opts.zoomTrack = 1;
opts.zProject = 0;

% process options
opts = processOptions(opts, varargin{:});
while length(opts.plotTitles) < 3;
    opts.plotTitles = [opts.plotTitles,''];
end

%% Image and coordinate acquisition

% open movie
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie),job.metadata);

% get coordinate system and plot channels
coordSysChan = job.options.coordSystemChannel;
plotChans = opts.plotChannels;
nChans = length(plotChans);

% get sister information
sisPair = opts.sisterPair;
sisterList = job.dataStruct{coordSysChan}.sisterList;

% get track information
trackIDs = sisterList(1).trackPairs(sisPair,1:2);
timePoint = opts.timePoint;

% get pixel resolution
pixelSize = job.metadata.pixelSize;

trackCoord = nan(nChans,3,2);

% accumulate track information by channel and sister
for c = plotChans
    for i = 1:2
        tk = trackIDs(i);
        track{c,tk} = job.dataStruct{c}.tracks(tk);
        trackList{c,tk} = job.dataStruct{c}.trackList(tk);
        
        startTime = track{c,tk}.seqOfEvents(1,1);
        if timePoint < startTime
            trackCoord(c,:,i) = nan(1,3);
        else
            trackCoord(c,:,i) = ...
                track{c,tk}.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
        end
    end
end

% calculate pair centre - check whether coordinate system is reference
if sum(plotChans == coordSysChan) > 0
    
    centrePoint{coordSysChan} = nanmean(trackCoord(coordSysChan,:,:),3);
    
    for c = plotChans(plotChans ~= coordSysChan)
        centrePoint{c} = centrePoint{coordSysChan} + job.options.chrShift.result{coordSysChan,c}(:,1:3);
    end

% if not, then use first available channel
else
    
    newRefChan = min(plotChans);
    centrePoint{newRefChan} = nanmean(trackCoord(newRefChan,:,:),3);
    
    for c = plotChans(plotChans ~= newRefChan)
        centrePoint{c} = centrePoint{newRefChan} + job.options.chrShift.result{coordSysChan,c}(:,1:3);
    end
    
end
% convert each centre coordinate to pixels
for c = plotChans
    centrePxl{c} = round(centrePoint{c}./pixelSize);
end

mapChans = opts.channelMap;
if opts.transpose
    dims = [1 2];
else
    dims = [2 1];
end

%% RGB image production

% produce image
rgbImg = zeros([job.cropSize(dims), 3]);
rgbCrpd= zeros([2*ceil(opts.zoomScale*(2/pixelSize(1)))+1 ...
                2*ceil(opts.zoomScale*(2/pixelSize(1)))+1 ...
                3]);
rgbImgCS = zeros([job.cropSize(dims)*opts.subpixelate, 3]);
rgbCrpdCS= zeros([opts.subpixelate*2*ceil(opts.zoomScale*(2/pixelSize(1)))+1 ...
                opts.subpixelate*2*ceil(opts.zoomScale*(2/pixelSize(1)))+1 ...
                3]);

% produce raw image
for c = plotChans
    % read stack
    img = kitReadImageStack(reader,md,timePoint,c,job.crop,0);
%     img([1 2]) = img(dims);
    % max project over three z-slices around point
    if opts.zProject == 1
        img = max(img, [], 3);
    elseif opts.zProject == -1
        img = max(img(:,:,centrePxl{c}(3)-1:centrePxl{c}(3)+1), [], 3);
    else
        img = img(:,:,centrePxl{c}(3));
    end
    if opts.transpose
        img = img';
    end
    
    % define contrast stretch
    irange(c,:)=stretchlim(img,opts.contrast{c});
    % contrast stretch
    rgbImg(:,:,mapChans(c)) = imadjust(img', irange(c,:), []);
    
end

% produce chromatic shifted image
first = 1;
removed = zeros(3,2,2);
for c = plotChans
    if c ~= coordSysChan && first
        
        [img1,img2,tempRem] = ... 
            chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(c)),job.options.chrShift.result{coordSysChan,c}, ...
            'subpixelate',opts.subpixelate);
        
        irange(coordSysChan,:) = stretchlim(img1,opts.contrast{coordSysChan});
        irange(c,:) = stretchlim(img2,opts.contrast{c});
        img1 = imadjust(img1, irange(coordSysChan,:), []);
        img2 = imadjust(img2, irange(c,:), []);
        
        rgbImgCS(:,:,mapChans(coordSysChan)) = img1;
        rgbImgCS(:,:,mapChans(c)) = img2;
        
        % removed in form (coords,[start,end],channel)
        removed(1:2,:,coordSysChan) = tempRem(1:2,:,1);
        removed(1:2,:,c) = tempRem(1:2,:,2);
        
        first = 0;
    
    elseif c ~= coordSysChan
        [~,img2,tempRem] = ...
            chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(c)),job.options.chrShift.result{coordSysChan,c}, ...
            'interpol',opts.subpixelate);
        
        irange(c,:) = stretchlim(img2,[0.1 1]);
        img2 = imadjust(img2, irange(c,:), []);
        
        rgbImgCS(:,:,mapChans(c)) = img2;
        
        % removed in form (coords,[start,end],channel)
        removed(1:2,:,c) = tempRem(1:2,:,2);
        
    end
end

% produce cropped image around track centre
if opts.zoomTrack
    
    xReg = zeros(2,nChans); xRegCS = zeros(2,nChans);
    yReg = zeros(2,nChans); yRegCS = zeros(2,nChans);
    
    for c = plotChans
        
        xReg(:,c) = [centrePxl{c}(1)-ceil(opts.zoomScale*(2/pixelSize(1)))+1 ...
                   centrePxl{c}(1)+ceil(opts.zoomScale*(2/pixelSize(1)))+1];
        yReg(:,c) = [centrePxl{c}(2)-ceil(opts.zoomScale*(2/pixelSize(2)))+1 ...
                   centrePxl{c}(2)+ceil(opts.zoomScale*(2/pixelSize(2)))+1];
        
        % adjust chrShift regions for the amount of image removed during
        % chrShifting of the image itself
        xRegCS(:,c) = xReg(:,c)*opts.subpixelate - removed(1,1,c);
        yRegCS(:,c) = yReg(:,c)*opts.subpixelate - removed(2,1,c);
            
        if opts.transpose
            rgbCrpd(:,:,mapChans(c)) = rgbImg(yReg(1,c):yReg(2,c),xReg(1,c):xReg(2,c),mapChans(c));
            rgbCrpdCS(:,:,mapChans(c)) = rgbImgCS(yRegCS(1,c):yRegCS(2,c),xRegCS(1,c):xRegCS(2,c),mapChans(c));
        else
            rgbCrpd(:,:,mapChans(c)) = rgbImg(xReg(1,c):xReg(2,c),yReg(1,c):yReg(2,c),mapChans(c));
            rgbCrpdCS(:,:,mapChans(c)) = rgbImgCS(xRegCS(1,c):xRegCS(2,c),yRegCS(1,c):yRegCS(2,c),mapChans(c));
        end
        
        % adjust back the regions for plotting of the coordinates in the
        % next step
        xRegCS(:,c) = xRegCS(:,c) + removed(1,1,c);
        yRegCS(:,c) = yRegCS(:,c) + removed(2,1,c);
        
    end
end


%% Find coordinates to plot for each RGB image
    
% transform track coordinates to centre of cropped region
coord = trackCoord(plotChans,:,:);

% adjust to pixels
for i = 1:2
    coord(:,:,i) = coord(:,:,i)./repmat(pixelSize,nChans,1);
end
coordCS = coord*opts.subpixelate;
for c = 1:nChans

    if plotChans(c) ~= coordSysChan
        chrShift = job.options.chrShift.result{coordSysChan,plotChans(c)}(:,1:3)./pixelSize;
    else
        chrShift = zeros(1,3);
    end

    if opts.zoomTrack
        coord(c,1,:) = coord(c,1,:) - (xReg(1,plotChans(c))+1) + chrShift(1);
        coord(c,2,:) = coord(c,2,:) - (yReg(1,plotChans(c))+1) + chrShift(2);
        coordCS(c,1,:) = coordCS(c,1,:) - (xRegCS(1,plotChans(c))+1) + chrShift(1)*opts.subpixelate;
        coordCS(c,2,:) = coordCS(c,2,:) - (yRegCS(1,plotChans(c))+1) + chrShift(2)*opts.subpixelate;
    else
        coord(c,1,:) = coord(c,1,:) + chrShift(1);
        coord(c,2,:) = coord(c,2,:) + chrShift(2);
        coordCS(c,1,:) = coordCS(c,1,:) + chrShift(1)*opts.subpixelate;
        coordCS(c,2,:) = coordCS(c,2,:) + chrShift(2)*opts.subpixelate;
    end

end


%% Producing the figure

if nChans == 1
    plotImg = rgb2gray(rgbImg);
    plotImgCrpd = rgb2gray(rgbCrpd);
    plotCoord = coord;
else
    plotImg = rgbImgCS;
    plotImgCrpd = rgbCrpdCS;
    plotCoord = coordCS;
end 

% prepare figure environment
if opts.newFig
    figure
else
    figure(1)
end
clf

C = [ 0  1  0;
      1  0  0;
      0  0  1];
plotStyle = ['x' '+' 'o'];

if opts.zoomTrack
    
    % plot individual channels
    bigImgInd = [];
    if nChans > 1
        for c = 1:nChans
            subplot(nChans,nChans+1,c*(nChans+1))
            hold on
            imshow(plotImgCrpd(:,:,mapChans(c)))
            title(opts.plotTitles{c},'FontSize',20,'Color',C(plotChans(c),:))
            bigImgInd = [bigImgInd, (c-1)*(nChans+1)+1 : c*(nChans+1)-1 ];
        end
    end
    % plot dual image
    if isempty(bigImgInd)
        bigImgInd = 1;
    end
    subplot(nChans,nChans+1,bigImgInd)
    imshow(plotImgCrpd)
    title('All channels','FontSize',20)

    hold on
    % plot coordinates
    for c = 1:nChans
        for i = 1:2
            subplot(nChans,nChans+1,bigImgInd)
            if opts.transpose
                plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                    'Color',C(plotChans(c),:),'Marker',plotStyle(plotChans(c)),'MarkerSize',15)
            else
                plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                    'Color',C(plotChans(c),:),'Marker',plotStyle(plotChans(c)),'MarkerSize',15)
            end
            if nChans > 1
                subplot(nChans,nChans+1,c*(nChans+1))
                if opts.transpose
                    plot(plotCoord(c,1,i),plotCoord(c,2,i),...
                        'Color','k','Marker',plotStyle(plotChans(c)),'MarkerSize',15)
                else
                    plot(plotCoord(c,2,i),plotCoord(c,1,i),...
                        'Color','k','Marker',plotStyle(plotChans(c)),'MarkerSize',15)
                end
            end
        end
    end
    
else
    % plot image
    imshow(plotImg)
    
    hold on
    % plot tracked channel's coordinates (points too close together on
    % full image)
    for i = 1:2
        if opts.transpose
            plot(plotCoord(coordSysChan,1,i),plotCoord(coordSysChan,2,i),'wo','MarkerSize',15)
        else
            plot(plotCoord(coordSysChan,2,i),plotCoord(coordSysChan,1,i),'w','Marker',plotStyle(i),'MarkerSize',15)
        end
    end
    
end

hold off

end