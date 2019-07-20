function dublShowSisterPair(job,varargin)
% DUBLSHOWSISTERPAIR Plots image of dual-channel movie with coordinates
% for a given sister pair.
%
%    DUBLSHOWSISTERPAIR(JOB,...) Plots coordinates in two channels over
%    the movie image for a given sister pair at a given timepoint.
%
%    Options, defaults in {}:-
%
%    axis: {xy}, xz or yz. Coordinate axes in which to show images.
%
%    channelLabels: {'mNeonGreen','tagRFP','JF-646'} or similar. Labels of
%           each plot denoting which channel is which.
%
%    channelMap: {[2 1 3]} or some perturbation. Order in which the titled
%    channels are presented in the figure, where typically:
%    1=red, 2=green, 3=blue.
%
%    contrast: {{[0.1 1],[0.1 1],[0.1 1]}}, 'help', or similar. Upper and
%       lower contrast limits for each channel. Values must be in range
%       [0 1]. 'help' outputs the minimum and maximum intensity values as 
%       guidance, then requests values from the user.
%       Tips: - Increase the lower limit to remove background noise.
%             - Decrease the upper limit to increase brightness.
%
%    newFig: {0} or 1. Whether or not to show the sister pair in a new
%           figure.
%
%    channels: {[1 2]} or some subset of [1 2 3]. Vector of channels
%           for plotting, where typically:
%               1=red, 2=green, 3=blue.
%
%    project: 0, 1 or {-1}. Whether or not to project in the third axis.
%           Setting equal to -1 will project 1 µm distance surrounding the
%           sister pair.
%
%    sisterPair: {1} or number. Sister pair within JOB being plotted.
%
%    subpixelate: {9} or odd number. Number of pixels of accuracy within
%           which to correct for chromatic shift.
%
%    timePoint: {1} or number. Timepoint at which to plot the sister pair.
%
%    transpose: {0} or 1. Whether to tranpose the image.
%
%    zoomScale: {1} or number. Magnification of a zoomed image as a
%           proportion of the default, so that 0.5 will zoom out, and 1.5
%           will zoom in.
%
%    zoom: 0 or {1}. Whether or not to zoom into the specific sister
%           pair. A value of 0 plots the whole cell, with just the
%           coordinates of the sister pair plotted.
%
% Copyright (c) 2019 C. A. Smith

%% Process options

if nargin<1
  error('Must supply a job.');
end

% set default options
opts.axis = 'xy';
opts.channelLabels = {'mNeonGreen','tagRFP','JF-646'};
opts.channelMap = [2 1 3]; % green, red, blue
opts.channels = 1:2;
opts.contrast = repmat({[0.1 1]},1,3);
opts.newFig = 0;
opts.sisterPair = 1;
opts.subpixelate = 1;
opts.timePoint = 1;
opts.transpose = 0;
opts.zoomScale = 1;
opts.zoom = 1;
opts.project = -1;

% process options
opts = processOptions(opts, varargin{:});
while length(opts.channelLabels) < 3
    opts.channelLabels = [opts.channelLabels,{''}];
end

% convert axes into numbers
switch opts.axis
    case 'xy'
        axis = [1 2];
    case 'xz'
        axis = [1 3];
    case 'yz'
        axis = [2 3];
    otherwise
        error('Axes requested not recognised. Please provide either: ''xy'',''xz'' or ''yz''.');
end

%% Image and coordinate acquisition

% open movie
[md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata);

% get coordinate system and plot channels
coordSysChan = job.options.coordSystemChannel;
showChans = sort(opts.channels);
nChans = length(showChans);
plotChans = showChans;

% get sister information
sisPair = opts.sisterPair;

% get track information
trackIDs = job.dataStruct{coordSysChan}.sisterList(1).trackPairs(sisPair,1:2);
timePoint = opts.timePoint;

% if only one channel, run the single channel version
if nChans == 1
    kitLog('Only one channel being shown. Running kitShowSisterPair instead.');
    kitShowSisterPair(job,'channel',plotChans,'contrast',opts.contrast{plotChans},'newFig',opts.newFig,...
        'sisterPair',sisPair,'timePoint',timePoint,'title',opts.channelLabels{1},'transpose',opts.transpose,...
        'withinFig',0,'zoomScale',opts.zoomScale,'zoom',opts.zoom,'zProject',opts.project);
    return
end

% get pixel size
pixelSize = job.metadata.pixelSize(:,1:3);
chrShift = job.options.chrShift.result;

% accumulate track information by channel and sister
coord = nan(3,3,2);
for c = unique([plotChans coordSysChan])
    
    if ~isfield(job.dataStruct{c},'tracks')
    	plotChans = setdiff(plotChans,c);
        continue
    end 
    
    for iSis = 1:2
        
        tk = trackIDs(iSis);
        track = job.dataStruct{c}.tracks(tk);
        
        startTime = track.seqOfEvents(1,1);
        if timePoint < startTime
            coord(c,:,iSis) = nan(1,3);
        else
            coord(c,:,iSis) = ...
                track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
            coord(c,:,iSis) = coord(c,:,iSis) + chrShift{coordSysChan,c}(:,1:3);
            coord(c,:,iSis) = coord(c,:,iSis)./pixelSize;
        end
    end
    
end
coord = coord(:,[2 1 3],:);

% calculate pair centre
centrePxl = nanmean(coord(coordSysChan,:,:),3);
centrePxl = round(centrePxl);

mapChans = opts.channelMap;
if opts.transpose
    axis = fliplr(axis);
end

%% RGB image production

cropSize = job.ROI.cropSize;
rgbImg = zeros([cropSize(axis), 3]);

projAxis= setdiff(1:3,axis);

% produce raw image
for c = showChans
    
    % read stack
    img = kitReadImageStack(reader,md,timePoint,c,job.ROI.crop,0);
    
    % max projection
    if opts.project == 1
        img = max(img, [], projAxis);
    else
        switch projAxis
            case 1
                if opts.project == -1
                    projRange = ceil(0.5/pixelSize(projAxis));
                    img = max(img(centrePxl(1)-projRange:centrePxl(1)+projRange,:,:), [], 1);
                else
                    img = img(centrePxl(1),:,:);
                end
            case 2
                if opts.project == -1
                    projRange = ceil(0.5/pixelSize(projAxis));
                    img = max(img(:,centrePxl(2)-projRange:centrePxl(2)+projRange,:), [], 2);
                else
                    img = img(:,centrePxl(2),:);
                end
            case 3
                if opts.project == -1
                    projRange = ceil(0.5/pixelSize(projAxis));
                    img = max(img(:,:,centrePxl(3)-projRange:centrePxl(3)+projRange), [], 3);
                else
                    img = img(:,:,centrePxl(3));
                end
        end
    end
    img = squeeze(img);
    if opts.transpose
        img = img';
    end
    rgbImg(:,:,mapChans(c)) = img;
end

% define contrast stretch
for c = showChans
    if iscell(opts.contrast)
        irange = stretchlim(rgbImg(:,:,mapChans(c)),opts.contrast{c});
    elseif strcmp(opts.contrast,'help')
        % get maximum and minimum intensities of image
        intensityRange(1) = min(rgbImg(:,:,mapChans(c)));
        intensityRange(2) = max(rgbImg(:,:,mapChans(c)));
        % output possible range 
        fprintf('Channel %i. Range in third coordinate: [%i %i].\n',c,intensityRange(1),intensityRange(2));
        % request input from user
        userRange = input('Please provide range of intensities to image: ');
        while userRange(1)>userRange(2) 
            userRange = input('The maximum cannot be smaller than the minimum. Please try again: ');
        end
        irange = [userRange(1) userRange(2)];
        fprintf('\n');
    else
        irange = opts.contrast(c,:);
    end
    rgbImg(:,:,mapChans(c)) = imadjust(rgbImg(:,:,mapChans(c)), irange, []);
end

%% Produce chromatic shift corrected image

cropSize = job.ROI.cropSize*opts.subpixelate;
rgbImgCS = zeros([cropSize(axis), 3]);

first = 1;
for c = setdiff(showChans,coordSysChan)
  
  % switch axes on chrShift
  chrShift{coordSysChan,c} = chrShift{coordSysChan,c}(:,[2 1 3]);
    
  if first

    [img1,img2,new_cs] = ... 
        chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(c)),chrShift{coordSysChan,c}, ...
        'coords',axis,'pixelSize',pixelSize,'subpixelate',opts.subpixelate);

    % give the subpixelated images to the image structure
    rgbImgCS(:,:,mapChans(coordSysChan)) = img1;
    first = 0;

  else
    [~,img2,new_cs] = ...
        chrsComputeCorrectedImage(rgbImg(:,:,mapChans(coordSysChan)),rgbImg(:,:,mapChans(c)),chrShift{coordSysChan,c}, ...
        'coords',axis,'pixelSize',pixelSize,'subpixelate',opts.subpixelate);
  end
  
  % give the new subpixelated images to the image structure
  rgbImgCS(:,:,mapChans(c)) = img2;
  chrShift{coordSysChan,c}(:,axis) = (new_cs(axis,1,2)-new_cs(axis,2,2))';
  
end

%% Produce cropped image around pair centre

if opts.zoom
    % calculate spread about the centre pixel
    cropSpread = opts.zoomScale*(2./pixelSize);    
    cropSpread = ceil(cropSpread)*opts.subpixelate;
            
	% chromatic shifted regions
    centrePxl = centrePxl*opts.subpixelate;
    xReg = [centrePxl(axis(1))-cropSpread(axis(1)) centrePxl(axis(1))+cropSpread(axis(1))];
    yReg = [centrePxl(axis(2))-cropSpread(axis(2)) centrePxl(axis(2))+cropSpread(axis(2))];
    
    % get cropped region
    rgbCrpd = rgbImgCS(xReg(1):xReg(2),yReg(1):yReg(2),:);
    
end
            
%% Find coordinates to plot for each RGB image

% get just the axes we are interested in
coord = coord(:,axis,:);
coordCS = coord*opts.subpixelate;
for c = setdiff(plotChans,coordSysChan)
    for iSis = 1:2
        coordCS(c,:,iSis) = coordCS(c,:,iSis) - chrShift{coordSysChan,c}(:,axis);
    end
end

% correct coordinates for region position
if opts.zoom
    coordCS(:,1,:) = coordCS(:,1,:) - (xReg(1)-1) - (opts.subpixelate-1)/2;
    coordCS(:,2,:) = coordCS(:,2,:) - (yReg(1)-1) - (opts.subpixelate-1)/2;
end

%% Producing the figure

if opts.zoom
    plotImg = rgbCrpd;
else
    plotImg = rgbImgCS;
end
plotCoord = coordCS(:,[2 1],:);

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

if opts.zoom
    
    % plot individual channels
    bigImgInd = setdiff(1:nChans*(nChans+1),nChans+1:nChans+1:nChans*(nChans+1));
    for c = 1:nChans
        pc = showChans(c);
        subplot(nChans,nChans+1,c*(nChans+1))
        hold on
        imshow(plotImg(:,:,mapChans(pc)))
        title(opts.channelLabels{pc},'FontSize',20,'Color',C(pc,:))
    end
    subplot(nChans,nChans+1,bigImgInd)
    imshow(plotImg)
    title('All channels','FontSize',20)

    hold on
    % plot coordinates
    for c = 1:nChans
        pc = showChans(c);
        if ~ismember(plotChans,pc); continue; end
        for i = 1:2
            subplot(nChans,nChans+1,bigImgInd)
            plot(plotCoord(pc,1,i),plotCoord(pc,2,i),...
                'Color',C(pc,:),'Marker',plotStyle(pc),'MarkerSize',15)
            if nChans > 1
                subplot(nChans,nChans+1,c*(nChans+1))
                plot(plotCoord(pc,1,i),plotCoord(pc,2,i),...
                        'Color','k','Marker',plotStyle(pc),'MarkerSize',15)
            end
        end
    end
    
else
    % plot image
    imshow(plotImg)
    
    hold on
    % plot each channel's coordinates
    for c = 1:nChans
        pc = showChans(c);
        if ~ismember(plotChans,pc); continue; end
        for i = 1:2
            plot(plotCoord(pc,2,i),plotCoord(pc,1,i),...
                'Color',C(pc,:),'Marker',plotStyle(pc),'MarkerSize',15)
        end
    end
end

hold off

end