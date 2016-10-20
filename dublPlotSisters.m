function dublPlotSisters(job, varargin)
%DUBLPLOTSISTERS Produces a set of plots as in kitBasicPlots but for both
%              channels.
%
%    DUBLPLOTSISTERS(JOB,...) Plots all sister pairs for two channels
%    overlaid. By default also plots intersister distances and intermarker
%    distances. Supply options as string/value pairs following JOB.
%
%    Options, defaults in {}:-
%
%    channels: {1:2} or two numbers. Channels for which tracks are plotted.
%
%    col: {green; red; blue; light grey} or [4x3]-matrix of RGB values.
%         Colours to plot first channel, second channel, delta values, and
%         overlaid plots, respectively.
%
%    filter: {0} or 1. Filter plots by either z-depth, or z-angle < 30deg.
%
%    plotCoord: {1}, 2 or 3. Coordinate for which trajectories are plotted:
%               1=x, 2=y, 3=z.
%
%    plots: {1} or some subset of the [1:3]. Designates which data to
%           plot:
%               1 = sister pair trajectories
%               2 = interkinetochore distances
%               3 = intrakinetochore distances (including their average)
%
%    subset: {all sisters} or some vector of sister numbers. Subset of
%            sisters for plotting.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2014 C. A. Smith

if nargin<1 || isempty(job)
    error('Need to supply a job.')
end

% set up default options
opts.channels = 1:2;
opts.col = [ 0 ,0.5, 0  ;
             1 , 0 , 0  ;
             0 , 0 , 1 ];
opts.filter = 1;
opts.filterAngle = 30;
opts.filterType = 'zDepth';
opts.plotCoord = 1;
opts.plots = 1;
opts.subset = [];
opts.dualOnly = 1;

% and process the input options
opts = processOptions(opts, varargin{:});

% get all job information
dS1 = job.dataStruct{opts.channels(1)};
dS2 = job.dataStruct{opts.channels(2)};
sisterList1 = dS1.sisterList;
sisterList2 = dS2.sisterList;
nSisters1 = length(sisterList1);
nSisters2 = length(sisterList2);
if nSisters1 ~= nSisters2
    error('Do not have equal number of sisters. Cannot overlay plots.')
else
    nSisters = nSisters1;
    clear nSisters1 nSisters2
end
if isempty(opts.subset)
    opts.subset = 1:nSisters;
elseif max(opts.subset) > nSisters
    opts.subset = 1:nSisters;
    warning('Supplied sister subset out of range for cell. Have taken full set.')
end

% convert filterType to 'none' if no filter required
if ~opts.filter
    opts.filterType = 'none';
end

pixelSize = job.metadata.pixelSize(1:3);
filterDepth = pixelSize(3)/2;

switch opts.plotCoord
    case 1
        yLabPlots = 'x position';
    case 2
        yLabPlots = 'y position';
    case 3
        yLabPlots = 'z position';
    otherwise
        yLabPlots = '';
end

% find which sisters have both channels
if opts.dualOnly
    newSubset = [];
    for iSis = opts.subset
        
        % calculate for individual scenarios
        failed1 = sum(~isnan(sisterList1(iSis).coords1(:,opts.plotCoord)))==0;
        failed2 = sum(~isnan(sisterList1(iSis).coords2(:,opts.plotCoord)))==0;
        failed3 = sum(~isnan(sisterList2(iSis).coords1(:,opts.plotCoord)))==0;
        failed4 = sum(~isnan(sisterList2(iSis).coords2(:,opts.plotCoord)))==0;
        
        % combine for ultimate failure
        failed = (failed1*failed2 + failed3*failed4)~=0;
        % if not failed, add sister to the list
        if ~failed
            newSubset = [newSubset,iSis];
        end
        
    end
    opts.subset = newSubset;
end

% Plot sister tracks in X.
if sum(opts.plots == 1) >= 1
    figure();
    t = job.metadata.frameTime(1,:);
    n = length(opts.subset);
    fig_n=ceil(sqrt(n));
    fig_m=ceil(n/fig_n);
    clf;
    for i=1:length(opts.subset)
        subplot(fig_m,fig_n,i);
        
        pair1 = sisterList1(opts.subset(i));
        pair2 = sisterList2(opts.subset(i));
        
        hold on
        
        switch opts.filterType
        
            case 'zContortion'
                % get planeFit for rotation
                planeFit = dS1.planeFit;

                % get coordinates in microscope coordinates
                rotCoordsChan1sis1 = changeCoordinateSystem(pair1.coords1(:,1:3),planeFit,'p2m');
                rotCoordsChan2sis1 = changeCoordinateSystem(pair2.coords1(:,1:3),planeFit,'p2m');
                rotCoordsChan1sis2 = changeCoordinateSystem(pair1.coords2(:,1:3),planeFit,'p2m');
                rotCoordsChan2sis2 = changeCoordinateSystem(pair2.coords2(:,1:3),planeFit,'p2m');

                % calculate z-contortion
                zRotDiff(1,:) = rotCoordsChan2sis1(:,3)-rotCoordsChan1sis1(:,3);
                xyRotDiff(1,:) = sqrt(sum(rotCoordsChan2sis1(:,1:2)-rotCoordsChan1sis1(:,1:2),2));
                zContor(1,:) = atan(zRotDiff(1,:)./xyRotDiff(1,:))*180/pi;

                zRotDiff(2,:) = rotCoordsChan2sis2(:,3)-rotCoordsChan1sis2(:,3);
                xyRotDiff(2,:) = sqrt(sum(rotCoordsChan2sis2(:,1:2)-rotCoordsChan1sis2(:,1:2),2));
                zContor(2,:) = atan(zRotDiff(2,:)./xyRotDiff(2,:))*180/pi;

                % find points being kept after filtering
                filtered(:,1) = +(abs(zContor(1,:))<opts.filterAngle);
                filtered(:,2) = +(abs(zContor(2,:))<opts.filterAngle);
                filtered(filtered==0) = NaN;
                
            case 'zDepth'

                % get planeFit for rotation
                planeFit = dS1.planeFit;

                % get coordinates in microscope coordinates
                rotCoordsChan1sis1 = changeCoordinateSystem(pair1.coords1(:,1:3),planeFit,'p2m');
                rotCoordsChan2sis1 = changeCoordinateSystem(pair2.coords1(:,1:3),planeFit,'p2m');
                rotCoordsChan1sis2 = changeCoordinateSystem(pair1.coords2(:,1:3),planeFit,'p2m');
                rotCoordsChan2sis2 = changeCoordinateSystem(pair2.coords2(:,1:3),planeFit,'p2m');

                % calculate z-contortion for filtering if required
                zRotDiff(1,:) = rotCoordsChan2sis1(:,3)-rotCoordsChan1sis1(:,3);

                zRotDiff(2,:) = rotCoordsChan2sis2(:,3)-rotCoordsChan1sis2(:,3);

                % find points being kept after filtering
                filtered(:,1) = +(abs(zRotDiff(1,:))<filterDepth);
                filtered(:,2) = +(abs(zRotDiff(2,:))<filterDepth);
                filtered(filtered==0) = NaN;

            otherwise
                
                % if no filtering required, produce vector of 1s
                filtered = ones(job.metadata.nFrames,2);

        end
        
        plot(t, pair1.coords1(:,opts.plotCoord).*filtered(:,1),'Color',opts.col(1,:));
        plot(t, pair1.coords2(:,opts.plotCoord).*filtered(:,2),'Color',opts.col(1,:));
        plot(t, pair2.coords1(:,opts.plotCoord).*filtered(:,1),'Color',opts.col(2,:));
        plot(t, pair2.coords2(:,opts.plotCoord).*filtered(:,2),'Color',opts.col(2,:));
        sisTit = sprintf('%i',opts.subset(i));
        title(sisTit)
        if mod(i-1,fig_n)==0
            ylabel(yLabPlots);
        end
        if i>(fig_m-1)*fig_n
            xlabel('Time (secs)');
        end
        xlim([min(t) max(t)]);
    end
end

% Plot sister separations in X.
if sum(opts.plots == 2) >= 1
    figure();
    clf;
    for i=1:length(opts.subset)
        subplot(fig_m,fig_n,i);
        sep1 = sisterList1(opts.subset(i)).distances;
        sep2 = sisterList2(opts.subset(i)).distances;
        hold on
        plot(t, sep1(:,opts.plotCoord),'Color',opts.col(1,:));
        plot(t, sep2(:,opts.plotCoord),'Color',opts.col(2,:));
        if mod(i-1,fig_n)==0
            ylabel('sisterSep');
        end
        if i>(fig_m-1)*fig_n
            xlabel('time, s');
        end
        xlim([min(t) max(t)]);
    end
end

% Plot marker separations in X.
if sum(opts.plots == 3) >= 1
    figure();
    clf;
    for i=1:length(opts.subset)
        subplot(fig_m,fig_n,i);
        pair1 = sisterList1(opts.subset(i));
        pair2 = sisterList2(opts.subset(i));

        avgDel = abs(sisterList2(opts.subset(i)).distances - sisterList1(opts.subset(i)).distances)/2;
        delta1 = abs(pair1.coords1(:,opts.plotCoord) - pair2.coords1(:,opts.plotCoord));
        delta2 = abs(pair1.coords2(:,opts.plotCoord) - pair2.coords2(:,opts.plotCoord));
        hold on
        plot(t, avgDel(:,1),'Color',opts.col(3,:));
        plot(t, delta1(:,1),'Color',opts.col(1,:));
        plot(t, delta2(:,1),'Color',opts.col(2,:));
        if mod(i-1,fig_n)==0
            ylabel('delta');
        end
        if i>(fig_m-1)*fig_n
            xlabel('time, s');
        end
        xlim([min(t) max(t)]);
    end
end


function rotCoords = changeCoordinateSystem(coords,planeFit,direction)

% preallocation of variables
nFrames = size(coords,1);
framesWiPlane = []; framesNoPlane = [];
rotCoords = coords;

% find frames with and without plane fits
for iFrame = 1:nFrames;
    if isempty(planeFit(iFrame).plane);
        framesNoPlane = [framesNoPlane,iFrame];
    else
        framesWiPlane = [framesWiPlane,iFrame];
    end
end

%first find coordinate system for each timepoint
if length(framesWiPlane) > 1
    %get first frame to rotate
    firstFrameRotate = framesWiPlane(1);

    %find frames without plane that are before the first frame with plane
    framesBefore1 = framesNoPlane(framesNoPlane < firstFrameRotate);
    framesAfter1 = setxor(framesNoPlane,framesBefore1); % the rest of the frames

    %get the coordinate system of each frame with a plane
    coordSystem = zeros(3,3,nFrames);
    coordSystem(:,:,framesWiPlane) = cat(3,planeFit.planeVectors);

    %assign the coordinate system of each frame without a plane
    %if they are before the first frame with a plane, take the coordinate
    %system of the plane after
    %if they are after the first frame with a plane, take the coordinate
    %system of the plane before
    for iTime = framesBefore1(end:-1:1)
        coordSystem(:,:,iTime) = coordSystem(:,:,iTime+1);
    end
    for iTime = framesAfter1
        coordSystem(:,:,iTime) = coordSystem(:,:,iTime-1);
    end

end
if strcmp(direction,'p2m')
    for iImage = framesWiPlane
        %rotate the vectors
        rotCoords(iImage,:) = (coordSystem(:,:,iImage)*rotCoords(iImage,:)')';
    end
elseif strcmp(direction,'m2p')
    for iImage = framesWiPlane
        %rotate the vectors
        rotCoords(iImage,:) = (coordSystem(:,:,iImage)/rotCoords(iImage,:)')';
    end
end


