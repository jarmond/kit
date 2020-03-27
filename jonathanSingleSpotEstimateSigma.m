function sigma1 = jonathanSingleSpotEstimateSigma(job,md,reader,all_coords,amplitude,varargin)
%
%Compare also JonathanEstimateSigma and kitShowSisterPair
%
%all_coords contains positions of all detected spots in given frame so that
%a theoretical image can be computated so that we can effectively remove
%spots that we know about when fitting the shape of a single spot
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
opts.sisID = []; %1 or 2

% process options
opts = processOptions(opts, varargin{:});

%% Image and coordinate acquisition

if nargin < 3
    % open movie
    [md,reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),job.metadata,0);
end

% get coordinate system and plot channels
chan = opts.channel;

% get how close in to zoom to
zoomRangeMicrons = opts.zoomRangeMicrons;

% get sister information
sisPair = opts.sisterPair;
chosenSisInd = opts.sisID;
sisterList = job.dataStruct{chan}.sisterList;

if sisPair > length(sisterList)
    kitLog('Sister pair ID provided too large. Only have %i sisters in movie. Quitting.',length(sisterList))
    return
end

% get track information
timePoint = opts.timePoint;
% get pixel resolution
pixelSize = job.metadata.pixelSize;

% refChan = job.options.coordSystemChannel;
% % get chromatic shift
% chrShift = job.options.chrShift.result{refChan,chan}(1:3);

%trackIDs = sisterList(1).trackPairs(sisPair,1:2);
% coords = nan(2,3); %coords x sister
%
% % accumulate track information by channel and sister
% for iSis = 1:2
%     tk = trackIDs(iSis);
%     track = job.dataStruct{chan}.tracks(tk);
%
%     startTime = track.seqOfEvents(1,1);
%     endTime   = track.seqOfEvents(2,1);
%     if timePoint < startTime || timePoint > endTime
%         coords(iSis,:) = nan(1,3);
%     else
%         coords(iSis,:) = ...
%             track.tracksCoordAmpCG(8*(timePoint-(startTime-1))-7:8*(timePoint-(startTime-1))-5);
%         coords(iSis,:) = coords(iSis,:) + chrShift;
%         coords(iSis,:) = coords(iSis,:)./pixelSize;
%     end
% end

if sum(isnan(all_coords(timePoint,:,:,sisPair)),'all') == 6
    kitLog('No coordinates found for sister pair %i at time point %i. Ignoring.',sisPair,timePoint);
    sigma1 = nan(3,3);
    return
end

% calculate pair centre and convert to pixels
% centrePxl = nanmean(coords);
% centrePxl = round(centrePxl);

%default uses sister 1 provided this is not missing, otherwise use what is
%specified
if isempty(chosenSisInd) && (~any(isnan(all_coords(timePoint,:,1,sisPair))))
    %choose sister 1 (check if is nan)
    chosenSisInd = 1;
elseif isempty(chosenSisInd) && (~any(isnan(all_coords(timePoint,:,2,sisPair))))
    %already checked if both are nan so sister 2 should not be nan
    chosenSisInd = 2;
elseif (~isempty(chosenSisInd)) && (~any(isnan(all_coords(timePoint,:,chosenSisInd,sisPair)))) 
%stick with chosen Sis Ind
else %chosen sister is missing in that frame
    kitLog('No coordinates found for sister pair %i at time point %i. Ignoring.',sisPair,timePoint);
    sigma1 = nan(3,3);
    return
end
centrePxl = round(all_coords(timePoint,:,chosenSisInd,sisPair));
%% Image production

% read stack
img = kitReadImageStack(reader,md,timePoint,chan,job.ROI.crop,0);
xReg = [centrePxl(1)-ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(1)))+1 ...
    centrePxl(1)+ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(1)))+1];
yReg = [centrePxl(2)-ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(2)))+1 ...
    centrePxl(2)+ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(2)))+1];
zReg = [centrePxl(3)-ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(3)))+1 ...
    centrePxl(3)+ceil(opts.zoomScale*(zoomRangeMicrons/pixelSize(3)))+1];
xReg(1) = max(1,xReg(1));
yReg(1) = max(1,yReg(1));
zReg(1) = max(1,zReg(1));
xReg(2) = min(xReg(2),job.metadata.frameSize(1));
yReg(2) = min(yReg(2),job.metadata.frameSize(2));
zReg(2) = min(zReg(2),job.metadata.frameSize(3));
if opts.transpose
    imgCrpd = img(xReg(1):xReg(2),yReg(1):yReg(2),zReg(1):zReg(2));
else
    xReg(2) = min(xReg(2),job.metadata.frameSize(2));
    yReg(2) = min(yReg(2),job.metadata.frameSize(1));
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

modelMovie = zeros(numel(imgCrpd),1);
psfSigma = job.dataStruct{chan}.dataProperties.FILTERPRM(1:3);
nPairs = size(sisterList(1).trackPairs,1);
for j = setdiff(1:nPairs,sisPair) %exclude the sister pair of interest
    for i = 1:2
        xDists = repmat(ones(yReg(2)-yReg(1)+1,1) * (xReg(1):xReg(2)),1,1,zReg(2)-zReg(1)+1);
        yDists = repmat((yReg(1):yReg(2))' * ones(1,xReg(2)-xReg(1)+1),1,1,zReg(2)-zReg(1)+1);
        zDists = zeros(size(xDists));
        for k=1:(zReg(2)-zReg(1)+1)
            zDists(:,:,k) = ones(yReg(2)-yReg(1)+1,xReg(2)-xReg(1)+1)*(zReg(1)-1+k);
        end
        pixelList = cat(2,xDists(:),yDists(:),zDists(:));
        if ~any(isnan(all_coords(timePoint,:,i,j)))
            modelMovie = modelMovie + amplitude(timePoint,:,i,j)/amplitude(timePoint,:,chosenSisInd,sisPair)*prod(normpdf(pixelList,all_coords(timePoint,:,i,j),psfSigma),2); %may want to think about using logs, but I think for this purpose overflow not important
        end
    end
end
%%also try sister of selected KT
if ~any(isnan(all_coords(timePoint,:,3-chosenSisInd,sisPair)))
    modelMovie = modelMovie + amplitude(timePoint,:,chosenSisInd,sisPair)/amplitude(timePoint,:,3-chosenSisInd,sisPair)*prod(normpdf(pixelList,all_coords(timePoint,:,3-chosenSisInd,sisPair),psfSigma),2);
end
if opts.verbose
%    psfSigma
%    reshape(modelMovie,size(imgCrpd))
    max(modelMovie,[],'all')
    max(imgCrpd,[],'all')
    figure; imshow(max(reshape(modelMovie,size(imgCrpd)),[],3),[],'InitialMagnification',1000);
end
scaling = max(imgCrpd,[],'all')/max(modelMovie,[],'all')/1.5;
imgCrpd = bsxfun(@minus,imgCrpd, scaling*reshape(modelMovie,size(imgCrpd)));
%[mu,sigma1] = jonathanEstimateSigma(imgCrpd,opts.verbose,0,opts.sigmaScale);

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

%%%%%%%%%%%%%%%%%%%%%
Movie = imgCrpd;
verbose = opts.verbose;
scaleSigma=opts.sigmaScale; %for spots, helpful to be able to scale the covariance

pixelList = zeros(numel(Movie),3);
ind = 0;
for k=1:size(Movie,3)
    for j=1:size(Movie,2)
        for i=1:size(Movie,1)
            ind = ind+1;
            pixelList(ind,:) = [i,j,k];
        end
    end
end

% fit a single gaussian
fitFun1 = @(x) singleGaussianObjectiveFun(x,Movie,pixelList,scaleSigma);
if scaleSigma == 1
    sigma0 = rand(3,3);
    sigma0 = sigma0*sigma0'; %ensure positive definite
else
    sigma0 = scaleSigma*ones(3,1);
end
movie_sz = size(Movie);
if size(Movie,3)==1
    movie_sz = [movie_sz, 1];
end
b0 = mode(Movie,'all'); %initial estimate of background from looking at image values
mu0 = movie_sz/2;
x0 = [b0,1+255*(scaleSigma==1),mu0,sigma0(:)'];
if isnan(singleGaussianObjectiveFun(x0,Movie,pixelList,scaleSigma))
fprintf('oops it is nan');
x0
Movie
pixelList
scaleSigma
scaling
max(modelMovie,[],'all')
end
%% perform optimization for initial model
%fprintf('Initializing ...\n');
if verbose
    options = optimoptions('fmincon','Display', 'notify');
else
    options = optimoptions('fmincon','Display', 'off');
end
num_unknowns = 8 + 6*(scaleSigma==1); %default to learning full covariance and initializing around scale 1ish
[x,~,~,~] = fmincon(fitFun1,x0,[],[],[],[],zeros(1,num_unknowns),[1,256,...
    movie_sz,10*ones(1,num_unknowns-5)],[],options);
%fprintf('Done\n')
if scaleSigma == 1
    sigma1 = reshape(x(6:end),3,3);
    sigma1 = sigma1*sigma1';
else
    sigma1 = diag(x(6:end));
end
mu = x(3:5);

if verbose
    figure; imshow(imgCrpd(:,:,round(size(imgCrpd,3)/2)),[],'InitialMagnification',5000)
    figure; subplot(1,2,1);
    imshow(Movie(:,:,round(size(Movie,3)/2)),[],'InitialMagnification',1000)
    set(gca,'YDir','normal')
    subplot(1,2,2);
    xx = 1:.1:size(Movie,1); %// x axis
    yy = 1:.1:size(Movie,2); %// y axis
    
    [X, Y] = meshgrid(xx,yy); %// all combinations of x, y
    Z = mvnpdf([X(:) Y(:), round(size(Movie,3)/2)*ones(size(X(:)))],mu,sigma1); %// compute Gaussian pdf
    Z = reshape(Z,size(X)); %// put into same size as X, Y
    contour(X,Y,Z), axis equal  %// contour plot; set same scale for x and y...
end
%fprintf('All done\n');

end

function residual = singleGaussianObjectiveFun(theta,psfMovie,pixelList,scaleSigma)
%% Suppose a model of the form b + a*exp(-(x-mu).^2/2*sigma^2)
b = theta(1); %background
a = theta(2); %amplitude
mu = theta(3:5); %mean: centre of gaussian

if scaleSigma==1
    sigma = reshape(theta(6:end),3,3); %covariance: psf
    sigma = sigma*sigma'; %ensure symmetric and positive definite
else
    sigma = diag(theta(6:end));
end

modelMovie = reshape(b + a*mvnpdf(pixelList,mu,sigma),size(psfMovie));
residual = sum(sum(sum(abs(bsxfun(@minus,psfMovie, modelMovie)))));
end
