function sigmaFitted = jonathanFitPSF(psfMovie,showPlots)
% Here we will load in an experimental PSF and fit a gaussian to it
% Fit a full 3D gaussian(no assumptions about symmetry or diagonal entries)
%
%psfMovie - 3D array containing image of measured bead
%showPlots - logical to determine whether to show plots of image and model
%Jonathan U Harrison 2019-02-26
%%%%%%%%%%%%%%%%%%%%
rng('default');
if nargin < 1
    % TODO: add feature to open file system to select file
    error('No file selected. Input path to experimental PSF or an array');
%     psf.movieDirectory = '../../Data/2019-02-26';
%     psf.movie = 'OS_LLSM_181205_MC139_bead_images_for_PSF.ome.tif';
%     psfMovie = fullfile(psf.movieDirectory,psf.movie)
end
if ischar(psfMovie)
%create reader to read movie
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [metadata, reader] = kitOpenMovie(psfMovie);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %load movie
    psfMovie = kitReadWholeMovie(reader,metadata,1,[],0,1);
elseif ~isa(psfMovie,'double')
    error('Unknown type for psfMovie. Try inputting a string or array');
end

if nargin < 2
    showPlots=1;
end

beadSize = 100; % in nm
pixelSize = [104,104,308*cos(57.2)];%100

%first draft: put in loop. TODO: vetorise this part (only done once though)
pixelList = zeros(numel(psfMovie),3);
ind = 0;
for k=1:size(psfMovie,3)
    for j=1:size(psfMovie,2)
        for i=1:size(psfMovie,1)
            ind = ind+1;
            pixelList(ind,:) = [i,j,k];
        end
    end
end

%we are first going to fit a single gaussian assuming sub pixel point for
%the bead. We will use this to inform an initial condition for a model
%where the bead has finite size not necessarily restricted to a single
%pixel

fitFun1 = @(x) singleGaussianObjectiveFun(x,psfMovie,pixelList);
sigma0 = rand(3,3);
sigma0 = sigma0*sigma0'; %ensure positive definite
b0 = 0.0028; %initial estimate of background from looking at image values
x0 = [b0,256,size(psfMovie)/2,sigma0(:)'];
%% perform optimization for initial model
fprintf('Initializing ...\n');
[x,~,~,~] = fmincon(fitFun1,x0,[],[],[],[],zeros(1,14),[1,256,...
    size(psfMovie),10*ones(1,9)]);
fprintf('Done\n')
sigma1 = reshape(x(6:end),3,3);
sigma1 = sigma1*sigma1';
sigma1 = diag(sigma1)'/4; %matlab filter functions need a diagonal matrix
%NB. to use non diagonal matrix for gaussian filter, could write something.
%But bead images fit well to a diagonal gaussian. (off diagonal entries 0
%to 2dp)
b1 = x(1); 
mu1 = x(3:5); %Starting guess based on previous result
x1 = [b1,sum(psfMovie(:)-b1),mu1,sigma1]; 
fitFun2 = @(x) GaussianFilterObjectiveFun(x,psfMovie,pixelList,...
    beadSize,pixelSize);
fltXYZ = roundOddOrEven(4*sigma1,'odd','inf');
if showPlots
    modelMovie = reshape(x1(2)*(sqrt(sum((repmat(pixelSize, ...
        numel(psfMovie),1).*(bsxfun(@minus, pixelList, repmat(x1(3:5), ...
        numel(psfMovie),1)))).^2,2)) < beadSize), size(psfMovie));
    modelMovie = imgaussfilt3(modelMovie,sigma1,'FilterSize', ...
        max(fltXYZ)*ones(1,3)) + x1(1);
    %plot before
    plotComparisonObservedAndModel(psfMovie,modelMovie,2);
    title('Initial guess');
end
%%%%%%%%%%%%
%% perform optimization
fprintf('Solving for PSF...\n');
[x,~,~,~] = fmincon(fitFun2,x1,[],[],[],[],zeros(1,8),[1,256,...
    size(psfMovie),20*ones(1,3)]);

sigmaFitted = x(6:end);
%this has been in units of pixels
%account for pixel size and convert to um for consistency with rest of KiT
sigmaFitted = sigmaFitted.*pixelSize/10^3; 
fprintf('Done: PSF is x:%f y:%f z:%f in units of um \n',sigmaFitted);

if showPlots
    %%%%%%%%%%%%%
    %plot after
    fltXYZ = roundOddOrEven(4*x(6:end),'odd','inf');
    modelMovie = reshape(x(2)*(sqrt(sum((repmat(pixelSize,...
        numel(psfMovie),1).*(bsxfun(@minus, pixelList, repmat(x(3:5), ...
        numel(psfMovie),1)))).^2,2)) < beadSize), size(psfMovie));
    modelMovie = imgaussfilt3(modelMovie,x(6:end),'FilterSize', ...
        max(fltXYZ)*ones(1,3)) + x(1);
    plotComparisonObservedAndModel(psfMovie,modelMovie,2);
    title('Final guess');
    plotComparisonObservedAndModel(psfMovie,modelMovie,3);
    title('Final guess');
end

%%%%%%%%%%%%%%%%

end
function residual = singleGaussianObjectiveFun(theta,psfMovie,pixelList)
%% Suppose a model of the form b + a*exp(-(x-mu).^2/sigma^2)
b = theta(1); %background
a = theta(2); %amplitude
mu = theta(3:5); %mean: centre of gaussian
sigma = reshape(theta(6:end),3,3); %covariance: psf
sigma = sigma*sigma'; %ensure symmetric and positive definite

modelMovie = reshape(b + a*mvnpdf(pixelList,mu,sigma),size(psfMovie));
residual = sum(sum(sum(abs(bsxfun(@minus,psfMovie, modelMovie)))));
end
function residual = GaussianFilterObjectiveFun(theta,psfMovie,pixelList,...
                                                beadSize, pixelSize)
%% Make a mask for the bead image and convolve with gaussian filter
b = theta(1); %background
a = theta(2); %amplitude
mu = theta(3:5); %mean: centre of bead
sigma = theta(6:8); %diagonal entries of psf

%Note: we rely on data about the size of the bead and the pixel size

%Take an indiator or binary mask for an object the same size as the bead is
%supposed to be. Put this at position mu (in pixel units).
%Multiply by some amplitude a. Make into the size of an image.
fltXYZ = roundOddOrEven(4*theta(1:3),'odd','inf');
modelMovie = reshape(a*(sqrt(sum((repmat(pixelSize,numel(psfMovie),1).* ...
    (bsxfun(@minus, pixelList, repmat(mu,numel(psfMovie),1)))).^2,2)) < ...
    beadSize), size(psfMovie));
%Convolve the model image with a gaussian corresponding to the proposed PSF
%parameters sigma. Add a small amount of background b
modelMovie = imgaussfilt3(modelMovie,sigma,'FilterSize',max(fltXYZ)*...
    ones(1,3)) + b;
%Compute a total L1 loss across all voxels. We will aim to minimise this.
residual = sum(sum(sum(abs(bsxfun(@minus,psfMovie, modelMovie)))));
end
function plotComparisonObservedAndModel(psfMovie,modelMovie,dim)
if nargin<3
    dim=3; %project in z. May want to project in other dimensions too
end
figure;
subplot(1,2,1);
imshow(squeeze(max(psfMovie,[],dim)));
subplot(1,2,2);
imshow(squeeze(max(modelMovie,[],dim)));
end