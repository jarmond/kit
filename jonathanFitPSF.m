function [sigmaFitted,x] = jonathanFitPSF(psfMovie,showPlots)
% Here we will load in an experimental PSF and fit a gaussian to it
% Fit a full 3D gaussian(no assumptions about symmetry or diagonal entries)
%
%psfMovie - 3D array containing image of measured bead
%showPlots - logical to determine whether to show plots of image and model
%Jonathan U Harrison 2019-02-26
%%%%%%%%%%%%%%%%%%%%
rng('default');
if nargin < 1
    %load in psfMovie
    psf.movieDirectory = '../../Data/2019-02-26';
    psf.movie = 'OS_LLSM_181205_MC139_bead_images_for_PSF.ome.tif';
    %create reader to read movie
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [psf.metadata, reader] = kitOpenMovie(fullfile(psf.movieDirectory,...
        psf.movie));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %load movie
    psfMovie = kitReadWholeMovie(reader,psf.metadata,1,[],0,1);
end

if nargin < 2
    showPlots=1;
end

beadSize = 100; % in nm
pixelSize = [104,104,308];

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

fitFun = @(x) singleGaussianObjectiveFun(x,psfMovie,pixelList,...
    beadSize,pixelSize);
sigma0 = rand(1,3); %matlab filter functions need a diagonal matrix
%NB. to use non diagonal matrix for gaussian filter, could write something.
%But bead images fit well to a diagonal gaussian. (off diagonal entries 0
%to 2dp)
x0 = [0.0028,256,64,64,50,sigma0]; %Starting guess near middle of stack
if showPlots
    modelMovie = reshape(x0(2)*(sum((repmat(pixelSize,numel(psfMovie),1).* ...
        (bsxfun(@minus, pixelList, repmat(x0(3:5),numel(psfMovie),1)))).^2,2) < ...
        beadSize), size(psfMovie));
    modelMovie = imgaussfilt3(modelMovie,sigma0) + x0(1);
    %plot before
    plotComparisonObservedAndModel(psfMovie,modelMovie,2);
    title('Initial guess');
end
%%%%%%%%%%%%
%% perform optimization
[x,fp,flp,op] = fmincon(fitFun,x0,[],[],[],[],zeros(1,8),[1,1000,...
    size(psfMovie),10*ones(1,3)]);
sigmaFitted = diag(x(6:end));


if showPlots
    %%%%%%%%%%%%%
    %plot after
    modelMovie = reshape(x(2)*(sum((repmat(pixelSize,numel(psfMovie),1).* ...
        (bsxfun(@minus, pixelList, repmat(x(3:5),numel(psfMovie),1)))).^2,2) < ...
        beadSize), size(psfMovie));
    modelMovie = imgaussfilt3(modelMovie,x(6:end)) + x(1);
    plotComparisonObservedAndModel(psfMovie,modelMovie,2);
    title('Final guess');
    plotComparisonObservedAndModel(psfMovie,modelMovie,3);
    title('Final guess');
end

%%%%%%%%%%%%%%%%

end
function residual = singleGaussianObjectiveFun(theta,psfMovie,pixelList,...
                                                beadSize, pixelSize)
%% Suppose a model of the form b + a*exp(-(x-mu).^2/sigma^2)
b = theta(1); %background
a = theta(2); %amplitude
mu = theta(3:5); %mean: centre of bead
sigma = theta(6:end); %diagonal entries of psf

%Note: we rely on data about the size of the bead and the pixel size

%Take an indiator or binary mask for an object the same size as the bead is
%supposed to be. Put this at position mu (in pixel units).
%Multiply by some amplitude a. Make into the size of an image.
modelMovie = reshape(a*(sum((repmat(pixelSize,numel(psfMovie),1).* ...
    (bsxfun(@minus, pixelList, repmat(mu,numel(psfMovie),1)))).^2,2) < ...
    beadSize), size(psfMovie));
%Convolve the model image with a gaussian corresponding to the proposed PSF
%parameters sigma. Add a small amount of background b
modelMovie = imgaussfilt3(modelMovie,sigma) + b;
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