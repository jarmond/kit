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

fitFun = @(x) singleGaussianObjectiveFun(x,psfMovie,pixelList);
sigma0 = rand(3,3);
sigma0 = sigma0*sigma0';
x0 = [0.0028,256,64,64,50,sigma0(:)']; %Starting guess near middle of stack
if showPlots
modelMovie = reshape(x0(1) + x0(2)*mvnpdf(pixelList,x0(3:5),...
    reshape(x0(6:end),3,3)),size(psfMovie));
%plot before
plotComparisonObservedAndModel(psfMovie,modelMovie,2);
title('Initial guess');
end
%%%%%%%%%%%%
%% perform optimization
[x,fp,flp,op] = fmincon(fitFun,x0,[],[],[],[],zeros(1,14),[1,256,...
    size(psfMovie),10*ones(1,9)]);
sigmaFitted = reshape(x(6:end),3,3);
sigmaFitted = sigmaFitted*sigmaFitted';
modelMovie = reshape(x(1) + x(2)*mvnpdf(pixelList,x(3:5),...
    sigmaFitted),size(psfMovie));
%%%%%%%%%%%%%
if showPlots
%plot after
plotComparisonObservedAndModel(psfMovie,modelMovie,2);
title('Final guess');
plotComparisonObservedAndModel(psfMovie,modelMovie,3);
title('Final guess');
end

%%%%%%%%%%%%%
%What we have determined thus far corresponds to the true bead convolved
%with the psf that we want to know about.
%We have to account for this.

beadSize = 100*ones(1,3); % in nm, each dimension
pixelSize = [104,104,308];
beadSizeInPixels = beadSize./pixelSize;

%define Bead to deconvolve with to leave PSF?
patchXYZ = roundOddOrEven(4*beadSizeInPixels,'odd','inf');
beadFilter = fspecial3('gaussian',patchXYZ,beadSizeInPixels);
deconBead = deconvlucy(psfMovie,beadFilter,10); %10 iterations
deconBead = min(deconBead,1);

plotComparisonObservedAndModel(psfMovie,deconBead,2);


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
function plotComparisonObservedAndModel(psfMovie,modelMovie,dim)
if nargin<3
    dim=3; %project in z. May want to project in other dimensions too
end
figure;
subplot(1,2,1)
imshow(squeeze(max(psfMovie,[],dim)));
subplot(1,2,2);
imshow(squeeze(max(modelMovie,[],dim)));
end