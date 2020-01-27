function [mu,sigma1] = jonathanEstimateSigma(Movie,verbose,writeToFile)
%%%%%%%%%%%%%%%%%%%%
rng('default');
if nargin < 1
    % TODO: add feature to open file system to select file
    error('No file selected. Input path to experimental PSF or an array');
    %     psf.movieDirectory = '../../Data/2019-02-26';
    %     psf.movie = 'OS_LLSM_181205_MC139_bead_images_for_PSF.ome.tif';
    %     psfMovie = fullfile(psf.movieDirectory,psf.movie)
end
if ischar(Movie)
    pathToMovie = Movie;
    %create reader to read movie
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [metadata, reader] = kitOpenMovie(pathToMovie);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %load movie
    Movie = kitReadWholeMovie(reader,metadata,1,[],0,1);
elseif isa(Movie,'double')
    pathToMovie=0; %no path provided, use default if writing to file
elseif ~isa(Movie,'double')
    error('Unknown type for psfMovie. Try inputting a string or array');
end

if nargin<2
    verbose=0;
end

if nargin<3
    writeToFile=0; %output to csv
end

%first draft: put in loop. TODO: vetorise this part (only done once though)
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

%we are first going to fit a single gaussian
fitFun1 = @(x) singleGaussianObjectiveFun(x,Movie,pixelList);
sigma0 = rand(3,3);
sigma0 = sigma0*sigma0'; %ensure positive definite
b0 = mode(Movie,'all'); %initial estimate of background from looking at image values
mu0 = size(Movie)/2;
x0 = [b0,256,mu0,sigma0(:)'];
%% perform optimization for initial model
fprintf('Initializing ...\n');
[x,~,~,~] = fmincon(fitFun1,x0,[],[],[],[],zeros(1,14),[1,256,...
    size(Movie),10*ones(1,9)])
fprintf('Done\n')
sigma1 = reshape(x(6:end),3,3);
sigma1 = sigma1*sigma1'

if writeToFile
   if ischar(pathToMovie) && (strcmp(pathToMovie(end-7:end),'.ome.tif'))
        fprintf('Saving estimated sigma...\n');
        csvwrite(sprintf('%s_sigma.csv',pathToMovie(1:end-8)),sigma1);
    else
        fprintf('Saving estimated sigma in default path...\n');
        csvwrite('~/sigma.csv',sigma1);
    end
end

if verbose
    figure; subplot(1,2,1);
    imshow(Movie(:,:,round(size(Movie,3)/2)),[],'InitialMagnification',1000)
    set(gca,'YDir','normal')
    subplot(1,2,2);
    mu = x(3:5); %// data
    x = 1:.1:size(Movie,1); %// x axis
    y = 1:.1:size(Movie,2); %// y axis
    
    [X Y] = meshgrid(x,y); %// all combinations of x, y
    Z = mvnpdf([X(:) Y(:), mu(3)*ones(size(X(:)))],mu,sigma1); %// compute Gaussian pdf
    Z = reshape(Z,size(X)); %// put into same size as X, Y
    contour(X,Y,Z), axis equal  %// contour plot; set same scale for x and y...    
end
fprintf('All done\n');

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
