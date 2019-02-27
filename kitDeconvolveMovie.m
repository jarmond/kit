function deconMovie = kitDeconvolveMovie(job,movie,channel,verbose)
%
% Jonathan U Harrison 2019-02-19
%%%%%%%%%%%%%%%%%%%%%

if ~isfield(job,'dataStruct')
    error(['Assumes that PSF is defined in dataStruct.', ...
        'Use kitMakeMakiDatastruct to construct this']);
end

is3D = job.metadata.is3D;
ndims = 2 + is3D;
%define PSF
filters = createFilters(ndims,job.dataStruct{channel}.dataProperties);
PSF = fspecial3('gaussian',filters.signalP(4:6),filters.signalP(1:3));

deconMovie = zeros(size(movie));
nFrames = job.dataStruct{1}.dataProperties.movieSize(ndims+1);
for i=1:nFrames %loop is faster here than deconvolving full movie together
    deconMovie(:,:,:,i) = deconvlucy(movie(:,:,:,i),PSF,10); %10 iterations
end

%Take transpose in x,y. Min ensures normalized after deconvolution
%deconMovie = min(permute(deconMovie,[2,1,3,4]),1); 
deconMovie = min(deconMovie,1);

savename = sprintf('%s/%s',job.movieDirectory,job.ROI.movie);
[extStart,extEnd] = regexp(savename,'.ome.tif','start','end');
savename(extStart:extEnd)=[];
fprintf('Saving the deconvolved movie ... %s \n',savename);

%The following assumes bioformats is available
addpath bfmatlab;
bfCheckJavaPath(1);
bfsave(deconMovie, sprintf('%sDeconvolved.ome.tif',savename),...
    'dimensionOrder', 'XYZTC','BigTiff',true,'Compression', 'LZW');

if verbose
    figure;
    subplot(1,2,1);
    img = movie(:,:,:,10);
    imshow(max(img,[],3));
    subplot(1,2,2);
    deconImg = deconMovie(:,:,:,10);
    imshow(max(deconImg,[],3));
    figure;
    subplot(1,2,1);
    imshow(img(:,:,40));
    subplot(1,2,2);
    imshow(deconImg(:,:,40));
end
fprintf('Deconvolution finished\n');
end