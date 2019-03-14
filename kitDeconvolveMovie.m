function deconMovie = kitDeconvolveMovie(job,movie,channel,...
                                         pathToExpmtPSF,saveoutput,...
                                         verbose)
%kitDeconvolveMovie should perform deconvolution using the Richardson-Lucy
%algorithm to reduce the blur on raw images. We can use either the
%theoretical of experimentally determined PSF for this. Better results are
%seen with experimental PSFs. 
%job - struct containing info about job
%movie - 4D double array
%channel - int, usually 1
%pathToExpmtPSF - char, path to PSF image,
%                                     
% Jonathan U Harrison 2019-02-19
%%%%%%%%%%%%%%%%%%%%%

if ~isfield(job,'dataStruct')
    error(['Assumes that PSF is defined in dataStruct.', ...
        'Use kitMakeMakiDatastruct to construct this']);
end

is3D = job.metadata.is3D;
ndims = 2 + is3D;

if ischar(movie) || isempty(movie)
        %user has provided a path to movie file
        if ischar(movie)
            pathToMovie = movie; 
        elseif isempty(movie)
            pathToMovie = fullfile(job.movieDirectory,job.ROI.movie);
        end
        [metadata,reader] = kitOpenMovie(pathToMovie);
        movie = kitReadWholeMovie(reader,metadata,channel,[],0,1);
end

if ischar(pathToExpmtPSF)
%    PSF = pathToExpmtPSF;
    PSF = jonathanFitPSF(pathToExpmtPSF,0);    
    fltXYZ = roundOddOrEven(4*PSF,'odd','inf');
    PSF = fspecial3('gaussian',max(fltXYZ)*ones(1,3),PSF);
else
    %define PSF based on theoretical properties
    filters = createFilters(ndims,job.dataStruct{channel}.dataProperties);
    PSF = fspecial3('gaussian',filters.signalP(4:6),filters.signalP(1:3));
end

deconMovie = zeros(size(movie));
nFrames = job.dataStruct{channel}.dataProperties.movieSize(ndims+1);
for i=1:nFrames %loop is faster here than deconvolving full movie together
    deconMovie(:,:,:,i) = deconvlucy(movie(:,:,:,i),PSF,10); %10 iterations
end

%Take transpose in x,y. Min ensures normalized after deconvolution
%deconMovie = min(permute(deconMovie,[2,1,3,4]),1); 
deconMovie = min(deconMovie,1);

if saveoutput
savename = sprintf('%s/%s',job.movieDirectory,job.ROI.movie);
[extStart,extEnd] = regexp(savename,'.ome.tif','start','end');
savename(extStart:extEnd)=[];
fprintf('Saving the deconvolved movie ... %s \n',savename);

%The following assumes bioformats is available
addpath bfmatlab;
bfCheckJavaPath(1);
bfsave(deconMovie, sprintf('%sDeconvolved.ome.tif',savename),...
    'dimensionOrder', 'XYZTC','BigTiff',true,'Compression', 'LZW');
end

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