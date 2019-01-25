
%minimal experimentation to test different detection methods on the LLSM
%data
%run via: result = runtests('jonathanDetectionTest'); 
%Perhaps it would actually be OK to relax the requirements on upper limits
%of biologically realistic numbers of spots, since here we are detecting
%CANDIDATE spots, which are later refined
%JH 2019-01-16
%%%%%%%%%%%%%%%%
%setup
clear all;
close all;

%load in deconvolved data/jobset for testing
%jobset = kitLoadJobset('../../Data/Lattice light sheet/kitjobset_190116_kitDev_v002.mat');
jobset = kitLoadJobset('../../Data/Lattice light sheet/kittracking002-kitjobset_190115_centroidv001-OS_181205_MC139_LatticeLightSheet_4_decon_5bg.mat');
%jobset = kitLoadJobset('../../Data/Lattice light sheet/kitjobset_190122_raw_adaptive_mmf_NA14.mat');
job = jobset;
job.index=1;
channel = 1;

if isfield(job,'metadata')
    if iscell(job.metadata)
        job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
end
job = kitSaveJob(job);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%load movie
movie = kitReadWholeMovie(reader,job.metadata,channel,job.ROI.crop,0,1);
[imageSizeX,imageSizeY,imageSizeZ,nFrames] = size(movie);
%setup properties
ds = kitMakeMakiDatastruct(job,1); %in channel 1 only
job.dataStruct{1}.dataProperties = ds.dataProperties;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% spots = cell(nFrames,1);
% for i=1:nFrames
%     img = movie(:,:,:,i);
%     spots{i} = histcutSpots(img,job.options,job.dataStruct{channel}.dataProperties);
% end
job.options.adaptiveLambda = 10; %lambda=10 for deconvolved, lambda=-20 for raw
spots = adaptiveSpots(movie,job.options.adaptiveLambda,1);

% job.options.wavelet.numLevels=10;
% job.options.wavelet.localMAD=10; % locally estimated MAD
% job.options.wavelet.backSub=10;  % background subtraction
% job.options.wavelet.prefilter=10; % denoise prefilter
% job.options.wavelet.minLevel=10;
% %thresh = waveletAdapt(movie,job.options); %2
% job.options.waveletLevelThresh = 2;
% spots = cell(nFrames,1);
% for i=1:nFrames
%     img = movie(:,:,:,i);
%     spots{i} = waveletSpots(img,job.options);
% end

%max project and show spots
figure;
imshow(max(movie(:,:,:,1),[],3));
hold on;
spots_max_project = max(spots{1},[],3);
scatter(spots_max_project(:,2),spots_max_project(:,1));

%project image to see line plot
projection = max(movie(:,29,10,1),[],3);
figure;
plot(projection)

%%%%%%%%%%%%%%%
% to write tests specify expectations of length of spots, and what spots
% should be
%%%%%%%%%%%%%%%%%%%
%% Test1: spots has the right dimensions

assert(length(spots)==size(movie,4),'Missing spots output for some frames')
assert(size(spots{1},2)==3,'Expecting a 3D movie to test on')
assert(size(spots{1},1)>0, ...
'Expect to find at least some kinetochores in the first frame')
    
%% Test2: Found roughly the right amount of spots
nSpotsPerFrame = zeros(size(spots));
for j = 1:50
    nSpotsPerFrame(j) = size(spots{j},1);
end
assert(min(nSpotsPerFrame)>job.options.minSpotsPerFrame, ...
    sprintf('Expect a biologically reasonable number of spots per frame, found min of %d',min(nSpotsPerFrame)));
assert(max(nSpotsPerFrame)<job.options.maxSpotsPerFrame, ...
    sprintf('Expect a biologically reasonable number of spots per frame, found max of %d',max(nSpotsPerFrame)));
assert(min(nSpotsPerFrame)>0.5*96, ...
    sprintf('Expect a biologically reasonable number of spots per frame, found min of %d',min(nSpotsPerFrame)));
assert(max(nSpotsPerFrame)<1.5*96, ...
    sprintf('Expect a biologically reasonable number of spots per frame, found max of %d',max(nSpotsPerFrame)));
%% Test3: does a spot seem like a genuine spot
assert(all(spots{1}(1,:)>0), 'Spot coords should be positive')
assert(all(spots{1}(1,:) <= ds.dataProperties.movieSize(1:3)), ...
    'Spot coords should lie within the movie')
assert(all(movie([spots{1}(1,:),1])>0), ...
    'First spot in first frame should have nonzero amplitude') 
assert(all(movie([spots{1}(1,:),1])>median(median(median(median(movie))))), ...
    'First spot in first frame should have amplitude grater than average')
