%more detailed tests of spot refinement
%simulate a test image

%%%%%%%%%%%%%%%
clear all;
close all;

%load from previous jobset
jobString = '../../Data/Lattice light sheet/kittracking001-kitjobset_190116_decon_adaptive_mmf_v003-OS_181205_MC139_LatticeLightSheet_4_decon_5bg.mat';
job = kitLoadJobset(jobString);
job.options.coordMode{1} = 'gaussian';
job.options.debug.mmfVerbose = 1;
job.index=1;
channel=1;

%%%%%%%%%%%%%%%
%extract objects from jobset
frame=11; %look at frame 11
cands = job.dataStruct{1}.initCoord(1).localMaxima(frame).cands;
cands(:,end+1) = NaN; %this would be for detections in a neighbouring channel
psfSigma = job.dataStruct{1}.dataProperties.psfSigma([1,1,2]);

%open movie and get a frame from it
if isfield(job,'metadata')
    if iscell(job.metadata)
        job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
end
movie = kitReadWholeMovie(reader,job.metadata,channel,job.ROI.crop,0,1);
[imageSizeX,imageSizeY,imageSizeZ,nFrames] = size(movie);
image = movie(:,:,:,frame);

%%%%%%%%%%%%%%%%%%
%perform mixture model fitting
fit = mixtureModelFit(cands,image,psfSigma,job.options,0);

%% Test1: refinement should give fewer spots
assert(size(cands,1)>=size(fit,1),'refinement should reduce num of spots');
assert((size(cands,2)==4) | size(cands,2)==3 , ...
    'expect certain size for candidate spots object depending on 2 or 3d');
assert((size(fit,2)==5) | size(fit,2)==6 , ...
    'expect certain size for refined spots depending on 2 or 3d');

