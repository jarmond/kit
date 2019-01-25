%test plane fit
%run via: result = runtests('jonathanInitialPlaneFitTest');
%JH 2019-01-24
%%%%%%%%%%%%%%%%
%setup
clear all;
close all;

%alternative approach would be to use simulated data

%load in deconvolved data/jobset for testing
jobset = kitLoadJobset('../../Data/Lattice light sheet/kittracking002-kitjobset_190115_centroidv001-OS_181205_MC139_LatticeLightSheet_4_decon_5bg.mat');
job = jobset;
job.options.spotMode{1} = 'adaptive'; %use adaptive and mmf methods
job.options.coordMode{1} = 'gaussian';
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
%add datastruct field
ds = kitMakeMakiDatastruct(job,channel);
job.dataStruct{channel}.dataProperties = ds.dataProperties;

%%%%%%%%%%%%%%%%%%%%

%detect particles; tested elsewhere
job = kitFindCoords(job,reader,channel);

%%%%%%%%%%%%%%%%%%%%%
%plane fit
% Fit plane in chosen channel.
planeChan = job.options.coordSystemChannel;
if strcmp(job.options.coordMode{planeChan}, 'none')
    % No spot tracking in plane channel so populate dataStruct.
    job.dataStruct{planeChan} = kitMakeMakiDatastruct(job, planeChan);
end
if strcmp(job.options.coordSystem, 'register')
    job.dataStruct{planeChan} = kitRegisterFrames(job,reader,job.dataStruct{planeChan});
else
    job.dataStruct{planeChan} = kitFitPlane(job,reader,job.dataStruct{planeChan},planeChan,0);
end

%%%%%%%%%%%%%%%%%%%%
%% Test1: check that a plane has been fit

assert(length(job.dataStruct{channel}.planeFit)>1, ...
    'Check there is something in plane fit field of the struct')
assert(length(job.dataStruct{channel}.planeFit)==job.metadata.nFrames, ...
    'Check output for each frame exists')

