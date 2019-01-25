
%test tracking/linkage of detected particles
%run via: result = runtests('jonathanTrackingTest'); 
%JH 2019-01-24
%%%%%%%%%%%%%%%%
%setup
clear all;
close all;

%note that an alternative (possibly more minimal) would be to use simulated
%data

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

%perform tracking, with verbose option on
job.dataStruct{channel} = kitGenerateTracks(job.dataStruct{channel},1);

%%%%%%%%%%%%%%%%%%%%
% to write tests specify expectations 
%%%%%%%%%%%%%%%%%%%

%% Test1: check if completed and some basic things
assert(~job.dataStruct{channel}.failed,'The job failed')
assert(length(job.dataStruct{channel}.dataProperties.psfSigma)==3, ...
    'psf has length 3')

%% Test2:
assert(length(job.dataStruct{1}.tracks)>=1, ...
    'At least one kinetichore tracked')



