
%test refinement of candidate spots
%data
%run via: result = runtests('jonathanRefinementTest');
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

%detect particles
job = kitFindCoords(job,reader,channel);
initCoord = job.dataStruct{1}.initCoord;

%%%%%%%%%%%%%%%
% to write tests specify expectations of length of spots, and what spots
% should be
%%%%%%%%%%%%%%%%%%%

%% Test1: check right number of refined spots
for i=1:job.metadata.nFrames
    assert(initCoord(i).nSpots>job.options.minSpotsPerFrame, ...
        sprintf('Expect a biologically reasonable number of spots per frame, found min of %d',initCoord(i).nSpots));
    assert(initCoord(i).nSpots<job.options.maxSpotsPerFrame, ...
        sprintf('Expect a biologically reasonable number of spots per frame, found max of %d',initCoord(i).nSpots));
    assert(initCoord(i).nSpots>0.5*96, ...
        sprintf('Expect a biologically reasonable number of spots per frame, found min of %d',initCoord(i).nSpots));
    assert(initCoord(i).nSpots<1.5*96, ...
        sprintf('Expect a biologically reasonable number of spots per frame, found max of %d',initCoord(i).nSpots));
end
%% Test2: check spots seem like possibly real spots
for i=1:job.metadata.nFrames
    assert(all(initCoord(i).allCoord(1:3)>0), 'Spot coords should be positive')
    assert(all(initCoord(i).allCoord(1:3) <= ds.dataProperties.movieSize(1:3)), ...
        'Spot coords should lie within the movie')
    assert(all(initCoord(i).amp(:) >= 0), ...
        'Spot amplitudes should be positive')
end
