function job = jonathanExamplePlaneFitTest(job,makeMovie)
% Perform and test plane fitting
%
%take output from jonathanExampleDetectionTest.m and
%jonathanExampleRefinementTest.m rather than redoing
% Jonathan U Harrison 2019-02-12
%%%%%%%%%%

%%%%% can get output from metaphase example by rerunning previous test
if nargin <1
    [spots, movie, job] = jonathanExampleDetectionTest();
    job = jonathanExampleRefinementTest(spots,movie,job);
end

if nargin<2
    makeMovie=0; %make a movie of the planes for each frame
end

%read in movie
if isfield(job,'metadata')
    if iscell(job.metadata)
        job.metadata = job.metadata{job.index};
    end
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie),'valid',job.metadata);
else
    [job.metadata, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.ROI.movie));
end

% Fit plane in chosen channel.
  planeChan = job.options.coordSystemChannel;
  job.options.debug.showPlaneFit = 2;
  job.options.debug.makePlaneFitMovie = makeMovie;
  kitLog('Fitting plane in channel %d', planeChan);
  if strcmp(job.options.coordMode{planeChan}, 'none')
    % No spot tracking in plane channel so populate dataStruct.
    job.dataStruct{planeChan} = kitMakeMakiDatastruct(job, planeChan);
  end
  if strcmp(job.options.coordSystem, 'register')
    job.dataStruct{planeChan} = kitRegisterFrames(job,reader,job.dataStruct{planeChan});
  else
    job.dataStruct{planeChan} = kitFitPlane(job,reader,job.dataStruct{planeChan},planeChan,0);
  end

  %% Test1: check if reasonable number of spots
assert(~job.dataStruct{planeChan}.failed, 'Make sure fit has not failed');
assert(size(job.dataStruct{1}.planeFit(1).plane,2)==4, ...
    'Check size of plane for first time slice');
assert(all(job.dataStruct{1}.planeFit(1).eigenValues>0), ...
    'Check all eigenvalues positive for first plane');
  
end