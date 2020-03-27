function jobset = jonathanPipelineJobsetSetup(movieDir,identifier,dt)
%%make custom kit jobset files for tracking in pipeline
%% in most cases not a replacement for using kitSetupJob
%% but helpful when we want to automate creating many standard jobset files
%%J.U. harrison 2020-01-28
%%
%%%%%%%%

if nargin<2
    identifier = 'DonaldDuck_auto_v101';
end

if nargin<3
    dt = 4.7;
end

jobset = kitDefaultOptions();
jobset.options.jobProcess='zandt';

if ~isfield(jobset,'jobsetVersion') || ...
        jobset.jobsetVersion < kitVersion(2)
    jobset = kitJobset(jobset);
end

if ~isfield(jobset,'ROI')
    jobset.ROI = [];
end

jobset.movieDirectory = movieDir; %need to get external input from bash script?
movieFiles = kitFindFiles(movieDir, 'flowdec_deconvolved.ome.tif',1,0,1); %will find deconvolved movies
shortenedMovieFiles = strrep(movieFiles,[movieDir filesep],''); %remove folder part from path. NB try to make sure it does not already end in a /
jobset.movieFiles = shortenedMovieFiles; 

for i=1:length(jobset.movieFiles)
    [md,~]=kitOpenMovie(movieFiles{i}); %full path
    crop = [1 1 md.frameSize(1:2)];
    cropSize = md.frameSize(1:3);
    r = length(jobset.ROI) + 1;
    jobset.ROI(r).movie = jobset.movieFiles{i}; %short path
    jobset.ROI(r).crop = crop;
    jobset.ROI(r).cropSize = cropSize;
end

%compute avaerage search radius
avgDisp = 3.6/60; %um/min -> um/s
r = zeros(2,1);
r(2) = 6*avgDisp*dt;
r(1) = 0.1*avgDisp*dt;
%set
jobset.options.autoRadiidt = dt;
jobset.options.autoRadiiAvgDisp = 3.6/60;
jobset.options.minSearchRadius = [r(1) 3*r(1) r(1)];
jobset.option.maxSearchRadius = [r(2) 3*r(2) r(2)]; %[0.75,3,0.75];

jobset.options.flatBackground = 1;
jobset.options.spotMode{1} = 'adaptive';
jobset.options.mmf.maxMmfTime = -1; %no maximum time

jobset.filename = ['kitjobset_' datestr(now,'yymmdd') '_' identifier '.mat'];
jobset.filename = fullfile(jobset.movieDirectory,jobset.filename);
kitSaveJobset(jobset);
end


