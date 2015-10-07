function kitRunJobs(jobset,varargin)
% KITRUNJOBS Runs tracking analysis on jobset
%
%    KITRUNJOBS(JOBSET,...) Runs tracking analysis on JOBSET. JOBSET should be
%    output from either KITSETUPJOBS or KITLOADJOBSET. Additional options can
%    be supplied as string/value pairs.
%
%    Options (default in {}):-
%
%    tasks: {1:7} or a vector. Specify tasks to perform on jobs:
%                              1: finding spots
%                              2: fitting plane
%                              3: tracking spots
%                              4: grouping sisters
%                              5: extracting tracks
%                              6: updating classes
%                              7: aligning
%                              8: intensity
%
%    existing: {0} or 1. Load existing jobs first.
%
%    parallel: {0} or 1. Set 1 to use multiple processors to parallel
%    process.
%
%    subset: {[]} or vector of job numbers. Set to a subset of indices of movies
%    to analysis, instead of processing them all.
%
%    errorfail: {0} or 1. Crash out if error occurs if 1.
%
% Copyright (c) 2013 Jonathan W. Armond

% Check minimum MATLAB version.
% FIXME Check minimum toolbox versions also.
if verLessThan('matlab','7.14')
  error('Minimum required MATLAB version is 7.14 (R2012a)');
end

% Download BioFormats, if required.
kitDownloadBioFormats();

% Upgrade jobset options, if required.
if ~isfield(jobset.options,'jobsetVersion') || ...
    jobset.options.jobsetVersion < kitVersion(2)
  jobset = kitJobset(jobset);
end

nROIs = length(jobset.ROI);

% Default options.
options.subset = 1:nROIs;
options.parallel = 0;
options.errorfail = 0;
options.tasks = 1:7;
options.existing = 0;
options.callback = [];
% Get user options.
options = processOptions(options, varargin{:});

% Check options.
if ~all(ismember(options.subset,1:nROIs))
  error('Subset values must be in range 1 to %d',nROIs);
end

% If using matlabpool for parallel computation, report workers.
[~,name] = fileparts(jobset.filename);
if options.parallel
  kitLog(['Running ' name ' in parallel']);
else
  kitLog(['Running ' name ' serially']);
end


if isfield(jobset,'variantName')
  fprintf('Jobset variant: %s\n',jobset.variantName);
end

% Copy out job info for each movie.
jobs = cell(nROIs,1);
if options.existing
  for i=1:nROIs
    jobs{i} = kitLoadJob(jobset,i);
  end
else
  for i=1:nROIs
    jobs{i} = jobset;
    jobs{i}.movie = jobset.ROI(i).movie;
    jobs{i}.index = i;
    jobs{i}.crop = jobset.ROI(i).crop;
    jobs{i}.cropSize = jobset.ROI(i).cropSize;
    % Update versions, may be different to jobset creator.
    jobs{i}.version = kitVersion();
    jobs{i}.matlabVersion = version;
    %jobs{i}.lociVersion = char(loci.formats.FormatTools.VERSION);
    % Record host.
    if ispc
      [~,jobs{i}.host] = system('echo %COMPUTERNAME%');
    else
      [~,jobs{i}.host] = system('hostname');
    end
  end
end


exceptions = [];
for i = options.subset
  if options.parallel
    kitLog('Submitting tracking job %d', i);
    batchJob{i} = batch(@kitTrackMovie, 1, {jobs{i},options.tasks});
  else
    try
      kitLog('Tracking movie %d', i);
      kitTrackMovie(jobs{i},options.tasks);
      if ~isempty(options.callback)
        options.callback(i);
      end
    catch me
      kitLog('Error in movie %d: %s',i,me.identifier);
      ex.me = me;
      ex.idx = i;
      exceptions = [exceptions ex];
      if options.errorfail
        disp(getReport(me));
        throw(me);
      end
    end
  end
end

if ~isempty(exceptions)
  disp('Errors occured:')
end
for i = 1:length(exceptions)
  ex = exceptions(i);
  fprintf('In movie %d, error %s:\n',ex.idx,ex.me.identifier);
  disp(getReport(ex.me));
end

if options.parallel
  % Wait for parallel tracking jobs.
  for i = options.subset
    kitLog('Waiting for %d tracking jobs to complete', length(options.subset)-i+1);
    b = batchJob{i};
    wait(b);
    diary(b);
    delete(b);
  end
end

kitLog('Tracking complete');

% Dump jobset diagnostics.
[pathstr,name,ext] = fileparts(jobset.filename);
diagfile = ['diags_' name];
if isfield(jobset,'variantName')
  diagfile = [diagfile '_' jobset.variantName];
end
diagfile = fullfile(pathstr,[diagfile '.txt']);
fid = fopen(diagfile,'wt');
C = onCleanup(@() fclose(fid));
for c = 1:length(jobset.options.coordMode)
  if ~strcmp(jobset.options.coordMode{c}, 'none')
    kitJobsetDiagnostics(jobset,c,0,fid);
  end
end
