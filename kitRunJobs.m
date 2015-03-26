function kitRunJobs(jobset,varargin)
% KITRUNJOBS Runs tracking analysis on jobset
%
%    KITRUNJOBS(JOBSET,...) Runs tracking analysis on JOBSET. JOBSET should be
%    output from either KITSETUPJOBS or KITLOADJOBSET. Additional options can
%    be supplied as string/value pairs.
%
%    Options (default in {}):-
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

nMovies = length(jobset.movieFiles);

% Default options.
options.subset = 1:nMovies;
options.parallel = 0;
options.errorfail = 0;
% Get user options.
options = processOptions(options, varargin{:});

% Check options.
if ~all(ismember(options.subset,1:nMovies))
  error('Subset values must be in range 1 to %d',nMovies);
end

% If using matlabpool for parallel computation, report workers.
if options.parallel
  kitLog('Running in parallel');
else
  kitLog('Running serially');
end

% Upgrade jobset options, if required.
if ~isfield(jobset.options,'jobsetVersion') || ...
    jobset.options.jobsetVersion < kitVersion(2)
  defJob = kitDefaultOptions();
  jobset = structCopyMissingFields(jobset,defJob);
end

% Copy out job info for each movie.
jobs = cell(nMovies,1);
for i=1:nMovies
  jobs{i} = jobset;
  jobs{i}.movie = jobset.movieFiles{i};
  jobs{i}.index = i;
  jobs{i}.crop = jobset.crop{i};
  jobs{i}.cropSize = jobset.cropSize{i};
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

exceptions = [];
for i = options.subset
  if options.parallel
    kitLog('Submitting tracking job %d', i);
    batchJob{i} = batch(@kitTrackMovie, 1, {jobs{i}});
  else
    try
      kitLog('Tracking movie %d', i);
      kitTrackMovie(jobs{i});
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
diagfile = fullfile(pathstr,['diags_' name '.txt']);
fid = fopen(diagfile,'wt');
for c = 1:length(jobset.options.coordMode)
  if ~strcmp(jobset.options.coordMode{c}, 'none')
    kitJobsetDiagnostics(jobset,c,0,fid);
  end
end
fclose(fid);
