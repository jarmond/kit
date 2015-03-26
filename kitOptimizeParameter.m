function diagnostics=kitOptimizeParameter(jobset,jobId,parameter,values,channel)
% KITOPTIMIZEPARAMETER Run tracking over sweep of parameter for optimization
%
%    DIAGNOSTICS = KITOPTIMIZEPARAMETER(JOBSET,JOBID,PARAMETER,VALUES,CHANNEL)
%    Run tracking for job numbered JOBID from JOBSET, over sweep of PARAMETER
%    across VALUES for optimization. PARAMETER should be a string representing a
%    valid field name of JOB.OPTIONS. VALUES should be a vector containing the
%    values for PARAMETER, or a matrix where each row is a value for PARAMETER.
%
%    CHANNEL Optional. Specify tracking channel, default = 1.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<5
  channel = 1; % tracking channel
end

if ~ischar(parameter) || ~isfield(jobset.options,parameter)
  error('PARAMETER should be a string representing field of JOB.OPTIONS');
end

% Ensure column vector, if vector.
if isvector(values) && size(values,2) > 1
  values = values';
end

% Upgrade jobset options, if required.
if ~isfield(jobset.options,'jobsetVersion') || ...
    jobset.options.jobsetVersion < kitVersion(2)
  defJob = kitDefaultOptions();
  jobset = structCopyMissingFields(jobset,defJob);
end

nValues = size(values,1);
for i = 1:nValues
  kitLog('Setting option %s to %s',parameter,num2str(values(i),'%g '));
  jobset.options.(parameter) = values(i,:);

  % Create job.
  job = jobset;
  job.movie = jobset.movieFiles{jobId};
  job.index = jobId;
  job.crop = jobset.crop{jobId};
  job.cropSize = jobset.cropSize{jobId};
  % Update versions, may be different to jobset creator.
  job.version = kitVersion();
  job.matlabVersion = version;
  % Inhibit saving, since all jobs would use same file.
  job.options.disableSave = 1;

  % Submit tracking job.
  kitLog('Submitting tracking job %d', i);
  batchJob(i) = batch(@kitTrackMovie, 1, {job});
end

% Wait for tracking jobs.
for i = 1:nValues
  kitLog('Waiting for %d tracking jobs to complete', nValues-i+1);
  b = batchJob(i);
  wait(b);
  diary(b);

  if strcmp(b.State,'finished')
    % Get result.
    r = fetchOutputs(b);
    job = r{1};
    diagnostics(i) = job.dataStruct{channel}.diagnostics;
  else
    kitLog('Job %d failed', i);
  end

  delete(b);
end

% Save output.
[~,jobsetName] = fileparts(job.filename);
[moviePath,movieName] = fileparts(job.movie);
outputName = fullfile(job.movieDirectory,moviePath,...
                      ['diagnostics-' jobsetName '-' parameter '.mat']);
save(outputName,'diagnostics','parameter','values');
