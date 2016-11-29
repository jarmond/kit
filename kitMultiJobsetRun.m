function kitMultiJobsetRun(jobsets,exec)
% KITMULTIJOBSETRUN Runs all jobs from multiple jobsets in parallel or serially.
%
% Accepts cell array of filenames or cell array of jobset structs.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2016 C. A. Smith

if nargin<2
  exec = 'serial';
end

% Ensure input is cell array.
if ~iscell(jobsets)
  jobsets = {jobsets};
end

nJobsets = length(jobsets);
kitLog('Running %d jobsets',nJobsets);
for i = 1:nJobsets
  jobset = jobsets{i};
  if ischar(jobset)
    jobset = kitLoadJobset(jobset);
  end
  kitLog('Running jobset: %s',jobset.filename);

  % Run it.
  kitRunJobs(jobset,'exec',exec);
end
