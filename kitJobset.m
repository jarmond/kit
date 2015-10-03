function jobset=kitJobset(jobset)
% KITJOBSET Create or upgrade a jobset.

if nargin<1
  jobset = kitDefaultOptions();
end

if jobset.jobsetVersion < 6
  % Upgrade to use ROIs.
  [mFiles,~,mIdxs] = unique(jobset.movieFiles);
  for i=1:length(jobset.crop)
    jobset.ROI(i).movieIdx = mIdxs(i);
    jobset.ROI(i).crop = jobset.crop{i};
    jobset.ROI(i).cropSize = jobset.cropSize{i};
  end
  jobset = rmfield(jobset,{'crop','cropSize'});
  jobset.movieFiles = mFiles;
end

% Copy any missing fields.
defJob = kitDefaultOptions();
jobset = structCopyMissingFields(jobset,defJob);

