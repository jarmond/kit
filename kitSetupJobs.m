function jobset=kitSetupJobs(jobset)
% KITSETUPJOBS Front end to kinetochore tracker
%
%    JOBSET=KITSETUPJOBS(JOBSET) Front end to kinetochore tracker. Queries
%    user for movies and parameters and returns JOBSET populated with data
%    required to track movies. A previous JOBSET may be supplied to change
%    existing options.
%
% Copyright (c) 2013 Jonathan W. Armond

% Download BioFormats, if required.
kitDownloadBioFormats();

if nargin<1 || isempty(jobset)
  % Populate default options.
  jobset = kitDefaultOptions();

  % Ask for options.
  jobset = kitAskUser(jobset);
  % TODO validate

  jobset.filename = [fullfile(jobset.movieDirectory,jobset.filename) '.mat'];
end

kitSaveJobset(jobset);

% Ask for cropping.
if ~isfield(jobset,'crop')
  mvs = jobset.movieFiles; % copy movieFiles.
  jobset.movieFiles = {};
  jobset.crop = {};
  jobset.cropSize = {};
  nMovies = length(mvs);
  for i=1:nMovies
    [crop, cropSize] = ...
        kitCropMovie(fullfile(jobset.movieDirectory,mvs{i}));
    % Add multiple ROIs.
    for j=1:size(crop,1)
      jobset.movieFiles{end+1} = mvs{i};
      jobset.crop{end+1} = crop(j,:);
      jobset.cropSize{end+1} = cropSize(j,:);
    end
  end
end

kitSaveJobset(jobset);