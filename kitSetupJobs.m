function jobset=kitSetupJobs(jobset,varargin)
% KITSETUPJOBS Front end to kinetochore tracker
%
%    JOBSET=KITSETUPJOBS(JOBSET) Front end to kinetochore tracker. Queries
%    user for movies and parameters and returns JOBSET populated with data
%    required to track movies. A previous JOBSET may be supplied to change
%    existing options.
%
%    Options (default in {}):-
%
%    forcecrop: {0} or 1. Set -1 to force redisplay of crop GUI for all movies, otherwise specify movie numbers.
%
% Copyright (c) 2015 Jonathan W. Armond

options.forcecrop = 0;
% Get user options.
options = processOptions(options, varargin{:});

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
nMovies = length(jobset.movieFiles);
if ~isfield(jobset,'crop') || options.forcecrop~=0 || numel(jobset.crop)<nMovies ...
    || any(cellfun(@isempty,jobset.crop))
  if ~isfield(jobset,'crop')
    jobset.crop = {};
    jobset.cropSize = {};
  end

  mvs = jobset.movieFiles;
  origCrops= jobset.crop;
  nOrigMovies = nMovies;
  for i=1:nOrigMovies
    if options.forcecrop == -1 || ismember(i,options.forcecrop) || numel(origCrops) < i || isempty(origCrops)
      [crop, cropSize] = ...
          kitCropMovie(fullfile(jobset.movieDirectory,mvs{i}));
      jobset.crop{i} = crop(1,:);
      jobset.cropSize{i} = cropSize(1,:);

      % For multiple ROIs, generate duplicate movieFiles.
      for j=2:size(crop,1)
        nMovies = nMovies+1;
        jobset.movieFiles{nMovies} = jobset.movieFiles{i};
        jobset.crop{nMovies} = crop(j,:);
        jobset.cropSize{nMovies} = cropSize(j,:);
      end
    end

    % Save progress each time.
    kitSaveJobset(jobset);
  end
end

kitSaveJobset(jobset);
