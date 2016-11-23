function [jobset,jobs] = chrsChangeChromaticShift(jobset,newShift,varargin)
%CHRSCHANGECHROMATICSHIFT Changes the chromatic shift measurements given to
%    a certain jobset.
%
%    CHRSCHANGECHROMATICSHIFT(JOBSET,NEWSHIFT,...) Gives a new JOBSET, and
%    JOBS if requested by the user, containing the chromatic shift
%    measurements given in NEWSHIFT.
%
%    Options, defaults in {}:-
%
%    chanVect: {[1 2]} or two numbers between 1 and 3. The two channels for
%         which the new chromatic shift measurement applies, directed from
%         the first to the second number of the vector.
%
%    andJobs: {0} or 1. Whether or not to correct all jobs already
%         processed.
%         NB: This is not recommended for significant changes in chromatic
%         shift measurements.
%
%
% Copyright (c) 2016 C. A. Smith

% default options
opts.chanVect = [];
opts.andJobs = 0;
% process user-defined options
opts = processOptions(opts,varargin{:});

%% Change jobset information

% check whether the chromatic shift provided is in cell form, otherwise
% ensure a channel vector is provided
if isempty(opts.chanVect)
    if ~iscell(newShift)
        error('If providing only one chromatic shift vector, need to specify channels.')
    else
        chrShift = newShift;
        source = jobset.options.chrShift.jobset;
        for fromChan=1:2
          for toChan=2:3
            if fromChan<toChan && sum(chrShift{fromChan,toChan})~=0
              source{fromChan,toChan} = 'Unknown source';
              source{toChan,fromChan} = 'Unknown source';
            else
              source{fromChan,toChan} = [];
              source{toChan,fromChan} = [];
            end
          end
        end
    end
    
% check whether the channel vector includes only allowed channels
elseif ~isempty(setdiff(opts.chanVect,1:3))
    error('Can only choose channels from 1 to 3.')
    
else % if channel vector provided, proceed as appropriate

    fromChan = opts.chanVect(1);
    toChan = opts.chanVect(2);
    
    % get original chromatic shift cell array
    chrShift = jobset.options.chrShift.result;
    source = jobset.options.chrShift.jobset;
    
    % replace the newShift with only the vector specified in chanVect
    if iscell(newShift)
        warning('Replacing only specified channel from new shift structure.')
        newShift = newShift{fromChan,toChan};
    end
    
    % ensure that only the coordinate measurements are replaced
    lenShift = length(newShift);
    chrShift{fromChan,toChan}(1:lenShift) = newShift;
    
    % replace the mirrored measurement (ie. from toChan to fromChan)
    revShift = newShift;
    revShift(1:3) = revShift(1:3)*-1;
    chrShift{toChan,fromChan}(1:lenShift) = revShift;
    
    % give a source for the chromatic shift calculation
    source{fromChan,toChan} = 'Unknown source';
    source{toChan,fromChan} = 'Unknown source';
    
end

% change chromatic shift information in the jobset and save
jobset.options.chrShift.result = chrShift;
jobset.options.chrShift.jobset = source;
kitSaveJobset(jobset);

%% Change jobs information if required

if opts.andJobs
    
    % get job structure
    jobs = kitGenerateMovieStructs(jobset);
    
    % for each job in the jobset, change the chromatic shift information
    for iMov = 1:length(jobs)
        jobs(iMov).options.chrShift.result = chrShift;
        jobs(iMov) = kitSaveJob(jobs(iMov));
    end
else
    jobs = NaN;
end

end




