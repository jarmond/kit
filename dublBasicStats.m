function dublBasicStats(intraStructure,varargin)
% DUBLBASICSTATS Gives a list of statistics for delta measurements.
%
%    DUBLBASICSTATS(COMPILEDDELTA,...) Calculates the median, mean and
%    standard error and deviation of measurements of delta collated in
%    INTRASTRUCTURE and outputs to the command line. INTRASTRUCTURE is
%    produced by dublIntraMeasurements.
%
%    Options, defaults in {}:-
%
%    coordSystem: {'plate'} or 'microscope'. Which coordinate system in
%       which to provide statistics. Note that some statistics are
%       coordinate-independent.
%
%    depthFilter: 0 or {1}. Whether or not to give depth-filtered
%       measurements of intra-measurements.
%
%    stat: {'delta3D'} or one of the following:
%           - 'delta3D'
%           - 'delta2D'
%           - 'delta1D'
%           - 'deltaXYZ'
%           - 'sisSep3D'
%           - 'sisSep2D'
%           - 'sisSepXYZ'
%           - 'twist3D'
%           - 'twistYZ'
%           - 'swivel3D'
%           - 'swivelYZ'
%           - 'swivelKMT'
%       The statistic to be printed to screen. If no statistic is provided,
%       the user will be prompted.
%
% Copyright (c) 2016 C. A. Smith


% default options
opts.coordSystem = 'plate';
opts.depthFilter = 1;
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

statsList = {'delta3D' ,'delta2D' ,'delta1D'  ,'deltaXYZ',...
             'sisSep3D','sisSep2D','sisSepXYZ',...
             'twist3D' ,'twistYZ' ,...
             'swivel3D','swivelYZ','swivelKMT'};
if ~ismember(opts.stat,statsList)
    fprintf('Output statistics options:\n')
    for iStat=1:length(statsList)
        fprintf('    %i) %s\n',iStat,statsList{iStat});
    end
    prompt = sprintf('Please type the number statistic you would like: ');
    result = input(prompt);
    opts.stat = statsList{result};
end

switch opts.stat
    
    case 'delta3D'
        if opts.depthFilter
            delta3D = intraStructure.plate.depthFilter.delta.threeD.all(:);
        else
            delta3D = intraStructure.plate.raw.delta.threeD.all(:);
        end
        
        median = nanmedian(delta3D)*1000;
        mean   = nanmean(delta3D)*1000;
        stdErr = nanserr(delta3D)*1000;
        stdDev = nanstd(delta3D)*1000;
        n      = min(sum(~isnan(delta3D)));

        fprintf('\n3D delta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'delta2D'
        if opts.depthFilter
            delta2D = intraStructure.plate.depthFilter.delta.twoD.all(:);
        else
            delta2D = intraStructure.plate.raw.delta.twoD.all(:);
        end   

        median = nanmedian(delta2D)*1000;
        mean   = nanmean(delta2D)*1000;
        stdErr = nanserr(delta2D)*1000;
        stdDev = nanstd(delta2D)*1000;
        n      = min(sum(~isnan(delta2D)));

        fprintf('\n2D delta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'deltaXYZ'
        if opts.depthFilter
            deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x.all(:);
            deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y.all(:);
            deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z.all(:);
        else
            deltaXYZ(:,1) = intraStructure.plate.raw.delta.x.all(:);
            deltaXYZ(:,2) = intraStructure.plate.raw.delta.y.all(:);
            deltaXYZ(:,3) = intraStructure.plate.raw.delta.z.all(:);
        end
    
        median = nanmedian(deltaXYZ)*1000;
        mean   = nanmean(deltaXYZ)*1000;
        stdErr = nanserr(deltaXYZ)*1000;
        stdDev = nanstd(deltaXYZ)*1000;
        n      = min(sum(~isnan(deltaXYZ)));

        fprintf('\nX, Y and Z delta measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f, %.2f] nm\n',median(1),median(2),median(3))
        fprintf('    mean    = [%.2f, %.2f, %.2f] nm\n',mean(1),mean(2),mean(3))
        fprintf('    std err = [%.2f, %.2f, %.2f] nm\n',stdErr(1),stdErr(2),stdErr(3))
        fprintf('    std dev = [%.2f, %.2f, %.2f] nm\n',stdDev(1),stdDev(2),stdDev(3))
        
    otherwise
        
        error('Statistic requested is not yet built into this version of dublBasicStats. See later release.')
    
end

fprintf('\n')

end