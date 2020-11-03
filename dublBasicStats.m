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
%    coordSystem: {'plate'} or 'microscope'. The coordinate system in which
%       to provide statistics. Note that some statistics are coordinate-
%       independent.
%
%    depthFilter: 0 or {1}. Whether or not to give depth-filtered
%       measurements of intra-measurements.
%
%    stat: {''} or one of the following:
%           - 'delta3D'
%           - 'delta2D'
%           - 'delta1D'
%           - 'deltaXYZ'
%           - 'sisSep3D'
%           - 'sisSep2D'
%           - 'sisSepXYZ'
%       The statistic to be printed to screen. If no statistic is provided,
%       the user will be prompted.
%
% Copyright (c) 2018 C. A. Smith


% default options
opts.coordSystem = 'microscope';
opts.depthFilter = 1;
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

statsList = {'delta3D' ,'delta2D' ,'delta1D'  ,'deltaXYZ',...
             'sisSep3D','sisSep2D','sisSepXYZ'};
if ~ismember(opts.stat,statsList)
    fprintf('Output statistics options:\n')
    for iStat=1:length(statsList)
        fprintf('    %i) %s\n',iStat,statsList{iStat});
    end
    prompt = sprintf('Please type the number statistic you would like: ');
    result = input(prompt);
    opts.stat = statsList{result};
end
isall = isfield(intraStructure.microscope.raw.delta.threeD,'all');

if strcmp(opts.coordSystem,'plate')
    if (isall && isempty(intraStructure.plate.raw.delta.threeD.all)) || (~isall && isempty(intraStructure.plate.raw.delta.threeD))
        kitLog('No plane fit in movies in intraMeasurements. Converting coordinate system to ''microscope''');
        opts.coordSystem = 'microscope';
    end
end

switch opts.stat

    case 'delta3D'
            
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        delta3D = intraStructure.plate.depthFilter.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.plate.depthFilter.delta.threeD(:);
                    end
                else
                    if isall
                        delta3D = intraStructure.plate.raw.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.plate.raw.delta.threeD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        delta3D = intraStructure.microscope.depthFilter.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.microscope.depthFilter.delta.threeD(:);
                    end
                else
                    if isall
                        delta3D = intraStructure.microscope.raw.delta.threeD.all(:);
                    else
                        delta3D = intraStructure.microscope.raw.delta.threeD(:);
                    end
                end
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
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        delta2D = intraStructure.plate.depthFilter.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.plate.depthFilter.delta.twoD(:);
                    end
                else
                    if isall
                        delta2D = intraStructure.plate.raw.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.plate.raw.delta.twoD(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        delta2D = intraStructure.microscope.depthFilter.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.microscope.depthFilter.delta.twoD(:);
                    end
                else
                    if isall
                        delta2D = intraStructure.microscope.raw.delta.twoD.all(:);
                    else
                        delta2D = intraStructure.microscope.raw.delta.twoD(:);
                    end
                end
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
        
    case 'delta1D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    delta1D = intraStructure.plate.depthFilter.delta.oneD(:);
                else
                    delta1D = intraStructure.plate.raw.delta.oneD(:);
                end   
            case 'microscope'
                if opts.depthFilter
                    delta1D = intraStructure.microscope.depthFilter.delta.oneD(:);
                else
                    delta1D = intraStructure.microscope.raw.delta.oneD(:);
                end
        end

        median = nanmedian(delta1D)*1000;
        mean   = nanmean(delta1D)*1000;
        stdErr = nanserr(delta1D)*1000;
        stdDev = nanstd(delta1D)*1000;
        n      = min(sum(~isnan(delta1D)));

        fprintf('\n1D delta measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'deltaXYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    if isall
                        deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z(:);
                    end
                else
                    if isall
                        deltaXYZ(:,1) = intraStructure.plate.raw.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.plate.raw.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.plate.raw.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.plate.raw.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.plate.raw.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.plate.raw.delta.z(:);
                    end
                end
            case 'microscope'
                if opts.depthFilter
                    if isall
                        deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z(:);
                    end
                else
                    if isall
                        deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x.all(:);
                        deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y.all(:);
                        deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z.all(:);
                    else
                        deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x(:);
                        deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y(:);
                        deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z(:);
                    end
                end
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
    
    case 'sisSep3D'
        
        switch opts.coordSystem
            case 'plate'
                sisSep3D = intraStructure.plate.sisSep.threeD(:);
            case 'microscope'
                sisSep3D = intraStructure.microscope.sisSep.threeD(:);
        end
        
        median = nanmedian(sisSep3D);
        mean   = nanmean(sisSep3D);
        stdErr = nanserr(sisSep3D);
        stdDev = nanstd(sisSep3D);
        n      = min(sum(~isnan(sisSep3D)));

        fprintf('\n3D sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f um\n',median)
        fprintf('    mean    = %.2f um\n',mean)
        fprintf('    std err = %.2f um\n',stdErr)
        fprintf('    std dev = %.2f um\n',stdDev)
        
    case 'sisSep2D'
            
        switch opts.coordSystem
            case 'plate'
                sisSep2D = intraStructure.plate.sisSep.twoD(:);
            case 'microscope'
                sisSep2D = intraStructure.microscope.sisSep.twoD(:);
        end
        
        median = nanmedian(sisSep2D);
        mean   = nanmean(sisSep2D);
        stdErr = nanserr(sisSep2D);
        stdDev = nanstd(sisSep2D);
        n      = min(sum(~isnan(sisSep2D)));

        fprintf('\n2D sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f um\n',median)
        fprintf('    mean    = %.2f um\n',mean)
        fprintf('    std err = %.2f um\n',stdErr)
        fprintf('    std dev = %.2f um\n',stdDev)
        
    case 'sisSepXYZ'
        
        switch opts.coordSystem
            case 'plate'
                sisSepXYZ(:,1) = intraStructure.plate.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.plate.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.plate.sisSep.z(:);
            case 'microscope'
                sisSepXYZ(:,1) = intraStructure.microscope.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.microscope.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.microscope.sisSep.z(:);
        end

        median = nanmedian(sisSepXYZ);
        mean   = nanmean(sisSepXYZ);
        stdErr = nanserr(sisSepXYZ);
        stdDev = nanstd(sisSepXYZ);
        n      = min(sum(~isnan(sisSepXYZ)));

        fprintf('\nX, Y and Z sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f, %.2f] um\n',median(1),median(2),median(3))
        fprintf('    mean    = [%.2f, %.2f, %.2f] um\n',mean(1),mean(2),mean(3))
        fprintf('    std err = [%.2f, %.2f, %.2f] um\n',stdErr(1),stdErr(2),stdErr(3))
        fprintf('    std dev = [%.2f, %.2f, %.2f] um\n',stdDev(1),stdDev(2),stdDev(3))

    otherwise
        
        error('Statistic requested is not yet built into this version of dublBasicStats. See later release.')
    
end

fprintf('\n')

end

%% Sub-functions

function outs = findoutliers(data)
  if nargin<1 || isempty(data)
    return
  end
  if verLessThan('matlab','9.2')
    nTests = length(data);
    outs = zeros(nTests,1);
    for iTest = 1:nTests
      outs(iTest) = ttest2(data,data(iTest),'alpha',0.0455);
    end    
  else
    outs = isoutlier(data,'mean');
  end
end
