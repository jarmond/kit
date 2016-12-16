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
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    delta3D = intraStructure.plate.depthFilter.delta.threeD.all(:);
                else
                    delta3D = intraStructure.plate.raw.delta.threeD.all(:);
                end
            case 'microscope'
                if opts.depthFilter
                    delta3D = intraStructure.microscope.depthFilter.delta.threeD.all(:);
                else
                    delta3D = intraStructure.microscope.raw.delta.threeD.all(:);
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
                    delta2D = intraStructure.plate.depthFilter.delta.twoD.all(:);
                else
                    delta2D = intraStructure.plate.raw.delta.twoD.all(:);
                end
            case 'microscope'
                if opts.depthFilter
                    delta2D = intraStructure.microscope.depthFilter.delta.twoD.all(:);
                else
                    delta2D = intraStructure.microscope.raw.delta.twoD.all(:);
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
                    deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x.all(:);
                    deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y.all(:);
                    deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z.all(:);
                else
                    deltaXYZ(:,1) = intraStructure.plate.raw.delta.x.all(:);
                    deltaXYZ(:,2) = intraStructure.plate.raw.delta.y.all(:);
                    deltaXYZ(:,3) = intraStructure.plate.raw.delta.z.all(:);
                end
            case 'microscope'
                if opts.depthFilter
                    deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x.all(:);
                    deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y.all(:);
                    deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z.all(:);
                else
                    deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x.all(:);
                    deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y.all(:);
                    deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z.all(:);
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
                sisSep3D = intraStructure.plate.sisSep.threeD.all(:);
            case 'microscope'
                sisSep3D = intraStructure.microscope.sisSep.threeD.all(:);
        end
        
        median = nanmedian(sisSep3D)*1000;
        mean   = nanmean(sisSep3D)*1000;
        stdErr = nanserr(sisSep3D)*1000;
        stdDev = nanstd(sisSep3D)*1000;
        n      = min(sum(~isnan(sisSep3D)));

        fprintf('\n3D sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'sisSep2D'
            
        switch opts.coordSystem
            case 'plate'
                sisSep2D = intraStructure.plate.sisSep.twoD.all(:);
            case 'microscope'
                sisSep2D = intraStructure.microscope.sisSep.twoD.all(:);
        end
        
        median = nanmedian(sisSep2D)*1000;
        mean   = nanmean(sisSep2D)*1000;
        stdErr = nanserr(sisSep2D)*1000;
        stdDev = nanstd(sisSep2D)*1000;
        n      = min(sum(~isnan(sisSep2D)));

        fprintf('\n2D sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'sisSepXYZ'
        
        switch opts.coordSystem
            case 'plate'
                sisSepXYZ(:,1) = intraStructure.plate.sisSep.x.all(:);
                sisSepXYZ(:,2) = intraStructure.plate.sisSep.y.all(:);
                sisSepXYZ(:,3) = intraStructure.plate.sisSep.z.all(:);
            case 'microscope'
                sisSepXYZ(:,1) = intraStructure.microscope.sisSep.x.all(:);
                sisSepXYZ(:,2) = intraStructure.microscope.sisSep.y.all(:);
                sisSepXYZ(:,3) = intraStructure.microscope.sisSep.z.all(:);
        end

        median = nanmedian(sisSepXYZ)*1000;
        mean   = nanmean(sisSepXYZ)*1000;
        stdErr = nanserr(sisSepXYZ)*1000;
        stdDev = nanstd(sisSepXYZ)*1000;
        n      = min(sum(~isnan(sisSepXYZ)));

        fprintf('\nX, Y and Z sister separation measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f, %.2f] nm\n',median(1),median(2),median(3))
        fprintf('    mean    = [%.2f, %.2f, %.2f] nm\n',mean(1),mean(2),mean(3))
        fprintf('    std err = [%.2f, %.2f, %.2f] nm\n',stdErr(1),stdErr(2),stdErr(3))
        fprintf('    std dev = [%.2f, %.2f, %.2f] nm\n',stdDev(1),stdDev(2),stdDev(3))

    case 'twist3D'
        
        switch opts.coordSystem
            case 'plate'
                twist3D = intraStructure.plate.twist.threeD.all(:);
            case 'microscope'
                error('Twist is not defined for a non-plate coordinate system.')
        end
        
        median = nanmedian(twist3D)*1000;
        mean   = nanmean(twist3D)*1000;
        stdErr = nanserr(twist3D)*1000;
        stdDev = nanstd(twist3D)*1000;
        n      = min(sum(~isnan(twist3D)));

        fprintf('\n3D twist measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'twistYZ'
        
        switch opts.coordSystem
            case 'plate'
                twistYZ(:,1) = intraStructure.plate.twist.y.all(:);
                twistYZ(:,2) = intraStructure.plate.twist.z.all(:);
            case 'microscope'
                error('Twist is not defined for a non-plate coordinate system.')
        end

        median = nanmedian(twistYZ)*1000;
        mean   = nanmean(twistYZ)*1000;
        stdErr = nanserr(twistYZ)*1000;
        stdDev = nanstd(twistYZ)*1000;
        n      = min(sum(~isnan(twistYZ)));

        fprintf('\nY and Z twist measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] nm\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] nm\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] nm\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] nm\n',stdDev(1),stdDev(2))

    case 'swivel3D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    swivel3D = intraStructure.plate.depthFilter.swivel.threeD.all(:);
                else
                    swivel3D = intraStructure.plate.raw.swivel.threeD.all(:);
                end
            case 'microscope'
                if opts.depthFilter
                    swivel3D = intraStructure.microscope.depthFilter.swivel.threeD.all(:);
                else
                    swivel3D = intraStructure.microscope.raw.swivel.threeD.all(:);
                end
        end
        
        median = nanmedian(swivel3D)*1000;
        mean   = nanmean(swivel3D)*1000;
        stdErr = nanserr(swivel3D)*1000;
        stdDev = nanstd(swivel3D)*1000;
        n      = min(sum(~isnan(swivel3D)));

        fprintf('\n3D swivel measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'swivelKMT'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    swivelKMT = intraStructure.plate.depthFilter.swivel.kMT(:);
                else
                    swivelKMT = intraStructure.plate.raw.swivel.kMT(:);
                end
            case 'microscope'
                if opts.depthFilter
                    swivelKMT = intraStructure.microscope.depthFilter.swivel.kMT(:);
                else
                    swivelKMT = intraStructure.microscope.raw.swivel.kMT(:);
                end
        end

        median = nanmedian(swivelKMT)*1000;
        mean   = nanmean(swivelKMT)*1000;
        stdErr = nanserr(swivelKMT)*1000;
        stdDev = nanstd(swivelKMT)*1000;
        n      = min(sum(~isnan(swivelKMT)));

        fprintf('\nKinetochore-to-microtubule angle measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f nm\n',median)
        fprintf('    mean    = %.2f nm\n',mean)
        fprintf('    std err = %.2f nm\n',stdErr)
        fprintf('    std dev = %.2f nm\n',stdDev)
        
    case 'swivelYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    swivelYZ(:,1) = intraStructure.plate.depthFilter.swivel.x.all(:);
                    swivelYZ(:,2) = intraStructure.plate.depthFilter.swivel.y.all(:);
                    swivelYZ(:,3) = intraStructure.plate.depthFilter.swivel.z.all(:);
                else
                    swivelYZ(:,1) = intraStructure.plate.raw.swivel.x.all(:);
                    swivelYZ(:,2) = intraStructure.plate.raw.swivel.y.all(:);
                    swivelYZ(:,3) = intraStructure.plate.raw.swivel.z.all(:);
                end
            case 'microscope'
                if opts.depthFilter
                    swivelYZ(:,1) = intraStructure.microscope.depthFilter.swivel.x.all(:);
                    swivelYZ(:,2) = intraStructure.microscope.depthFilter.swivel.y.all(:);
                    swivelYZ(:,3) = intraStructure.microscope.depthFilter.swivel.z.all(:);
                else
                    swivelYZ(:,1) = intraStructure.microscope.raw.swivel.x.all(:);
                    swivelYZ(:,2) = intraStructure.microscope.raw.swivel.y.all(:);
                    swivelYZ(:,3) = intraStructure.microscope.raw.swivel.z.all(:);
                end
        end
    
        median = nanmedian(swivelYZ)*1000;
        mean   = nanmean(swivelYZ)*1000;
        stdErr = nanserr(swivelYZ)*1000;
        stdDev = nanstd(swivelYZ)*1000;
        n      = min(sum(~isnan(swivelYZ)));

        fprintf('\nY and Z swivel measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] nm\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] nm\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] nm\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] nm\n',stdDev(1),stdDev(2))
        
    otherwise
        
        error('Statistic requested is not yet built into this version of dublBasicStats. See later release.')
    
end

fprintf('\n')

end