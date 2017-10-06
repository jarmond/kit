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
%    coordSystem: {'plate'}, 'microscope' or 'sisters'. The coordinate system in which
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
%           - 'rawInts'
%           - 'normInts'
%       The statistic to be printed to screen. If no statistic is provided,
%       the user will be prompted.
%
% Copyright (c) 2017 C. A. Smith


% default options
opts.coordSystem = 'plate';
opts.depthFilter = 1;
opts.stat = '';
% get user options
opts = processOptions(opts,varargin{:});

statsList = {'delta3D' ,'delta2D' ,'delta1D'  ,'deltaXYZ',...
             'sisSep3D','sisSep2D','sisSepXYZ',...
             'twist3D' ,'twistYZ' ,...
             'swivel3D','swivelYZ','swivelKMT',...
             'rawInts','normInts'};
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
                    delta3D = intraStructure.plate.depthFilter.delta.threeD(:);
                else
                    delta3D = intraStructure.plate.raw.delta.threeD(:);
                end
            case 'microscope'
                if opts.depthFilter
                    delta3D = intraStructure.microscope.depthFilter.delta.threeD(:);
                else
                    delta3D = intraStructure.microscope.raw.delta.threeD(:);
                end
            case 'sisters'
                if opts.depthFilter
                    delta3D = intraStructure.sisters.depthFilter.delta.threeD(:);
                else
                    delta3D = intraStructure.sisters.raw.delta.threeD(:);
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
                    delta2D = intraStructure.plate.depthFilter.delta.twoD(:);
                else
                    delta2D = intraStructure.plate.raw.delta.twoD(:);
                end
            case 'microscope'
                if opts.depthFilter
                    delta2D = intraStructure.microscope.depthFilter.delta.twoD(:);
                else
                    delta2D = intraStructure.microscope.raw.delta.twoD(:);
                end
            case 'sisters'
                if opts.depthFilter
                    delta2D = intraStructure.sisters.depthFilter.delta.twoD(:);
                else
                    delta2D = intraStructure.sisters.raw.delta.twoD(:);
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
            case 'sisters'
                if opts.depthFilter
                    delta1D = intraStructure.sisters.depthFilter.delta.oneD(:);
                else
                    delta1D = intraStructure.sisters.raw.delta.oneD(:);
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
                    deltaXYZ(:,1) = intraStructure.plate.depthFilter.delta.x(:);
                    deltaXYZ(:,2) = intraStructure.plate.depthFilter.delta.y(:);
                    deltaXYZ(:,3) = intraStructure.plate.depthFilter.delta.z(:);
                else
                    deltaXYZ(:,1) = intraStructure.plate.raw.delta.x(:);
                    deltaXYZ(:,2) = intraStructure.plate.raw.delta.y(:);
                    deltaXYZ(:,3) = intraStructure.plate.raw.delta.z(:);
                end
            case 'microscope'
                if opts.depthFilter
                    deltaXYZ(:,1) = intraStructure.microscope.depthFilter.delta.x(:);
                    deltaXYZ(:,2) = intraStructure.microscope.depthFilter.delta.y(:);
                    deltaXYZ(:,3) = intraStructure.microscope.depthFilter.delta.z(:);
                else
                    deltaXYZ(:,1) = intraStructure.microscope.raw.delta.x(:);
                    deltaXYZ(:,2) = intraStructure.microscope.raw.delta.y(:);
                    deltaXYZ(:,3) = intraStructure.microscope.raw.delta.z(:);
                end
            case 'sisters'
                if opts.depthFilter
                    deltaXYZ(:,1) = intraStructure.sisters.depthFilter.delta.x(:);
                    deltaXYZ(:,2) = intraStructure.sisters.depthFilter.delta.y(:);
                    deltaXYZ(:,3) = intraStructure.sisters.depthFilter.delta.z(:);
                else
                    deltaXYZ(:,1) = intraStructure.sisters.raw.delta.x(:);
                    deltaXYZ(:,2) = intraStructure.sisters.raw.delta.y(:);
                    deltaXYZ(:,3) = intraStructure.sisters.raw.delta.z(:);
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
            case 'sisters'
                sisSep3D = intraStructure.sisters.sisSep.threeD(:);
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
                sisSep2D = intraStructure.plate.sisSep.twoD(:);
            case 'microscope'
                sisSep2D = intraStructure.microscope.sisSep.twoD(:);
            case 'sisters'
                sisSep2D = intraStructure.sisters.sisSep.twoD(:);
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
                sisSepXYZ(:,1) = intraStructure.plate.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.plate.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.plate.sisSep.z(:);
            case 'microscope'
                sisSepXYZ(:,1) = intraStructure.microscope.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.microscope.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.microscope.sisSep.z(:);
            case 'sisters'
                sisSepXYZ(:,1) = intraStructure.sisters.sisSep.x(:);
                sisSepXYZ(:,2) = intraStructure.sisters.sisSep.y(:);
                sisSepXYZ(:,3) = intraStructure.sisters.sisSep.z(:);
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
                twist3D = intraStructure.plate.twist.threeD(:);
            otherwise
                error('Twist is not defined for a non-plate coordinate system.')
        end
        
        median = nanmedian(twist3D);
        mean   = nanmean(twist3D);
        stdErr = nanserr(twist3D);
        stdDev = nanstd(twist3D);
        n      = min(sum(~isnan(twist3D)));

        fprintf('\n3D twist measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f deg\n',median)
        fprintf('    mean    = %.2f deg\n',mean)
        fprintf('    std err = %.2f deg\n',stdErr)
        fprintf('    std dev = %.2f deg\n',stdDev)
        
    case 'twistYZ'
        
        switch opts.coordSystem
            case 'plate'
                twistYZ(:,1) = intraStructure.plate.twist.y(:);
                twistYZ(:,2) = intraStructure.plate.twist.z(:);
            otherwise
                error('Twist is not defined for a non-plate coordinate system.')
        end

        median = nanmedian(twistYZ);
        mean   = nanmean(twistYZ);
        stdErr = nanserr(twistYZ);
        stdDev = nanstd(twistYZ);
        n      = min(sum(~isnan(twistYZ)));

        fprintf('\nY and Z twist measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] deg\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] deg\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] deg\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] deg\n',stdDev(1),stdDev(2))

    case 'swivel3D'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    swivel3D = intraStructure.plate.depthFilter.swivel.threeD(:);
                else
                    swivel3D = intraStructure.plate.raw.swivel.threeD(:);
                end
            case 'microscope'
                if opts.depthFilter
                    swivel3D = intraStructure.microscope.depthFilter.swivel.threeD(:);
                else
                    swivel3D = intraStructure.microscope.raw.swivel.threeD(:);
                end
            otherwise
                error('Swivel can only be defined for plate or microscope coordinate systems.')
        end
        
        median = nanmedian(swivel3D);
        mean   = nanmean(swivel3D);
        stdErr = nanserr(swivel3D);
        stdDev = nanstd(swivel3D);
        n      = min(sum(~isnan(swivel3D)));

        fprintf('\n3D swivel measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f deg\n',median)
        fprintf('    mean    = %.2f deg\n',mean)
        fprintf('    std err = %.2f deg\n',stdErr)
        fprintf('    std dev = %.2f deg\n',stdDev)
        
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
            otherwise
                error('Swivel can only be defined for plate or microscope coordinate systems.')
        end

        median = nanmedian(swivelKMT);
        mean   = nanmean(swivelKMT);
        stdErr = nanserr(swivelKMT);
        stdDev = nanstd(swivelKMT);
        n      = min(sum(~isnan(swivelKMT)));

        fprintf('\nKinetochore-to-microtubule angle measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f deg\n',median)
        fprintf('    mean    = %.2f deg\n',mean)
        fprintf('    std err = %.2f deg\n',stdErr)
        fprintf('    std dev = %.2f deg\n',stdDev)
        
    case 'swivelYZ'
        
        switch opts.coordSystem
            case 'plate'
                if opts.depthFilter
                    swivelYZ(:,1) = intraStructure.plate.depthFilter.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.plate.depthFilter.swivel.z(:);
                else
                    swivelYZ(:,1) = intraStructure.plate.raw.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.plate.raw.swivel.z(:);
                end
            case 'microscope'
                if opts.depthFilter
                    swivelYZ(:,1) = intraStructure.microscope.depthFilter.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.microscope.depthFilter.swivel.z(:);
                else
                    swivelYZ(:,1) = intraStructure.microscope.raw.swivel.y(:);
                    swivelYZ(:,2) = intraStructure.microscope.raw.swivel.z(:);
                end
            otherwise
                error('Swivel can only be defined for plate or microscope coordinate systems.')
        end
    
        median = nanmedian(swivelYZ);
        mean   = nanmean(swivelYZ);
        stdErr = nanserr(swivelYZ);
        stdDev = nanstd(swivelYZ);
        n      = min(sum(~isnan(swivelYZ)));

        fprintf('\nY and Z swivel measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] deg\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] deg\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] deg\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] deg\n',stdDev(1),stdDev(2))
        
    case 'rawInts'
        
        ints(:,1) = intraStructure.intensity.mean.inner(:);
        ints(:,2) = intraStructure.intensity.mean.outer(:);
        if opts.depthFilter
            filt = isnan(intraStructure.plate.depthFilter.delta.threeD(:));
            ints(repmat(filt,1,2)) = NaN;
        end
        
        median = nanmedian(ints)*1000;
        mean   = nanmean(ints)*1000;
        stdErr = nanserr(ints)*1000;
        stdDev = nanstd(ints)*1000;
        n      = min(sum(~isnan(ints)));

        fprintf('\nRaw inner- and outer-marker intensity measurements (n = %i):\n\n',n)
        fprintf('    median  = [%.2f, %.2f] au\n',median(1),median(2))
        fprintf('    mean    = [%.2f, %.2f] au\n',mean(1),mean(2))
        fprintf('    std err = [%.2f, %.2f] au\n',stdErr(1),stdErr(2))
        fprintf('    std dev = [%.2f, %.2f] au\n',stdDev(1),stdDev(2))
        
    case 'normInts'
        
        ints(:,1) = intraStructure.intensity.mean.inner(:);
        ints(:,2) = intraStructure.intensity.mean.outer(:);
        if opts.depthFilter
            filt = isnan(intraStructure.plate.depthFilter.delta.threeD(:));
            ints(repmat(filt,1,2)) = NaN;
        end
        normInts = ints;
        normInts(:,2) = normInts(:,2)./ints(:,1);
        normInts(:,1) = normInts(:,1)./ints(:,2);
        
        median = nanmedian(normInts);
        mean   = nanmean(normInts);
        stdErr = nanserr(normInts);
        stdDev = nanstd(normInts);
        n      = min(sum(~isnan(normInts)));

        fprintf('\nInner-normalised outer-marker intensity measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f au\n',median(2))
        fprintf('    mean    = %.2f au\n',mean(2))
        fprintf('    std err = %.2f au\n',stdErr(2))
        fprintf('    std dev = %.2f au\n',stdDev(2))
        fprintf('\n');
        fprintf('\nOuter-normalised inner-marker intensity measurements (n = %i):\n\n',n)
        fprintf('    median  = %.2f au\n',median(1))
        fprintf('    mean    = %.2f au\n',mean(1))
        fprintf('    std err = %.2f au\n',stdErr(1))
        fprintf('    std dev = %.2f au\n',stdDev(1))
        
    otherwise
        
        error('Statistic requested is not yet built into this version of dublBasicStats. See later release.')
    
end

fprintf('\n')

end