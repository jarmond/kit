function filters=createFilters(nDims,dataProperties)

backFilterParms = dataProperties.FT_SIGMA * 14; % 14 is 15-1
backFilterParms(4:6) = roundOddOrEven(backFilterParms,'odd','inf');
% Constrain filter by number of planes in stack.
backFilterParms(6) = min(backFilterParms(6), (dataProperties.movieSize(3)-1)*2+1);

switch nDims
  case 3
    backgroundFilter = GaussMask3D(backFilterParms(1:3),...
        backFilterParms(4:6),[],1,[],[],1);
    signalFilter = GaussMask3D(dataProperties.FILTERPRM(1:3),...
        dataProperties.FILTERPRM(4:6),[],1,[],[],1);
    signalFilterP = dataProperties.FILTERPRM(4:6);
    backFilterP = backFilterParms(4:6);
    borderMode = 2;
    % make separated noise mask
    noiseMask = {...
        ones(dataProperties.FILTERPRM(4),1,1)./dataProperties.FILTERPRM(4),...
        ones(1,dataProperties.FILTERPRM(5),1)./dataProperties.FILTERPRM(5),...
        ones(1,1,dataProperties.FILTERPRM(6))./dataProperties.FILTERPRM(6),...
        };
  case 2
    backgroundFilter = GaussMask2D(backFilterParms(1:2),...
        backFilterParms(4:5),[],1,[]);
    signalFilter = GaussMask2D(dataProperties.FILTERPRM(1:2),...
        dataProperties.FILTERPRM(4:5),[],1,[]);
    signalFilterP = dataProperties.FILTERPRM(4:5);
    backFilterP = backFilterParms(4:5);
    borderMode = 1;
    % make separated noise mask
    noiseMask = {...
        ones(dataProperties.FILTERPRM(4),1,1)./dataProperties.FILTERPRM(4),...
        ones(1,dataProperties.FILTERPRM(5),1)./dataProperties.FILTERPRM(5)};
  otherwise
    error('Unsupported dimensions: %d',ndims);
end

filters.background = backgroundFilter;
filters.signal = signalFilter;
filters.backgroundP = backFilterP;
filters.signalP = signalFilterP;
filters.border = borderMode;
filters.noise = noiseMask;