function intraStruct = dublMakeIntraStructure(paired)
% DUBLMAKEINTRASTRUCTURE Produces an empty structure into which
% intra-kinetochore measurements are compiled using dublIntraMeasurements.
%
%    DUBLMAKEINTRASTRUCTURE() Produces a structure allowing for both inter-
%    and intra-kinetochore measurements to be compiled in order to draw
%    population-scale analyses. No input is required.
%
%
% Copyright (c) 2016 C. A. Smith

if nargin<1 || isempty(paired)
  paired = 1;
end

intraStruct.dublVersion = 'DubL 1.1.0';
intraStruct.exptLabel = [];
intraStruct.options.channels = [];
intraStruct.options.paired = paired;
intraStruct.options.plate = 0;

% make substructures for all pair-derived measurements
direction = struct('P',[],'AP',[],'N',[],'S',[]);
sisSep = struct('x',[],'y',[],'z',[],...
               'twoD',[],'threeD',[]);
subint = struct('inner',[],'outer',[]);
intensity = struct('mean',subint,'max',subint,'bg',subint);

% make substructure for each microscope and plate coordinate systems
coords = struct('x',[],'y',[],'z',[]);
delta  = struct('x',[],'y',[],'z',[],...
               'oneD',[],'twoD',[],'threeD',[]);
swivel = struct('y',[],'z',[],'threeD',[],'kMT',[]);
% those only appropriate for plate coordinate system
twist  = struct('y',[],'z',[],'threeD',[]);
sisterCentreSpeed = [];
plateThickness = [];

% produce generic substructure for microscope coordinate system
microscope.coords             = coords;
if paired
microscope.sisSep             = sisSep;
end
  microscope.raw.delta          = delta;
  microscope.depthFilter.delta  = delta;
if paired
  microscope.raw.swivel         = swivel;
  microscope.depthFilter.swivel = swivel;
end
% produce same again for plate coordinate system
plate = microscope;
if paired
  plate.raw.twist         = twist;
  plate.depthFilter.twist = twist;
  plate.sisterCentreSpeed = sisterCentreSpeed;
  plate.plateThickness    = plateThickness;
  
    intraStruct.direction = direction;
end
    intraStruct.microscope = microscope;
    intraStruct.plate = plate;
    intraStruct.sisters = microscope;
    intraStruct.intensity = intensity;
    
end
    
    
    
    