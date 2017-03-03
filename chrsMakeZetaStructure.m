function intraStruct = chrsMakeZetaStructure()
% CHRSMAKEZETASTRUCTURE Produces an empty structure into which
% inter-marker measurements are compiled using chrsZetaMeasurements.
%
%    CHRSMAKEZETASTRUCTURE() Produces a structure allowing for inter-marker
%    measurements to be compiled in order to draw population-scale
%    analyses. No input is required.
%
%
% Copyright (c) 2017 C. A. Smith

intraStruct.dublVersion = 'ChrS 1.0.0';
intraStruct.channelVector = [1 2];

coords = struct('x',[],'y',[],'z',[]);
zeta  = struct('x',[],'y',[],'z',[],'twoD',[],'threeD',[]);

intraStruct.coords = coords;
intraStruct.zeta = zeta;
    
    
    
    