function fmtregexp=kitSupportedFormats(asfilterspec)
% KITSUPPORTEDFORMATS Regular expression matching all supported formats.
%
%    FMTREGEXP = KITSUPPORTEDFORMATS() Returns regular expression matching all
%    supported formats.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<1
  asfilterspec = 0;
end

% TODO use bfGetFileExtensions

% Support DV and OME-TIFFs.
if asfilterspec == 0
  fmtregexp = '(dv|ome.tiff|ome.tif|r3d|d3d)$';
else
  fmtregexp = {'*.*', 'All files';
               '*.dv;*.r3d;*.d3d', 'DeltaVision files';...
               '*.ome.tiff;*.ome.tif', 'OME files'};
end