function [metadata,reader]=kitOpenMovie(movieFileName)
% KITOPENMOVIE Open movie and extract metadata from BioFormats metadata
%
%    [METADATA,READER] = KITOPENMOVIE(MOVIEFILENAME) Open movie with filename
%    MOVIEFILENAME and extract METADATA from BioFormats metadata. Returns a
%    READER Java object which can be used to extract image data.
%
% Copyright (c) 2013 Jonathan W. Armond

kitLog('Opening movie: %s', movieFileName);
if ~exist(movieFileName,'file')
  error('Could not find file: %s',movieFileName);
end

addpath bfmatlab;
bfCheckJavaPath(1);

reader = bfGetReader(movieFileName);

% Read basic image metadata.
md.nFrames = reader.getSizeT();
md.nChannels = reader.getSizeC();
md.nPlanes = reader.getImageCount();
pixelType = reader.getPixelType();
md.nBytesPerPixel = loci.formats.FormatTools.getBytesPerPixel(pixelType);
md.isSigned = loci.formats.FormatTools.isSigned(pixelType);
md.isFloatingPoint = loci.formats.FormatTools.isFloatingPoint(pixelType);
md.isLittleEndian = reader.isLittleEndian();

switch md.nBytesPerPixel
  case 1
    if md.isSigned, md.dataType = 'int8'; else md.dataType = 'uint8'; end
  case 2
    if md.isSigned, md.dataType = 'int16'; else md.dataType = 'uint16'; end
  case 4
    md.dataType = 'single';
  case 8
    md.dataType = 'double';
  otherwise
    error('Only 8- or 16-bit integer, or 32- or 64-bit floating point data supported');
end

md.frameSize = [reader.getSizeX() reader.getSizeY() reader.getSizeZ()];


% Read additional metadata
metaTable = reader.getMetadataStore();
try

% Wavelengths. Needed to estimate PSFs.
numWvs = metaTable.getChannelCount(0);
md.wavelength = [507 610 461]/1000; % Default assumes EGFP, mCherry,
                                    % DAPI. FIXME Ask user.
warnWv = 0;
for i=1:numWvs
  chWv = metaTable.getChannelEmissionWavelength(0,i-1);
  if ~isempty(chWv)
    md.wavelength(i) = chWv.getValue()/1000;
  else
    warnWv = 1;
  end
end
if warnWv
  warning('Missing metadata: Assuming wavelengths %d, %d, %d nm',...
          1000*md.wavelength);
end

% Timepoints per plane
nZPlanes = md.frameSize(3);
nTimepoints = md.nFrames;
idx = 0;
defDt = 2; % Default assume every 2 sec. FIXME Ask user.
warnT = 0;
for i=1:nTimepoints
  defT = (i-1)*defDt;
  for j=1:nZPlanes
    try
      md.frameTime(j,i) = metaTable.getPlaneDeltaT(0, idx).doubleValue;
    catch
      % Use default, if missing metadata.
      md.frameTime(j,i) = defT;
      warnT = 1;
    end
    idx = idx+1;
  end
end
if warnT
  warning('Missing metadata: Assuming %d s per frame.', ...
          defDt);
end
dT = diff(md.frameTime(1,:));
warndT = 0.5;
if max(dT)-min(dT) > warndT;
  warning('Frame time deviates by more than %f',warndT);
end

% Is 3D image?
md.is3D = nZPlanes > 1;

% Numerical aperture
try
  md.na = metaTable.getObjectiveLensNA(0,0).doubleValue;
catch
  warning('Missing metadata: Assuming NA = 1.4');
  md.na = 1.4; % Assume default.
end

% Physical pixel size
md.pixelSize = [
  metaTable.getPixelsPhysicalSizeX(0).getValue(),...
  metaTable.getPixelsPhysicalSizeY(0).getValue(),...
  metaTable.getPixelsPhysicalSizeZ(0).getValue()];
if any(md.pixelSize < 0.001) || any(md.pixelSize(1:2) > 1) || ...
    md.pixelSize(3) > 5
  warning('Pixel sizes are strange: %f x %f x %f',md.pixelSize);
end


catch e
  reader.close();
  clear reader;
  error(['Missing required metadata: ' e.message]);
end

% Rename output.
metadata = md;
