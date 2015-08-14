function kitDownloadBioFormats()
% KITDOWNLOADBIOFORMATS Download BioFormats, if required.
%
%     Download BioFormats, if required. Requires current directory to be KiT
%     directory.
%
% Copyright (c) 2015 Jonathan W. Armond

filename = 'bfmatlab';
url = ['http://downloads.openmicroscopy.org/bio-formats/5.0.6/artifacts/' ...
       'bfmatlab.zip'];

% Check for existing jar.
download = 0;
if exist(filename,'dir') == 0
  % Directory doesn't exist, check for zip.

  if exist([filename '.zip'],'file') == 0
    % Zip doesn't exist. Download it.
    kitLog('Downloading BioFormats...');
    % Check we are in KiT directory.
    if exist(fullfile(pwd,'kitDownloadBioFormats.m'),'file') ~= 2
      error('Must be in KiT directory to download Loci Tools');
    end

    unzip(url);
    kitLog('Download complete.');
  else
    % Unzip.
    unzip([filename '.zip']);
  end
end
