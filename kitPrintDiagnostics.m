function kitPrintDiagnostics(dataStruct,process)
% KITPRINTDIAGNOSTICS Print out basic tracking diagnostics
%
% Created by: Jonathan W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2018 C. A. Smith

diag = dataStruct.diagnostics;

fprintf('Elapsed compute time: %s\n',timeString(diag.elapsedTime));
if diag.failed
  fprintf('Job failed');
else
  fprintf('Particles per frame: %.1f\n',diag.nSpotsPerFrame);
  if strcmp(process,'zandt')
    error('Option zandt not available for KiD. Try using a later version of KiT for tracking'); 
  end
  if ~strcmp(process,'chrshift')
    fprintf('Frames with a plane fit: %.2f%%\n',diag.percentWithPlane);
  end
end
