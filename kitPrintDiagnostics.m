function kitPrintDiagnostics(dataStruct)
% KITPRINTDIAGNOSTICS Print out basic tracking diagnostics
%
% Copyright (c) 2013 Jonathan W. Armond

diag = dataStruct.diagnostics;

fprintf('Elasped compute time: %s\n',timeString(diag.elapsedTime));
fprintf('Spots per frame: %.1f\n',diag.nSpotsPerFrame);
fprintf('# sister pairs tracked: %d\n',diag.nSisters);
fprintf('# individual KTs tracked: %d\n',diag.nTracks);
fprintf('Average sister pair track length: %.1f\n',diag.avgSisterTrackLength);
fprintf('Average KT track length: %.1f\n',diag.avgTrackLength);
fprintf('# long sisters (75%% length): %d\n',diag.nLongSisters);
fprintf('Frames with a plane fit: %.2f%%\n',diag.percentWithPlane);
