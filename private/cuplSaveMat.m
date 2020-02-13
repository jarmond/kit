function cuplSaveMat(analysis)
%CUPLSAVEMAT  Save analysis structure to mat-file.
%
%   CUPLSAVEMAT(ANALYSIS) Save analysis structure to mat-file.
%
% Copyright (c) 2010 Elina Vladimirou
% Copyright (c) 2013 Jonathan Armond

save(fullfile(analysis.outputDirectory,analysis.outputFilename),...
     'analysis','-v7.3');
