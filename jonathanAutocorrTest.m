function auotcorrs = jonathanAutocorrTest(jobset,ch)
  analysis = loadAnalysis(0,jobset,ch);
  fprintf('Computing autocorrelations');
  analysis = cuplAutocorrelation(analysis);
  fprintf('\n');
  autocorrs = analysis.autocorrs;

end
%end of main function

  function analysis = loadAnalysis(analysisLoaded,jobset,ch)
  persistent analysisP;
  if ~analysisLoaded || isempty(analysisP)
    analysisP.options.percentNan = 0.1;
    analysisP.options.minSistersPerCell = 1;
    analysisP.options.byDistNumBins = 12;
    analysisP.options.byDistMaxWidth = 12;
    analysisP.options.neighbourThreshold = 3;
    analysisP.options.monotelicAngle = 8;
    analysisP.options.trackChannel = 1;
    analysisP.options.channel = 1;
    analysisP.options.correctDrift = 0;
    analysisP.options.keepAllData = 1;
    analysisP.options.doTracks = 1;
    analysisP.options.poleCutoff = 4;
    analysisP.stages = {};

    fprintf('Loading data');
    analysisP = cuplLoadData(analysisP,jobset,ch);
    fprintf('Preprocessing data');
    analysisP = cuplPreprocess(analysisP);
    fprintf('\n');
    analysisLoaded = 1;
  end
  analysis = analysisP;
end