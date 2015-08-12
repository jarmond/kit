function kitJobsetDiagnostics(jobset,channel,short,fid)
% KITJOBSETDIAGNOSTICS Prints a summary of tracking results for all jobs in jobset.
%
% Copyright (c) 2013 Jonathan W. Armond

if nargin<2
  channel = 1;
end
if nargin<3
  short = 1;
end
if nargin<4
  fid = 1;
end

fields = {'nSpotsPerFrame','nTracks','nSisters','avgSisterTrackLength', ...
          'percentWithPlane','sisterVar','nLongSisters'};
fieldNames = {'job','# spots','# tracks','# sisters','sister length','plane fits',...
             'variance','# long sisters'};
fieldFmt = {'d','.1f','d','d','.1f','d','.4f','d'};
if short == 0
  fields = [fields {'elapsedTime'}];
  fieldNames = [fieldNames {'cputime'}];
  fieldFmt = [fieldFmt {'.0f'}];
end

% Read in diagnostics from each job.
nFields = length(fields);
nJobs = length(jobset.movieFiles);
stats = nan(nJobs,nFields+1);
for i=1:nJobs
  try
    job = kitLoadJob(jobset,i);
  catch
    continue;
  end
  if ~isfield(job,'dataStruct')
    continue
  end

  ds = job.dataStruct{channel};
  if isempty(ds)
    continue
  end

  if ~isfield(ds,'diagnostics')
    continue
  end
  diag = ds.diagnostics;
  stats(i,1) = i;
  for f=1:nFields
    field = fields{f};
    if isfield(diag,field)
      stats(i,f+1) = diag.(field);
    else
      warning('Missing field: %s',field);
    end
  end
end

% Print out diagnostics table.
printHeader(strjoin(fieldNames,'|'),fid);
l = cellfun(@(x) length(x), fieldNames);
for i=1:size(stats,1)
  for j=1:size(stats,2)
    fmt = sprintf('%%%d%s ',l(j),fieldFmt{j});
    fprintf(fid,fmt,stats(i,j));
  end
  fprintf(fid,'\n');
end

% Print totals
fprintf(fid,'sum ');
for j=2:size(stats,2)
  fmt = sprintf('%%%d.1f ',l(j));
  fprintf(fid,fmt,nansum(stats(:,j)));
end

% Print averages
fprintf(fid,'avg ');
for j=2:size(stats,2)
  fmt = sprintf('%%%d.1f ',l(j));
  fprintf(fid,fmt,nanmean(stats(:,j)));
end
fprintf(fid,'\n');
