function kitPlateCheck(job,varargin)
% KITPLATECHECK Visualize plate in xy-plane for coord-sys verification
%
%    KITPLATECHECK(JOB,...) Visualize plate in xy-plane for
%    coordinate-system verification. JOB may be either a job struct, loaded via
%    kitLoadJob, or a sisterList.
%
%    'channel': Optional, specify tracked channel to display.
%
%    'colors': Optional, set to 1 to show multiple colors to distinguish
%    tracks. Otherwise uses colors to distinguish sisters

opts.channel = 1;
opts.colors = 0;
opts = processOptions(opts,varargin{:});

if isfield(job,'dataStruct')
  % Is a job struct.
  sisterList = job.dataStruct{channel}.sisterList;
else
  % Assume a sisterList.
  sisterList = job;
end

coords1 = horzcat(sisterList.coords1);
coords2 = horzcat(sisterList.coords2);
figure;
if opts.colors
  plot(coords1(:,1:6:end),coords1(:,2:6:end),...
       coords2(:,1:6:end),coords2(:,2:6:end));
else
  plot(coords1(:,1:6:end),coords1(:,2:6:end),'b-',...
       coords2(:,1:6:end),coords2(:,2:6:end),'g-');
end
xlabel('x');
ylabel('y');
