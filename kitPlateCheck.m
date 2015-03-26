function kitPlateCheck(job,channel)
% KITPLATECHECK Visualize plate in xy-plane for coord-sys verification
%
%    KITPLATECHECK(JOB,CHANNEL) Visualize plate in xy-plane for
%    coordinate-system verification. JOB may be either a job struct, loaded via
%    kitLoadJob, or a sisterList.
%
%    CHANNEL Optional, specify tracked channel to display.

if nargin<2
  channel = 1;
end

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
plot(coords1(:,1:6:end),coords1(:,2:6:end),'b-',...
     coords2(:,1:6:end),coords2(:,2:6:end),'g-');
xlabel('x');
ylabel('y');
