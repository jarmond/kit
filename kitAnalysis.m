function kitAnalysis(jobset,ch)
% KITBASICPLOTS Produces a set of basic plots for a job
%
%   kitAnalysis(jobset,channel)
%
%  jobset: Struct containing tracking job setup options.
%
%  channel: Channel to plot. Defaults to 1.
%
% Copyright (c) 2015 Jonathan W. Armond

if nargin<2
  ch = 1;
end

% Setup GUI.
roiChanged = 1;
handles = createControls(jobset);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls(jobset)
  w = 30;
  h = 50;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 w+2 h]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiT ' kitVersion(1)];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';
  movegui(hs.fig,'center');
  figpos = hs.fig.Position;

  x = 1;
  y = h - 2;
  h = 2;
  lh = 1.5; % label height
  lw = 0.75*w; % label width
  ex = lw + 1; % edit box start
  ew = 0.25*w - 1; % edit box width

  [~,name] = fileparts(jobset.filename);
  t=label(hs.fig,['Jobset: ' name],[x y w lh]);
  t.FontWeight='bold';
  y=y-h;
  t=label(hs.fig,'Per cell analysis',[x y w lh]);
  t.FontWeight='bold';

  % Get ROI names.
  maxMovLen = 32;
  movieFiles = jobset.movieFiles;
  ROIString = {};
  if ~isempty(movieFiles)
    for i=1:length(jobset.ROI)
      ROIString{i} = [strshorten(movieFiles{jobset.ROI(i).movieIdx},maxMovLen) ' [' ...
                          num2str(round(jobset.ROI(i).crop),'%d ') ']'];
    end
  end
  y=y-h;
  label(hs.fig,'ROIs:',[x y w lh]);
  y=y-h;
  hs.ROI = popup(hs.fig,ROIString,[x y w h],@ROICB);

  y=y-h;
  label(hs.fig,'Min length %',[x y w lh]);
  hs.minLength = editbox(hs.fig,'75',[ex y ew h]);
  y=y-h;
  hs.sisters = button(hs.fig,'Sister trajectories',[x y w h],@sistersCB);
  y=y-h;
  hs.plate = button(hs.fig,'Metaphase plate',[x y w h],@plateCB);

  y=y-2*h;
  t=label(hs.fig,'Experiment analysis',[x y w lh]);
  t.FontWeight='bold';
  y=y-h;
  hs.diags = button(hs.fig,'Print diagnostics',[x y w h],@diagsCB);
end

function ROICB(hObj,event)
  roiChanged = 1;
end

function sistersCB(hObj,event)
  job = loadActiveJob();
  kitPlotSisters(job,'minLength',str2double(handles.minLength.String)/100,'channel',ch);
end

function plateCB(hObj,event)
  job = loadActiveJob();
  kitPlateCheck(job,'channel',ch);
end

function diagsCB(hObj,event)
  kitJobsetDiagnostics(jobset,ch);
end

function job = loadActiveJob()
  persistent jobP
  if roiChanged
    idx = handles.ROI.Value;
    jobP = kitLoadJob(jobset,idx);
    roiChanged = 0;
  end
  job = jobP;
end

end

%% LOCAL FUNCTIONS

function m=expand(m,rows)
% Expand m to have at least rows.
  sizeDiff = nFrames - size(m,1);
if ~isempty(m) && sizeDiff > 0
  m = [m; nan(sizeDiff,size(m,2))];
end
end
