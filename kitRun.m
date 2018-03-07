function kitRun(desc)
% KITRUN Display GUI to enable running detection and analysis of spot intensities.
%
%  Run kitRun('help') for a more descriptive run screen.
%
% Copyright (c) 2017 C. A. Smith

% Check whether user has asked for help.
if nargin<1
  desc=0;
else
  desc=1;
end

% Setup GUI.
dockStatus = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal');
handles = createControls(desc);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);
set(0,'DefaultFigureWindowStyle',dockStatus);

%% NESTED FUNCTIONS
function hs = createControls(desc)
  
  % Create figure.
  figw = 70;
  figh = 21+desc*9;
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiT ' kitVersion(1)];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  % Define font sizes.
  largefont = 16;
  medfont = 14;
  smallfont = 12;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 0.5; %small horizontal shift
  toplabely = figh; %top-most point
  
  % Set up initial positions.
  x = dx;
  w = figw-2*dx;
  y = toplabely-lh;
  
  % Home console title.
  labw = w;
  t = label(hs.fig,'Welcome to Kinetochore Tracker',[x y labw h],largefont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Running tools sub-title.
  labw = w;
  t = label(hs.fig,'1. Making and running jobs',[x y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Brief explanation of use.
  if desc
  labw = w; nLines=2;
  y = y-(nLines-1)*h;
  label(hs.fig,'Make a new job for kinetochore tracking, then run it by pressing ''Run''. Alternatively, load and edit an existing job.',[x y labw nLines*h],smallfont);
  y = y-lh;
  end
  
  % Run buttons.
  btnw = 12; btnh = 2;
  hs.newBtn = button(hs.fig,'New run',[x y btnw btnh],@newRunCB);
  x = x+(btnw+ddx);
  hs.loadBtn = button(hs.fig,'Load run',[x y btnw btnh],@loadRunCB);
  y = y-h; x = dx;
  btnw = 10; labw = w-(2*btnw+dx+ddx);
  hs.jobsetLab = editbox(hs.fig,'',[x y labw h]);
  hs.jobsetLab.Enable = 'inactive';
  x = figw-(2*btnw+dx+ddx);
  hs.editBtn = button(hs.fig,'Edit',[x y btnw btnh],@editRunCB);
  x = figw-(btnw+dx);
  hs.runBtn = button(hs.fig,'Run',[x y btnw btnh],@executeCB);
  y = y-lh;
  x = dx; btnw = 12;
  hs.multirunBtn = button(hs.fig,'Multi-job run',[x y btnw btnh],@multiRunCB,10);
  y = y-lh;
  
  % Running tools sub-title.
  x = dx; labw = w;
  t = label(hs.fig,'2. Manually check data',[x y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Brief explanation of use.
  if desc
  labw = w; nLines=1;
  y = y-(nLines-1)*h;
  label(hs.fig,'Deselect erroneous spots in order to clean up analysis.',[x y labw nLines*h],smallfont);
  y = y-lh;
  end
  
  % Manual spot and cell filter buttons.
  btnw = 15; btnh = 2;
  hs.spotFilterBtn = button(hs.fig,'Spot filtering',[x y btnw btnh],@filterSpotsCB);
  x = x+btnw+ddx;
  hs.cellFilterBtn = button(hs.fig,'Cell filtering',[x y btnw btnh],@filterCellsCB);
  x = x+btnw+ddx;
  hs.sisterPairingBtn = button(hs.fig,'Sister pairing',[x y btnw btnh],@sisterPairingCB);
  x = x+btnw+ddx;
  hs.categoriesBtn = button(hs.fig,'Categorise spots',[x y btnw btnh],@categoriseCB);
  hs.categoriesBtn.Enable = 'off';
  y = y-lh;
  
  % Analysis tools sub-title.
  x = dx; labw = w;
  t = label(hs.fig,'3. Analysis of results',[x y labw h],medfont);
  t.FontWeight = 'bold';
  y = y-lh;
  
  % Brief explanation of use.
  if desc
  labw = w; nLines=2;
  y = y-(nLines-1)*h;
  label(hs.fig,'Once one or more jobs have been run, run a full analysis. Multiple jobs can be loaded during analysis.',[x y labw nLines*h],smallfont);
  y = y-lh;
  end
  
  % CupL button.
  btnw = 15;
  hs.cuplBtn = button(hs.fig,'CupL analysis',[x y btnw btnh],@cuplCB);
  x = x+(btnw+ddx);
  hs.dublBtn = button(hs.fig,'DubL analysis',[x y btnw btnh],@dublCB);
  hs.dublBtn.Enable = 'off';
  x = x+(btnw+ddx);
  hs.intensityBtn = button(hs.fig,'Intensity analysis',[x y btnw btnh],@intensityCB);
  hs.intensityBtn.Enable = 'off';
  
  % Close button.
  x = figw-btnw-dx; y = 1;
  hs.closeBtn = button(hs.fig,'Close',[x y btnw btnh],@closeCB);
  
  % KiT logo. Match its position and dimensions to the save button.
  logow = 12.5; logox = figw-(logow+dx); y=toplabely-(logow+h)/3;
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[logox y logow logow/3]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kitlogo.png'),pos([4 3])));

  movegui(hs.fig,'center');
  
end

function newRunCB(hObj,event)
  kitLog('Opening jobset setup window')
  hs = handles;
  % run GUI to make new jobset
  hs.jobset = [];
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  % Small GUI to ask user which type of job to set up.
  mode = chooseMode;
  
  jobset = kitSetupJob(mode);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
  % Check whether user cancelled setup.
  if isfield(jobset,'cancel')
    return
  end
  
  % push jobset name to label
  filename = jobset.filename;
  maxMovLen = 55;
  filename = strshorten(filename,maxMovLen);
  hs.jobsetLab.String = filename;
  hs.runBtn.String = 'Run';
  % back up all information
  hs.jobset = jobset;
  handles = hs;
end

function loadRunCB(hObj,event)
  kitLog('Loading jobset')
  hs = handles;
  % run GUI to make new jobset
  hs.jobset = [];
  [filename,pathname] = uigetfile('*.mat','Select a KiT jobset file');
  jobset = kitLoadJobset(fullfile(pathname,filename));
  % push jobset name to label
  filename = jobset.filename;
  maxMovLen = 55;
  filename = strshorten(filename,maxMovLen);
  hs.jobsetLab.String = filename;
  hs.runBtn.String = 'Run';
  % back up all information
  hs.jobset = jobset;
  handles = hs;
end

function editRunCB(hObj,event)
  hs = handles;
  % check if a jobset is loaded
  if ~isfield(hs,'jobset') || isempty(hs.jobset)
    errorbox('No job yet created or loaded.')
    return
  end
  kitLog('Opening jobset setup window with loaded job')
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  hs.jobset = kitJobset(hs.jobset);
  jobset = kitSetupJob(hs.jobset.options.jobProcess,hs.jobset);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
  % push jobset name to label
  filename = jobset.filename;
  maxMovLen = 55;
  filename = strshorten(filename,maxMovLen);
  hs.jobsetLab.String = filename;
  hs.runBtn.String = 'Run';
  % back up all information
  hs.jobset = jobset;
  handles = hs;
end

function multiRunCB(hObj,event)
  
  % Get the currently-loaded jobset if exists.
  if isfield(handles,'jobset')
    jobset = handles.jobset;
  else
    jobset = [];
  end
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  kitMultiJobsetRun(jobset);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');

end

function filterSpotsCB(hObj,event)
  hs = handles;
  % check if a jobset is loaded
  if ~isfield(hs,'jobset') || isempty(hs.jobset)
    errorbox('No job yet created or loaded.')
    return
  end
  kitLog('Starting manual spot filtering')
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  kitFilterSpots(hs.jobset);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
end

function filterCellsCB(hObj,event)
  hs = handles;
  % check if a jobset is loaded
  if ~isfield(hs,'jobset') || isempty(hs.jobset)
    errorbox('No job yet created or loaded.')
    return
  end
  kitLog('Starting manual cell filtering')
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  kitFilterCells(hs.jobset);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
end

function sisterPairingCB(hObj,event)
  hs = handles;
  % check if a jobset is loaded
  if ~isfield(hs,'jobset') || isempty(hs.jobset)
    errorbox('No job yet created or loaded.')
    return
  end
  kitLog('Starting manual pairing of sister kinetochores')
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  kitManualPairSisters(hs.jobset);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
end

function executeCB(hObj,event)
  if ~isfield(handles,'jobset') || isempty(handles.jobset)
    errorbox('No job yet created or loaded.')
    return
  end
  
  % Turn off the GUI during processing.
  handles.runBtn.String = 'Running...';
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  pause(1);
  
  % Run the job.
  kitRunJob(handles.jobset,'tasks',handles.jobset.options.tasks);
  
  % We turn back on the interface
  set(guiObj,'Enable','on');
  handles.runBtn.String = 'Re-run';
end

function cuplCB(hObj,event)
  kitLog('Starting CupL analysis')
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  cuplAnalysis(handles.jobset);
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
end

function dublCB(hObj,event)
  kitLog('Starting DubL analysis')
  
  % Turn off the GUI during processing.
  guiObj=findobj(handles.fig,'Enable','on');
  set(guiObj,'Enable','inactive');
  
  dublAnalysis;
  
  % Re-activate GUI.
  set(guiObj,'Enable','on');
  
end

function closeCB(hObj,event)
  close(gcf);
end

%% Other functions

function mode = chooseMode()
  
  % List all modes.
  allModes = {'Tracking','Single timepoint','Chromatic shift'};
  allModesJS = {'zandt','zonly','chrshift'};
    
  % Create figure.
  figw = 50;
  figh = 5;
  f = figure('Resize','off','Units','characters','Position',[100 35 figw figh]);
  f.DockControls = 'off';
  f.MenuBar = 'none';
  f.Name = 'Choose a job type';
  f.NumberTitle = 'off';
  f.IntegerHandle = 'off';
  f.ToolBar = 'none';
  
  % Define font sizes.
  smallfont = 12;
  
  % Set some standard positions and distances.
  h = 1.5; %height
  lh = 1.5*h; %large height
  dx = 2.5; %horizontal shift
  ddx = 0.5; %small horizontal shift
  toplabely = figh; %top-most point
  
  % Set up initial positions.
  x = dx;
  w = figw-2*dx;
  y = toplabely-lh;
  
  % Print choices.
  labw = w;
  label(f,'Choose which type of job you wish to set up:',[x y labw h],smallfont);
  y = y-lh;
  btnw = (w-2*ddx)/3; btnh = 2;
  for i=1:3
    button(f,allModes{i},[x y btnw btnh],@chooseModeCB,smallfont);
    x = x+(btnw+ddx);
  end
  
  movegui(f,'center');
  uiwait(f);
  
  % Get the mode for the jobset.
  mode = allModesJS{mode};
    
  function chooseModeCB(hObj,event)
    mode = find(cellfun(@(x) strcmp(x,hObj.String),allModes));
    uiresume(f);
    close(f);
  end
  
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

end % kitSetupJob
