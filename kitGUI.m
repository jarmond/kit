function jobset=kitGUI(jobset)
% KITGUI Display GUI to setup and run tracking.
%
% Copyright (c) 2015 Jonathan W. Armond

% Download BioFormats, if required.
kitDownloadBioFormats();

if nargin<1 || isempty(jobset)
  jobset = kitDefaultOptions();
end

if ~isfield(jobset,'ROI')
  jobset.ROI = [];
end

coordSystemValues = {'Plate','Image'}; % in GUI
coordSystemValuesJS = lower(coordSystemValues); % in jobset
spotDetectValues = {'Histogram','Adaptive','None'};
spotDetectValuesJS = {'histcut','adaptive','none'};
spotRefineValues = {'Centroid','MMF','None'};
spotRefineValuesJS = {'centroid','gaussian','none'};

% Setup GUI. TODO feed in default jobset options
handles = createControls(jobset);
populateROIBox();
spotModeCB();
refineModeCB();
autoRadiiCB();
handles.fig.Visible = 'on';

uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls(jobset)
  opts = jobset.options;

  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 138 44]);
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiT ' kitVersion(1)];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  movegui(hs.fig,'center');
  figpos = hs.fig.Position;
  colwidth = [53 42 35];

  w=25; h=8;
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[figpos(3)-w-3 1 w h]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kitlogo.png'),pos([4 3])));

  %% ROI selection
  w = colwidth(1);
  toplabely = figpos(4)-2;
  % Movies
  hs.movies = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[2.5 23 w 15.15],'Max',inf,'Min',0);
  hs.movies.String = jobset.movieFiles;
  hs.movieDirectory = editbox(hs.fig,'',[2.5 40.1 w 1.7]);
  hs.movieDirectory.Enable = 'inactive';
  hs.movieDirectory.String = jobset.movieDirectory;
  hs.selectDirectory = button(hs.fig,'Select directory',[38.8 38.1 17.5 2],@selectDirectoryCB);
  label(hs.fig,'Available movies:',[2.5 38.3 20 1.5]);
  label(hs.fig,'ROI selection:',[2.5 toplabely 20 1.5],14);

  % ROIs
  hs.ROIs = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[2.5 3 w 18.15]);
  label(hs.fig,'ROIs:',[2.5 21.2 20 1.5]);
  hs.addROI = button(hs.fig,'Add ROI',[2.5 1 13 2],@addROICB);
  hs.viewROI = button(hs.fig,'View ROI',[17 1 13 2],@viewROICB);
  hs.deleteROI = button(hs.fig,'Delete ROI',[32 1 13 2],@deleteROICB);

  %% Tracking setup
  x = colwidth(1) + 3;
  w = colwidth(2);
  t = label(hs.fig,'Tracking setup',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  % Coordinate system
  label(hs.fig,'Coordinate system',[x 40 20 1.5]);
  hs.coordSys = popup(hs.fig,coordSystemValues,[74 40.1 20 1.5]);
  hs.coordSys.Value = mapStrings(opts.coordSystem,coordSystemValuesJS);
  label(hs.fig,'Coordinate system channel',[x 38 30 1.5]);
  hs.coordSysCh = editbox(hs.fig,num2str(opts.coordSystemChannel),[87 38 6 1.5]);

  % Channel modes
  b = 32;
  h = 5.5;
  for i=1:3
    p = uipanel(hs.fig,'Units','characters','Position',[x b-(i-1)*h w h],'FontSize',12,'Title',['Channel ' num2str(i)]);
    hs.spotMode{i} = popup(p,spotDetectValues,[22 2.5 0.45*w 1.5],@spotModeCB);
    hs.spotMode{i}.Value = mapStrings(opts.spotMode{i},spotDetectValuesJS);
    label(p,'Spot detection',[1 2.5 20 1.5]);

    hs.refineMode{i} = popup(p,spotRefineValues,[22 0.5 0.45*w 1.5],@refineModeCB);
    hs.refineMode{i}.Value = mapStrings(opts.coordMode{i},spotRefineValuesJS);
    label(p,'Spot refinement',[1 0.5 20 1.5]);
  end

  b = b-2*h;
  % Tasks
  % tasks = {'Spot finding','Plane fitting','Tracking','Sister grouping','Intensity measurement'};
  % h = 1.75;
  % panelh = h*(1+length(tasks));
  % p = uipanel(hs.fig,'Units','characters','Position',[x b-panelh w panelh],'FontSize',12,'Title','Tasks');
  % for i=1:length(tasks)
  %    hs.tasks{i} = checkbox(p,tasks{i},[1 panelh-h*(i+1) 30 h]);
  %    hs.tasks{i}.Value = hs.tasks{i}.Max;
  % end
  % b = b-panelh;

  %% Execution
  y = b-h;
  labelw = 0.5*w;
  t = label(hs.fig,'Execution',[x y labelw 1.5],14);
  t.FontWeight = 'bold';
  h = 2;
  y = y-h;
  label(hs.fig,'Jobset name',[x y labelw h],12);
  hs.filename = editbox(hs.fig,'jobset',[x+w-(w-labelw) y (w-labelw) h]);
  btnw = 0.5*w;
  bx = x + w - btnw;
  y = y-2*h;
  hs.save = button(hs.fig,'Save',[bx y btnw h],@saveCB);
  y = y-h;
  hs.execute = button(hs.fig,'Execute',[bx y btnw h],@executeCB);


  %% Options
  x = sum(colwidth(1:2)) + 4;
  w = 35;
  t = label(hs.fig,'Tracking options',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  h = 1.5;
  lh = 1.5*h; % large height
  y = toplabely-h;
  hs.autoRadii = checkbox(hs.fig,'Calculate search radii from dt',[x y w h],@autoRadiiCB,10);
  y = y-h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  t = label(hs.fig,'Frame dt',[x y labelw h],10);
  hs.autoRadiidt = editbox(hs.fig,num2str(opts.autoRadiidt),[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Est. avg. disp. of spots (um)',[x y labelw h],10);
  % Assume mean absolute displacment of sisters is about 0.06 μm/s as default.
  hs.autoRadiiAvgDisp = editbox(hs.fig,num2str(0.06),[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Min search radius (um)',[x y labelw h],10);
  hs.minSearchRadius = editbox(hs.fig,num2str(opts.minSearchRadius(1)),[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Max search radius (um)',[x y labelw h],10);
  hs.maxSearchRadius = editbox(hs.fig,num2str(opts.maxSearchRadius(1)),[editx y editw h],10);
  if isempty(opts.autoRadiidt)
    hs.autoRadii.Value = hs.autoRadii.Min; % Off
    hs.autoRadiidt.Enable = 'off';
    hs.autoRadiiAvgDisp.Enable = 'off';
    hs.minSearchRadius.Enable = 'on';
    hs.maxSearchRadius.Enable = 'on';
  else
    hs.autoRadii.Value = hs.autoRadii.Max; % On
    hs.autoRadiidt.Enable = 'on';
    hs.autoRadiiAvgDisp.Enable = 'on';
    hs.minSearchRadius.Enable = 'off';
    hs.maxSearchRadius.Enable = 'off';
  end

  y = y-h;
  hs.useSisterAlignment = checkbox(hs.fig,'Use sister alignment',[x y w h],@useSisterAlignmentCB,10);
  hs.useSisterAlignment.Value = opts.useSisterAlignment;
  y = y-h;
  % Adjust text box pos for multiple lines.
  t = label(hs.fig,'Max angle between sisters and plate normal (deg)',[x y-h/2 labelw lh],10);
  hs.maxSisterAlignmentAngle = editbox(hs.fig,num2str(opts.maxSisterAlignmentAngle),[editx y editw h],10);
  if ~hs.useSisterAlignment.Value
    hs.maxSisterAlignmentAngle.Enable = 'off';
  end
  y = y-lh;
  t = label(hs.fig,'Max average distance between sisters (um)',[x y-h/2 labelw lh],10);
  hs.maxSisterDist = editbox(hs.fig,num2str(opts.maxSisterSeparation),[editx y editw h],10);
  y = y-lh;
  t = label(hs.fig,'Min overlap between sister tracks',[x y-h/2 labelw lh],10);
  hs.minSisterTrackOverlap = editbox(hs.fig,num2str(opts.minSisterTrackOverlap),[editx y editw h],10);
  y = y-lh;
  t = label(hs.fig,'Min spots per frame',[x y labelw h],10);
  hs.minSpotsPerFrame = editbox(hs.fig,num2str(opts.minSpotsPerFrame),[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Weight for spot count',[x y labelw h],10);
  hs.adaptiveLambda = editbox(hs.fig,num2str(opts.adaptiveLambda),[editx y editw h],10);
  y = y-h;
  hs.mmfAddSpots = checkbox(hs.fig,'Resolve sub-resolution spots',[x y w h],'',10);
  hs.mmfAddSpots.Value = opts.mmfAddSpots;
  y = y-h;
  t = label(hs.fig,'Max MMF time per frame (min)',[x y labelw h],10);
  hs.maxMmfTime = editbox(hs.fig,num2str(opts.maxMmfTime),[editx y editw h],10);
  y = y-h;
  hs.deconvolve = checkbox(hs.fig,'Deconvolve movies',[x y w h],@deconvolveCB,10);
  y = y-h;
  hs.psfBtn = button(hs.fig,'Select PSF file',[x y labelw/2 h],@psfBtnCB,10);
  if opts.deconvolve
    hs.deconvolve.Value = hs.deconvolve.Max; % On
    hs.psfBtn.Enable = 'on';
  else
    hs.deconvolve.Value = hs.deconvolve.Min; % Off.
    hs.psfBtn.Enable = 'off';
  end

end

function selectDirectoryCB(hObj,event)
  if ~isempty(get(handles.ROIs,'String'))
    r = questdlg('Selecting movie directory will clear existing ROIs. Select?','Warning','Yes','No','No');
    if strcmp(r,'No')
      return
    end
  end

  dirName = uigetdir([], 'Select directory tree containing movies');
  if ~isempty(dirName)
    set(handles.movieDirectory, 'String', dirName);
    % Find movie files.
    movieFiles = kitFindFiles(dirName, kitSupportedFormats());
    % Strip search directory from filenames.
    for i=1:length(movieFiles)
      movieFiles{i} = strrep(movieFiles{i},[dirName filesep],'');
    end
    set(handles.movies,'String',movieFiles,'Value',1:length(movieFiles));
    set(handles.ROIs,'String',[]);
    ROI = [];
  end
end

function addROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  for i=1:length(v)
    [crop,cropSize] = kitCropMovie(fullfile(movieDir,movieFiles{v(i)}));
    for j=1:size(crop,1)
      r = length(jobset.ROI) + 1;
      jobset.ROI(r).movieIdx = v(i);
      jobset.ROI(r).crop = crop(j,:);
      jobset.ROI(r).cropSize = cropSize(j,:);
    end
  end
  populateROIBox();
end

function deleteROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    r = questdlg('Delete selected ROI?','Warning','Yes','No','No');
    if strcmp(r,'Yes')
      jobset.ROI(v) = [];
    end
  end
  populateROIBox();
end

function viewROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    movieFiles = handles.movies.String;
    movieDir = handles.movieDirectory.String;
    kitMovieProj(fullfile(movieDir,movieFiles{jobset.ROI(v).movieIdx}),[],jobset.ROI(v).crop);
  end
end

function spotModeCB(hObj,event)
  if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotDetectValues),'Adaptive'),handles.spotMode))
    handles.adaptiveLambda.Enable = 'on';
  else
    handles.adaptiveLambda.Enable = 'off';
  end
end

function refineModeCB(hObj,event)
  if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotRefineValues),'MMF'),handles.refineMode))
    handles.mmfAddSpots.Enable = 'on';
    handles.maxMmfTime.Enable = 'on';
  else
    handles.mmfAddSpots.Enable = 'off';
    handles.maxMmfTime.Enable = 'off';
  end
end

function autoRadiiCB(hObj,event)
  if handles.autoRadii.Value
    handles.autoRadiidt.Enable = 'on';
    handles.autoRadiiAvgDisp.Enable = 'on';
    handles.minSearchRadius.Enable = 'off';
    handles.maxSearchRadius.Enable = 'off';
    r = computeSearchRadii(str2double(handles.autoRadiidt.String),str2double(handles.autoRadiiAvgDisp.String));
    handles.minSearchRadius.String = num2str(r(1));
    handles.maxSearchRadius.String = num2str(r(2));
  else
    handles.autoRadiidt.Enable = 'off';
    handles.autoRadiiAvgDisp.Enable = 'off';
    handles.minSearchRadius.Enable = 'on';
    handles.maxSearchRadius.Enable = 'on';
  end
end

function useSisterAlignmentCB(hObj,event)
  if handles.useSisterAlignment.Value
    handles.maxSisterAlignmentAngle.Enable = 'on';
  else
    handles.maxSisterAlignmentAngle.Enable = 'off';
  end
end


function exectuteCB(hObj,event)
  updateJobset();
  kitSaveJobset(jobset);
  % Ask which tasks to run.
  taskStrs = {'Spot finding','Plane fitting','Tracking','Sister grouping','Intensity measurement'};
  [tasks,ok] = listdlg('ListString',taskStrs,'InitialValue',1:length(taskStrs),'Prompt','Select tasks to execute','Name','Select tasks...');
  if ok
    % Map to task numbers and add defaults.
    taskMap = [1 2 3 4 8];
    tasks = [taskMap(tasks) 5:7];

    progh = waitbar(0,sprintf('Tracking progress (%d/%d)',0,length(jobset.ROI)));
    kitRunJobs(jobset,'callback',@trackProgress,'tasks',tasks);
    delete(progh);
    uiresume(gcf);
  end

  function trackProgress(idx)
    waitbar(idx/length(jobset.ROIS),progh,sprintf('Tracking progress (%d/%d)',idx,length(jobset.ROI)));
  end
end

function saveCB(hObj,event)
  updateJobset();
  kitSaveJobset(jobset);
  uiresume(gcf);
end

function deconvolveCB(hObj,event)
  if handles.deconvolve.Value
    handles.psfBtn.Enable = 'on';
  else
    handles.psfBtn.Enable = 'off';
  end
end

function psfBtnCB(hObj,event)
  file = hObj.Value
  if exist(file,'file')
    [~,~,ext] = fileparts(file);
    if strcmp(ext,'mat')
      data = load(file);
      % Assume PSF is only variable.
      f = fieldnames(data);
      if length(f) > 1
        msgbox(sprintf('Multiple variables in MAT-file. Using .%s for PSF.',f{1}),'Warning','Warning');
      end
      jobset.psf = data.(f{1});
    else
      msgbox(sprintf('Not MAT-file: %s',file),'Error','Error');
    end
  else
    msgbox(sprintf('File not found: %s',file),'Error','Error');
  end
end

function populateROIBox()
  handles.ROIs.String=[];
  maxMovLen = 32;
  movieFiles = handles.movies.String;
  if ~isempty(movieFiles)
    for i=1:length(jobset.ROI)
      handles.ROIs.String{i} = [strshorten(movieFiles{jobset.ROI(i).movieIdx},maxMovLen) ' [' ...
                          num2str(round(jobset.ROI(i).crop),'%d ') ']'];
    end
  end
  if handles.ROIs.Value < length(jobset.ROI)
    handles.ROIs.Value = 1;
  end
end

function ret=mapStrings(inp,vals)
  if ischar(inp)
    % Map string to index.
    ret = find(cellfun(@(x) ~isempty(x),strfind(vals,inp)));
    assert(~isempty(ret));
  else
    % Map index to index.
    assert(inp > 0 && inp <= length(vals));
    ret = vals{inp};
  end
end

function updateJobset()
  jobset.movieDirectory = handles.movieDirectory.String;
  jobset.movieFiles = handles.movies.String;
  jobset.filename = fullfile(jobset.movieDirectory,[handles.filename.String '.mat']);

  opts = jobset.options;
  opts.coordSystem = mapStrings(handles.coordSys.Value,coordSystemValuesJS);
  opts.coordSystemChannel = str2double(handles.coordSysCh.String);
  for i=1:3
    opts.spotMode{i} = spotDetectValuesJS{handles.spotMode{i}.Value};
    opts.coordMode{i} = spotRefineValuesJS{handles.refineMode{i}.Value};
  end
  if handles.autoRadii.Value
    opts.autoRadiidt = str2double(handles.autoRadiidt.String);
    opts.autoRadiiAvgDisp = str2double(handles.autoRadiiAvgDisp.String);
    r = computeSearchRadii(opts.autoRadiidt,opts.autoRadiiAvgDisp);
  else
    opts.autoRadiidt = [];
    r = zeros(2,1);
    r(1) = str2double(handles.minSearchRadius.String);
    r(2) = str2double(handles.maxSearchRadius.String);
  end
  r = computeUnalignedLaggingRadii(r);
  opts.minSearchRadius = r(1,:);
  opts.maxSearchRadius = r(2,:); % in um
  opts.useSisterAlignment = handles.useSisterAlignment.Value;
  opts.maxSisterAlignmentAngle = str2double(handles.maxSisterAlignmentAngle.String);
  opts.maxSisterSeparation = str2double(handles.maxSisterDist.String);
  opts.minSisterTrackOverlap = str2double(handles.minSisterTrackOverlap.String);
  opts.minSpotsPerFrame = str2double(handles.minSpotsPerFrame.String);
  opts.adaptiveLambda = str2double(handles.adaptiveLambda.String);
  opts.mmfAddSpots = handles.mmfAddSpots.Value;
  opts.maxMmfTime = str2double(handles.maxMmfTime.String);
  opts.deconvole = handles.deconvolve.Value;
  jobset.options = opts;
end

function r=computeSearchRadii(dt,avgDisp)
  r = zeros(2,1);
  r(2) = 6*avgDisp*dt;
  r(1) = 0.1*avgDisp*dt;
end

function r=computeUnalignedLaggingRadii(r)
% Assume unaligned move 3x faster, lagging same speed.
  r = [r 3*r r];
end


function h=label(p,s,pos,sz)
% Create text label
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','text','String',s,'Units','characters','Position',pos,'FontSize',sz,'HorizontalAlignment','left');
end

function h=checkbox(p,s,pos,cb,sz)
% Create checkbox
  if nargin<4
    cb = '';
  end
  if nargin<5
    sz=12;
  end
  h = uicontrol(p,'Style','checkbox','String',s,'Units','characters','Position',pos,'FontSize',sz,'Callback',cb);
end

function h=editbox(p,s,pos,sz)
% Create edit box
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','edit','String',s,'Units','characters','Position',pos,'FontSize',sz,'HorizontalAlignment','left');
end

function h=popup(p,s,pos,cb,sz)
% Create popup menu
  if nargin<4
    cb = '';
  end
  if nargin<5
    sz=12;
  end
  h = uicontrol(p,'Style','popupmenu','String',s,'Units','characters','Position',pos,'FontSize',sz,'Callback',cb);
end

function h=button(p,s,pos,cb,sz)
% Create button
  if nargin<4
    cb = '';
  end
  if nargin<5
    sz=12;
  end
  h = uicontrol(p,'String',s,'Units','characters','Position',pos,'FontSize',sz,'Callback',cb);
end

end % kitGUI
