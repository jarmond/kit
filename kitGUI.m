function jobset=kitGUI(jobset)
% KITGUI Display GUI to setup and run tracking.
%
% Copyright (c) 2015 Jonathan W. Armond

if nargin<1 || isempty(jobset)
  jobset = kitDefaultOptions();
end

% Setup GUI. TODO feed in default jobset options
handles = createControls();
handles.fig.Visible = 'on';

ROIs = struct('movieIdx',[],'crop',[],'cropSize',[]);


function hs = createControls()
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 138 44]);
  movegui(hs.fig,'center');
  figpos = hs.fig.Position

  w=25; h=8;
  %hs.logo = uicontrol(hs.fig,'Units','characters','Position',[figpos(3)-w figpos(4)-h w h]);
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[figpos(3)-w-3 1 w h]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kitlogo.png'),pos([4 3])));

  %% ROI selection
  w = 53;
  toplabely = figpos(4)-2;
  % Movies
  hs.movies = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[2.5 23 w 15.15],'Max',inf,'Min',0);
  hs.movieDirectory = editbox(hs.fig,'',[2.5 40.1 w 1.7]);
  hs.movieDirectory.Enable = 'inactive';
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
  x = w + 3;
  w = 42;
  uicontrol(hs.fig,'Style','text','String','Tracking setup','Units','characters','Position',[x toplabely 20 1.5],'FontSize',14,'FontWeight','bold','HorizontalAlignment','left');
  % Coordinate system
  label(hs.fig,'Coordinate system',[x 40 20 1.5]);
  hs.coordSys = popup(hs.fig,{'plate','poles','image'},[74 40.1 20 1.5]);
  label(hs.fig,'Coordinate system channel',[x 38 30 1.5]);
  hs.coordSysCh = editbox(hs.fig,'1',[87 38 6 1.5]);

  % Channel modes
  b = 32;
  h = 5.5;
  for i=1:3
    p = uipanel(hs.fig,'Units','characters','Position',[x b-(i-1)*h w h],'FontSize',12,'Title',['Channel ' num2str(i)]);
    hs.spotMode{i} = popup(p,{'Histogram','Wavelet','None'},[22 2.5 0.45*w 1.5]);
    label(p,'Spot detection',[1 2.5 20 1.5]);

    hs.refineMode{i} = popup(p,{'Centroid','MMF','None'},[22 0.5 0.45*w 1.5]);
    label(p,'Spot refinement',[1 0.5 20 1.5]);
  end

  b = b-2*h;
  % Tasks
  tasks = {'Spot finding','Plane fitting','Tracking','Sister grouping','Intensity measurement'};
  h = 1.75;
  panelh = h*(1+length(tasks));
  p = uipanel(hs.fig,'Units','characters','Position',[x b-panelh w panelh],'FontSize',12,'Title','Tasks');
  for i=1:length(tasks)
     hs.tasks{i} = checkbox(p,tasks{i},[1 panelh-h*(i+1) 30 h]);
     hs.tasks{i}.Value = hs.tasks{i}.Max;
  end

  %% Options
  x = x+w+2;
  w = 35;
  t = label(hs.fig,'Tracking options',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  h = 1.5;
  lh = 1.5*h; % large height
  y = toplabely-h;
  hs.autoRadii = checkbox(hs.fig,'Calculate search radii from dt',[x y w h],10);
  hs.autoRadii.Value = hs.autoRadii.Max;
  y = y-h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  t = label(hs.fig,'Frame dt',[x y labelw h],10);
  hs.autoRadiidt = editbox(hs.fig,'2',[editx y editw h],10);
  y = y-h;
  t = label(hs.fig,'Min search radius',[x y labelw h],10);
  hs.minSearchRadius = editbox(hs.fig,'0.1',[editx y editw h],10);
  hs.minSearchRadius.Enable = 'off';
  y = y-h;
  t = label(hs.fig,'Max search radius',[x y labelw h],10);
  hs.maxSearchRadius = editbox(hs.fig,'0.75',[editx y editw h],10);
  hs.maxSearchRadius.Enable = 'off';
  y = y-h;
  % Adjust text box pos for multiple lines.
  t = label(hs.fig,'Max angle between sisters and plate normal (deg)',[x y-h/2 labelw lh],10);
  hs.maxSisterAlignmentAngle = editbox(hs.fig,'30',[editx y editw h],10);
  t.Extent
  y = y-lh;
  t = label(hs.fig,'Max average distance between sisters (um)',[x y-h/2 labelw lh],10);
  hs.maxSisterDist = editbox(hs.fig,'30',[editx y editw h],10);
  y = y-lh;
  t = label(hs.fig,'Min overlap between sister tracks',[x y-h/2 labelw lh],10);
  hs.minSisterTrackOverlap = editbox(hs.fig,'30',[editx y editw h],10);
  y = y-lh;
  t = label(hs.fig,'Min spots per frame',[x y labelw h],10);
  hs.minSpotsPerFrame = editbox(hs.fig,'30',[editx y editw h],10);
  y = y-h;
  hs.waveletPrefilter = checkbox(hs.fig,'Prefilter before detecting spots',[x y w h],10);
  y = y-h;
  hs.mmfAddSpots = checkbox(hs.fig,'Resolve sub-resolution spots',[x y w h],10);
  y = y-h;
  t = label(hs.fig,'Max MMF time per frame (min)',[x y labelw h],10);
  hs.maxMmfTime = editbox(hs.fig,'30',[editx y editw h],10);

  %% Execution
  y = y-2*h;
  t = label(hs.fig,'Execution',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  h = 2;
  y = y-h;
  t = label(hs.fig,'Jobset name',[x y w h])
  t.FontWeight = 'bold';
  y = y-h;
  editbox(hs.fig,'jobset',[x y w h]);
  btnw = 0.5*w;
  x = x + 0.5*(w-btnw);
  y = y-2*h;
  hs.save = button(hs.fig,'Save',[x y btnw h]);
  y = y-h;
  hs.execute = button(hs.fig,'Execute',[x y btnw h]);


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
    ROIs = [];
  end
end

function addROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  for i=1:length(v)
    [crop,cropSize] = kitCropMovie(fullfile(movieDir,movieFiles{v(i)}));
    for j=1:size(crop,1)
      r.movieIdx = v(i);
      r.crop = crop(j,:);
      r.cropSize = cropSize(j,:);
      ROIs = [ROIs; r];
    end
  end
  populateROIBox();
end

function deleteROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    r = questdlg('Delete selected ROI?','Warning','Yes','No','No');
    if strcmp(r,'Yes')
      ROIs(v) = [];
    end
  end
  populateROIBox();
end

function viewROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    movieFiles = handles.movies.String;
    movieDir = handles.movieDirectory.String;
    kitMovieProj(fullfile(movieDir,movieFiles{ROIs(v).movieIdx}),[],ROIs(v).crop);
  end
end

function populateROIBox()
  handles.ROIs.String=[];
  maxMovLen = 30;
  movieFiles = handles.movies.String;
  for i=1:length(ROIs)
    handles.ROIs.String{i} = [strshorten(movieFiles{ROIs(i).movieIdx},maxMovLen) ' [' ...
                        num2str(round(ROIs(i).crop),'%d ') ']'];
  end
  if handles.ROIs.Value < length(ROIs)
    handles.ROIs.Value = 1;
  end
end

function h=label(p,s,pos,sz)
% Create text label
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','text','String',s,'Units','characters','Position',pos,'FontSize',sz,'HorizontalAlignment','left');
end

function h=checkbox(p,s,pos,sz)
% Create checkbox
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','checkbox','String',s,'Units','characters','Position',pos,'FontSize',sz);
end

function h=editbox(p,s,pos,sz)
% Create edit box
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','edit','String',s,'Units','characters','Position',pos,'FontSize',sz,'HorizontalAlignment','left');
end

function h=popup(p,s,pos,sz)
% Create popup menu
  if nargin<4
    sz=12;
  end
  h = uicontrol(p,'Style','popupmenu','String',s,'Units','characters','Position',pos,'FontSize',sz);
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
