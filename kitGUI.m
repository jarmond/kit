function jobset=kitGUI(jobset)
% KITGUI Display GUI to setup and run tracking.
%
% Created by: J. W. Armond
% Modified by: C. A. Smith
% Copyright (c) 2017 C. A. Smith

% Download BioFormats, if required.
kitDownloadBioFormats();

if nargin<1 || isempty(jobset)
  jobset = kitDefaultOptions();
end

% Upgrade jobset, if required.
if ~isfield(jobset,'jobsetVersion') || ...
  jobset.jobsetVersion < kitVersion(2)
  jobset = kitJobset(jobset);
end

if ~isfield(jobset,'ROI')
  jobset.ROI = [];
end

jobProcessValues = {'2D/3D only','Chromatic shift'}; % in GUI
jobProcessValuesJS = {'zonly','chrshift'}; % in jobset
coordSystemValues = {'Metaphase plate','Centre of mass'};
coordSystemValuesJS = {'plate','com'};
spotDetectValues = {'Histogram','Neighbour','None'};
spotDetectValuesJS = {'histcut','neighbour','none'};
spotRefineValues = {'MMF','None'};
spotRefineValuesJS = {'gaussian','none'};
maskValues = {'Circle','Semi-circle','Cone'};
maskValuesJS = {'circle','semicirc','cone'};

% Setup GUI.
handles = createControls();
updateControls(jobset);
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);

%% NESTED FUNCTIONS
function hs = createControls()
  colwidth = [55 42 35 35];
  figw = sum(colwidth)+19;
  figh = 47;
  figpos = [100 35 figw figh];
  hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',figpos);
  set(0,'DefaultFigureWindowStyle','normal');
  hs.fig.DockControls = 'off';
  hs.fig.MenuBar = 'none';
  hs.fig.Name = ['KiD ' kitVersion(1)];
  hs.fig.NumberTitle = 'off';
  hs.fig.IntegerHandle = 'off';
  hs.fig.ToolBar = 'none';

  % define some font sizes
  headfont = 16;
  medfont = 14;
  smallfont = 12;

  w=25; h=8;
  hs.logo = uicontrol(hs.fig,'Units','characters','Position',[figw-w-6 1 w h]);
  pos = getpixelposition(hs.logo);
  set(hs.logo,'cdata',imresize(imread('private/kidlogo.png'),pos([4 3])));

  %% ROI selection
  x = 2.5;
  w = colwidth(1);
  h = 1.5;
  lh = 1.5*h;
  toplabely = figh-2;
  % Movies
  y = toplabely-1;
  hs.openBtn = button(hs.fig,'Open existing jobset',[x y 25 lh],@openExistingCB);
  y = y-lh;
  t = label(hs.fig,'1. Movie and ROI selection',[x y w 1.5],14);
  t.FontWeight = 'bold';
  hs.movies = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[x 0.52*figh w 0.285*figh],'Max',inf,'Min',0,'FontSize',smallfont);
  hs.movies.String = jobset.movieFiles;
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.movies);
  jScrollPane.setHorizontalScrollBarPolicy(32);

  y = y-lh;
  hs.movieDirectory = editbox(hs.fig,'',[20 y w-(20-x) h]);
  hs.movieDirectory.Enable = 'inactive';
  hs.movieDirectory.String = jobset.movieDirectory;
  hs.selectDirectory = button(hs.fig,'Select directory',[x y 17.5 lh],@selectDirectoryCB);
  y = y-h;
  hs.labelAvail = label(hs.fig,'Available movies:',[x y 20 h]);

  % ROIs
  label(hs.fig,'ROIs:',[x 0.48*figh 20 1.5]);
  hs.cropROI = button(hs.fig,'Add crop',[x 1+0.415*figh 13.5 2],@addROICB);
  hs.fullROI = button(hs.fig,'Add full',[x+13.9 1+0.415*figh 13.5 2],@skipROICB);
  hs.deleteROI = button(hs.fig,'Delete',[x+13.9*2 1+0.415*figh 13.5 2],@deleteROICB);
  hs.viewROI = button(hs.fig,'View',[x+13.9*3 1+0.415*figh 13.5 2],@viewROICB);
  hs.ROIs = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[x 1 w 0.415*figh],'Callback',@roisCB);
  hs.ROIs.Min = 1;
  % undocumentedmatlab.com hack to add horizontal scrollbars
  jScrollPane = findjobj(hs.ROIs);
  jScrollPane.setHorizontalScrollBarPolicy(32);

  %% Process setup
  x = colwidth(1) + 5;
  w = colwidth(2);
  h = 2;
  t = label(hs.fig,'2. Process setup',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  y = toplabely-2;
  % Job process
  label(hs.fig,'Process type',[x y 14 1.5]);
  hs.jobProc = popup(hs.fig,jobProcessValues,[x+21 y 22 1.5],@jobProcessCB);
  y = y-h;
  % Coordinate system
  hs.coordSysText = label(hs.fig,'Coordinate system',[x y 20 1.5]);
  hs.coordSys = popup(hs.fig,coordSystemValues,[x+21 y 22 1.5],@neighbourOptionsCB);
  y = y-h;
  label(hs.fig,'Coordinate system channel',[x y 30 1.5]);
  for i=1:3
    hs.coordSysCh{i} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String',num2str(i),'Position',[x+21+5.3*i y 6 1.5],'Callback',@coordSysChCB);
  end
  hs.coordSysCh{1}.Value = 1;
  hs.coordSysChNum = 1;
  % Channel modes
  h = 6;
  b = y-h-0.5;
  for i=1:3
    
    p = uipanel(hs.fig,'Units','characters','Position',[x b-(i-1)*h w h],'FontSize',12,'Title',['Channel ' num2str(i)]);
    hs.spotMode{i} = popup(p,spotDetectValues,[22 3 0.45*w 1.5],@spotModeCB);
    hs.spotModeText{i} = label(p,'Spot detection',[1 3 20 1.5]);

    hs.refineMode{i} = popup(p,spotRefineValues,[22 1 0.45*w 1.5],@refineModeCB);
    hs.refineModeText{i} = label(p,'Spot refinement',[1 1 20 1.5]);
    
  end
  %% Options - General options, column 1
  x = sum(colwidth(1:2)) + 8;
  w = 35;
  t = label(hs.fig,'3. Options',[x toplabely 20 1.5],14);
  t.FontWeight = 'bold';
  h = 1.5;
  lh = 1.5*h; % large height
  y = toplabely-lh;
  t = label(hs.fig,'Spot detection & refinement',[x y 30 1.5],12);
  t.FontWeight = 'bold';
  y = y-h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  label(hs.fig,'Min spots per frame',[x y labelw h],10);
  hs.minSpotsPerFrame = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  label(hs.fig,'Max spots per frame',[x y labelw h],10);
  hs.maxSpotsPerFrame = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  y = y-h;
  hs.maxMmfTimeText = label(hs.fig,'Max MMF time per frame (min)',[x y labelw h],10);
  hs.maxMmfTime = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.alphaAText = label(hs.fig,'Weight for intensity restriction:',[x y labelw h],10);
  hs.alphaAchText{1} = label(hs.fig,'Ch.1',[x+labelw-3 y-0.15 4 h],10);
  hs.alphaA{1} = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.alphaAchText{2} = label(hs.fig,'Ch.2',[x+labelw-3 y-0.15 4 h],10);
  hs.alphaA{2} = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.alphaAchText{3} = label(hs.fig,'Ch.3',[x+labelw-3 y-0.15 4 h],10);
  hs.alphaA{3} = editbox(hs.fig,[],[editx y editw h],10);
  y = y-lh;
  
  %% Execution
  h = 2;
  y = 1+4.5*h;
  labelw = 0.5*w;
  t = label(hs.fig,'4. Execution',[x y labelw 1.5],14);
  t.FontWeight = 'bold';
  lh = 1.5*h;
  y = y-lh;
  label(hs.fig,'Jobset name',[x y labelw h],12);
  hs.filename = editbox(hs.fig,'jobset',[x+w-(w-labelw) y (w-labelw) h]);
  btnw = 0.5*w;
  bx = x + w - btnw;
  y = y-h;
  hs.validateMetadata = button(hs.fig,'Validate metadata',[bx y btnw h],@validateCB);
  y = y-h;
  hs.save = button(hs.fig,'Save',[bx y btnw h],@saveCB);
  y = y-h;
  hs.execute = button(hs.fig,'Execute',[bx y btnw h],@executeCB);
  
  %% Options - Chromatic shift options, column 2
  x = sum(colwidth(1:3)) + 11;
  h = 1.5;
  lh = 1.5*h;
  y = toplabely-lh;
  hs.chrShiftOptionsText = label(hs.fig,'Chromatic shift',[x y w 1.5],12);
  hs.chrShiftOptionsText.FontWeight = 'bold';
  y = y-h;
  labelw = 0.75*w;
  editw = 0.2*w;
  editx = x+w-editw;
  hs.chromaticShift = checkbox(hs.fig,'Provide chromatic shift correction',[x y w h],@chromaticShiftCB,10);
  hs.chrShiftPanel = uipanel(hs.fig,'Units','characters','Position',[x y-5.75*h w 5.75*h],'FontSize',10,'Title','Chromatic shift jobsets');
  p = hs.chrShiftPanel;
  y = y-lh;
  editxl = w-4.25*(editw-1);
  editxc = w-1.5*(editw-1);
  editxr = w-editw/3-2;
  edity = 4.25*h;
  hs.chrShiftChanVectText = label(p,'Channel vector',[1 edity-h labelw/3 h*1.5],10);
  hs.chrShiftChanVectText.FontWeight = 'bold';
  hs.chrShiftChanVectText.HorizontalAlignment = 'left';
  hs.chrShiftJobsetText = label(p,'Jobset',[editxl edity-h/2 labelw h],10);
  hs.chrShiftJobsetText.FontWeight = 'bold';
  hs.chrShiftChanOrderText = label(p,'Channel order',[editxc edity-h labelw/3 h*1.5],10);
  hs.chrShiftChanOrderText.FontWeight = 'bold';
  hs.chrShiftChanOrderText.HorizontalAlignment = 'left';
  edity = edity-2*h;
  hs.ch1to2Text = label(p,'1  ->  2',[1 edity labelw h],10);
  hs.ch1to2Arrow = label(p,'->',[(editxc+editxr+0.75)/2 edity editw h],10);
  hs.ch1to2 = button(p,'-',[editxl edity 2*editw h],@ch1to2CB,10);
  hs.ch1to2_ch1num = editbox(p,[],[editxc edity editw/3 h],10);
  hs.ch1to2_ch2num = editbox(p,[],[editxr edity editw/3 h],10);
  edity = edity-h;
  hs.ch1to3Text = label(p,'1  ->  3',[1 edity labelw h],10);
  hs.ch1to3Arrow = label(p,'->',[(editxc+editxr+0.75)/2 edity editw h],10);
  hs.ch1to3 = button(p,'-',[editxl edity 2*editw h],@ch1to3CB,10);
  hs.ch1to3_ch1num = editbox(p,[],[editxc edity editw/3 h],10);
  hs.ch1to3_ch3num = editbox(p,[],[editxr edity editw/3 h],10);
  edity = edity-h;
  hs.ch2to3Text = label(p,'2  ->  3',[1 edity labelw h],10);
  hs.ch2to3Arrow = label(p,'->',[(editxc+editxr+0.75)/2 edity editw h],10);
  hs.ch2to3 = button(p,'-',[editxl edity 2*editw h],@ch2to3CB,10);
  hs.ch2to3_ch2num = editbox(p,[],[editxc edity editw/3 h],10);
  hs.ch2to3_ch3num = editbox(p,[],[editxr edity editw/3 h],10);
  y = y-5.5*h;
  hs.chrShiftFilter = checkbox(hs.fig,'Filter chromatic shift spots',[x y w h],@chrShiftFilterCB,10);
  y = y-h;
  hs.chrShiftamplitudeText = label(hs.fig,'Min spot intensity (% of max)',[x y labelw h],10);
  hs.chrShiftamplitude = editbox(hs.fig,[],[editx y editw h],10);
  y = y-h;
  hs.chrShiftnnDistText = label(hs.fig,'Min spot separation (um)',[x y labelw h],10);
  hs.chrShiftnnDist = editbox(hs.fig,[],[editx y editw h],10);
  y = y-lh;
  
  %% Options - Neighbour spot-detection options, column 2
  t = label(hs.fig,'Neighbour spot detection',[x y w 1.5],12);
  t.FontWeight = 'bold';
  y = y-h;
  hs.neighbourMaskShapeText = label(hs.fig,'Mask shape',[x y labelw h],10);
  hs.neighbourMaskShape = popup(hs.fig,maskValues,[editx-10 y editw+11 h],@neighbourOptionsCB,10);
  y = y-h;
  hs.neighbourMaskRadiusText = label(hs.fig,'Mask radius (um)',[x y labelw h],10);
  hs.neighbourMaskRadius = editbox(hs.fig,[],[editx y editw h],10);
  hs.neighbourOrientPanel = uipanel(hs.fig,'Units','characters','Position',[x y-3*h w 3*h],'FontSize',10,'Title','Neighbour orientation');
  p = hs.neighbourOrientPanel;
  editxl = 10;
  editxr = w-10-editw/3;
  editxc = (editxl+editxr)/2;
  hs.neighbourChanNumText = label(p,'Channel number',[(w-labelw)/2 edity+h labelw h],10);
  hs.neighbourChanNumText.FontWeight = 'bold';
  hs.neighbourChanNumText.HorizontalAlignment = 'center';
  hs.neighbourInnerText = label(p,'inner kchore',[1 edity-h*2/3 labelw/3 2*h],10);
  hs.neighbourOuterText = label(p,'outer kchore',[w-labelw*2/5 edity-h*2/3 labelw/3 2*h],10);
  hs.neighbourOuterText.HorizontalAlignment = 'right';
  hs.neighbourOrient{1} = editbox(p,[],[editxl edity editw/3 h],10);
  hs.neighbourOrient{2} = editbox(p,[],[editxc edity editw/3 h],10);
  hs.neighbourOrient{3} = editbox(p,[],[editxr edity editw/3 h],10);
  y = y-3*h-lh;
  
  %% Options - Intensity options, column 2
  t = label(hs.fig,'Intensity measurement',[x y w 1.5],12);
  t.FontWeight = 'bold';
  y = y-h;
  label(hs.fig,'Measure in channels...',[x y labelw h],10);
  for i=1:3
    hs.intensityExecute{i} = checkbox(hs.fig,num2str(i),[x+w-(4-i)*5 y w h],@intensityOptionsCB,10);
  end
  y = y-h;
  hs.intensityMaskShapeText = label(hs.fig,'Mask shape',[x y labelw h],10);
  hs.intensityMaskShape = popup(hs.fig,maskValues,[editx-10 y editw+11 h],[],10);
  y = y-h;
  hs.intensityMaskRadiusText = label(hs.fig,'Mask radius (um)',[x y labelw h],10);
  hs.intensityMaskRadius = editbox(hs.fig,[],[editx y editw h],10);

  movegui(hs.fig,'center');

end

% Update control status based on contents of jobset.
function updateControls(jobset)
  hs = handles;
  opts = jobset.options;

  if isfield(jobset,'movieDirectory')
    handles.movieDirectory.String = jobset.movieDirectory;
  end
  if isfield(jobset,'filename')
    [~,file] = fileparts(jobset.filename);
    hs.filename.String = file;
  end
  hs.jobProc.Value = mapStrings(opts.jobProcess,jobProcessValuesJS);
  hs.coordSys.Value = mapStrings(opts.coordSystem,coordSystemValuesJS);
  hs.coordSysChNum = opts.coordSystemChannel;
  for i=1:3
    hs.coordSysCh{i}.Value = (i==opts.coordSystemChannel);
    hs.spotMode{i}.Value = mapStrings(opts.spotMode{i},spotDetectValuesJS);
    hs.refineMode{i}.Value = mapStrings(opts.coordMode{i},spotRefineValuesJS);
    if ~isempty(opts.neighbourSpots.zSlices{i})
      hs.zStackMin{i}.String = num2str(opts.neighbourSpots.zSlices{i}(1));
      hs.zStackMax{i}.String = num2str(opts.neighbourSpots.zSlices{i}(end));
    end
    if ~isempty(opts.neighbourSpots.timePoints{i})
      hs.tPointsMin{i}.String = num2str(opts.neighbourSpots.timePoints{i}(1));
      hs.tPointsMax{i}.String = num2str(opts.neighbourSpots.timePoints{i}(end));
    end 
  end
  hs.maxMmfTime.String = num2str(opts.mmf.maxMmfTime);
  for iChan=1:3
    hs.alphaA{iChan}.String = num2str(opts.mmf.alphaA(iChan));
  end
  if isfield(jobset,'psfFile')
    hs.psfFile.String = jobset.psfFile;
  end
  hs.chromaticShift.Value = any(~cellfun('isempty',opts.chrShift.jobset(:)));
  hs.minChrShiftSpots.String = num2str(opts.chrShift.minSpots);
  if ~isempty(opts.chrShift.jobset{1,2})
    hs.ch1to2.String = opts.chrShift.jobset{1,2};
    if ~exist(hs.ch1to2.String,'file')
      kitLog('Chromatic shift jobset for channel 1 to 2 not found: %s',hs.ch1to2.String);
      kitLog('Please relocate this jobset:');
      ch1to2CB();
    end
    hs.ch1to2_ch1num.Enable = 'on';
    hs.ch1to2_ch2num.Enable = 'on';
    hs.ch1to2_ch1num.String = num2str(opts.chrShift.chanOrder{1,2}(1));
    hs.ch1to2_ch2num.String = num2str(opts.chrShift.chanOrder{1,2}(2));
  else
    hs.ch1to2.String = '-';
    hs.ch1to2_ch1num.String = '';
    hs.ch1to2_ch2num.String = '';
    hs.ch1to2_ch1num.Enable = 'off';
    hs.ch1to2_ch2num.Enable = 'off';
  end
  if ~isempty(opts.chrShift.jobset{1,3})
    hs.ch1to2.String = opts.chrShift.jobset{1,3};
    if ~exist(hs.ch1to3.String,'file')
      kitLog('Chromatic shift jobset for channel 1 to 3 not found: %s',hs.ch1to2.String);
      kitLog('Please relocate this jobset:');
      ch1to3CB();
    end
    hs.ch1to3_ch1num.Enable = 'on';
    hs.ch1to3_ch3num.Enable = 'on';
    hs.ch1to3_ch1num.String = num2str(opts.chrShift.chanOrder{1,3}(1));
    hs.ch1to3_ch3num.String = num2str(opts.chrShift.chanOrder{1,3}(2));
  else
    hs.ch1to3.String = '-';
    hs.ch1to3_ch1num.String = '';
    hs.ch1to3_ch3num.String = '';
    hs.ch1to3_ch1num.Enable = 'off';
    hs.ch1to3_ch3num.Enable = 'off';
  end
  if ~isempty(opts.chrShift.jobset{2,3})
    hs.ch2to3.String = opts.chrShift.jobset{2,3};
    if ~exist(hs.ch2to3.String,'file')
      kitLog('Chromatic shift jobset for channel 2 to 3 not found: %s',hs.ch1to2.String);
      kitLog('Please relocate this jobset:');
      ch2to3CB();
    end
    hs.ch2to3_ch2num.Enable = 'on';
    hs.ch2to3_ch3num.Enable = 'on';
    hs.ch2to3_ch2num.String = num2str(opts.chrShift.chanOrder{2,3}(1));
    hs.ch2to3_ch3num.String = num2str(opts.chrShift.chanOrder{2,3}(2));
  else
    hs.ch2to3.String = '-';
    hs.ch2to3_ch2num.String = '';
    hs.ch2to3_ch3num.String = '';
    hs.ch2to3_ch2num.Enable = 'off';
    hs.ch2to3_ch3num.Enable = 'off';
  end
  hs.chrShiftFilter.Value = opts.chrShift.filtering;
  hs.chrShiftamplitude.String = num2str(opts.chrShift.intensityFilter);
  hs.chrShiftnnDist.String = num2str(opts.chrShift.neighbourFilter);
  
  hs.neighbourMaskShape.Value = mapStrings(opts.neighbourSpots.maskShape,maskValuesJS);
  hs.neighbourMaskRadius.String = num2str(opts.neighbourSpots.maskRadius);
  for iChan=1:3
    hs.neighbourOrient{iChan}.String = num2str(opts.neighbourSpots.channelOrientation(iChan));
  end
  
  for iChan=1:3
    hs.intensityExecute{iChan}.Value = opts.intensity.execute(iChan);
  end
  hs.intensityMaskShape.Value = mapStrings(opts.intensity.maskShape,maskValuesJS);
  hs.intensityMaskRadius.String = num2str(opts.intensity.maskRadius);
  hs.minSpotsPerFrame.String = num2str(opts.minSpotsPerFrame);
  hs.maxSpotsPerFrame.String = num2str(opts.maxSpotsPerFrame);

  populateMovieBox();
  populateROIBox();
  jobProcessCB();
  coordSysChCB();
  spotModeCB();
  refineModeCB();
  chromaticShiftCB();
  chrShiftFilterCB();
  neighbourOptionsCB();
  intensityOptionsCB();
  
  handles = hs;
  
end

% Check controls for consistent input.
function tf=checkControls()
  hs = handles;
  v = hs.coordSysChNum;
  if ~isfinite(v) || v<1 || v>3
    errorbox('Invalid channel number for coordinate system channel. Should be a number between 1 and 3');
    tf = false;
    return
  end

  v = str2double(hs.minSpotsPerFrame.String);
  if ~isfinite(v) || v < 0
    errorbox('Invalid value for min spots per frame. Should be a positive number.')
    tf = false;
    return
  end
  v(2) = str2double(hs.maxSpotsPerFrame.String);
  if ~isfinite(v(2)) || v(2) < 0
    errorbox('Invalid value for maximum spots per frame. Should be a positive number.')
    tf = false;
    return
  elseif diff(v)<=0
    errorbox('Invalid values for spots per frame. Maximum number should be larger than the minimum.')
    tf = false;
    return
  end

  if isempty(hs.filename.String)
    errorbox('Jobset name is a required field.');
    tf = false;
    return
  end

  tf = true;
end

function roisCB(hObj,event)
  if isempty(handles.ROIs.Value)
    handles.ROIs.Value = 1;
  end
end

function openExistingCB(hObj,event)
  if ~isempty(get(handles.ROIs,'String'))
    r = questdlg('Selecting existing jobset will clear existing ROIs. Select?','Warning','Yes','No','No');
    if strcmp(r,'No')
      return
    end
  end

  [filename,pathname] = uigetfile('*.mat','Select existing jobset');
  if ~isempty(filename)
    filename = fullfile(pathname,filename);
    try
      jobset = kitLoadJobset(filename);
      if isempty(jobset) || ~isfield(jobset,'kit') || ~jobset.kit
        error('Jobset file corrupt');
      end
      % Upgrade jobset, if required.
      if ~isfield(jobset,'jobsetVersion') || ...
          jobset.jobsetVersion < kitVersion(2)
        jobset = kitJobset(jobset);
      end

      if ~isfield(jobset,'ROI')
        jobset.ROI = [];
      end

      updateControls(jobset);
      h=msgbox('Successfully loaded jobset','Success','Help','modal');
      uiwait(h);
    catch me
      errorbox(sprintf('Error loading jobset %s: %s',filename,me.message));
    end
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
    populateMovieBox();
    set(handles.ROIs,'String',[]);
  end
end

function addROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  if isempty(v)
    errorbox('Must select movies first to add ROIs');
    return
  end
  for i=1:length(v)
    [crop,cropSize] = kitCropMovie(fullfile(movieDir,movieFiles{v(i)}));
    if ~isempty(crop)
      for j=1:size(crop,1)
        r = length(jobset.ROI) + 1;
        jobset.ROI(r).movie = handles.movies.String{v(i)};
        jobset.ROI(r).crop = crop(j,:);
        jobset.ROI(r).cropSize = cropSize(j,:);
      end
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

function skipROICB(hObj,event)
  movieFiles = handles.movies.String;
  movieDir = handles.movieDirectory.String;
  v = handles.movies.Value;
  if isempty(v)
    errorbox('Must select movies first to skip ROIs');
    return
  end
  for i=1:length(v)
    [md,~]=kitOpenMovie(fullfile(movieDir,movieFiles{v(i)}));
    crop = [1 1 md.frameSize(1:2)];
    cropSize = md.frameSize(1:3);
    r = length(jobset.ROI) + 1;
    jobset.ROI(r).movie = handles.movies.String{v(i)};
    jobset.ROI(r).crop = crop;
    jobset.ROI(r).cropSize = cropSize;
  end
  populateROIBox();
end

function viewROICB(hObj,event)
  v = handles.ROIs.Value;
  if ~isempty(v)
    movieDir = handles.movieDirectory.String;
    kitMovieProj(fullfile(movieDir,jobset.ROI(v).movie),[],jobset.ROI(v).crop);
  end
end

function jobProcessCB(hObj,event)
  if strcmp(mapStrings(handles.jobProc.Value,jobProcessValues),'Chromatic shift')
    handles.coordSys.Value = 2;
    handles.coordSysText.Enable = 'off';
    handles.coordSys.Enable = 'off';
    handles.chromaticShift.Value = 0;
    chromaticShiftCB();
    handles.chromaticShift.Enable = 'off';
  else
    handles.coordSysText.Enable = 'on';
    handles.coordSys.Enable = 'on';
    handles.chromaticShift.Enable = 'on';
  end
  neighbourOptionsCB();
end

function coordSysChCB(hObj,event)
  if exist('hObj')
    chan = str2double(hObj.String);
    handles.coordSysChNum = chan;
    handles.spotMode{chan}.Value = 1;
    handles.refineMode{chan}.Value = 1;
    for notChan = setdiff(1:3,chan)
      handles.coordSysCh{notChan}.Value = 0;
      handles.spotMode{notChan}.Value = 3;
    end
  end
  spotModeCB();
  neighbourOptionsCB();
end


function spotModeCB(hObj,event)
  for i=1:3
    if strcmp(mapStrings(handles.spotMode{i}.Value,spotDetectValues),'None')
      handles.refineModeText{i}.Enable = 'off';
      handles.refineMode{i}.Value = 2;
      handles.refineMode{i}.Enable = 'off';
      handles.zSliceText{i}.Enable = 'off';
      handles.zStackMinText{i}.Enable = 'off';
      handles.zStackMin{i}.Enable = 'off';
      handles.zStackMaxText{i}.Enable = 'off';
      handles.zStackMax{i}.Enable = 'off';
      handles.tPointsText{i}.Enable = 'off';
      handles.tPointsMinText{i}.Enable = 'off';
      handles.tPointsMin{i}.Enable = 'off';
      handles.tPointsMaxText{i}.Enable = 'off';
      handles.tPointsMax{i}.Enable = 'off';
    else
      handles.refineModeText{i}.Enable = 'on';  
      handles.refineMode{i}.Enable = 'on';
      handles.zSliceText{i}.Enable = 'on';
      handles.zStackMinText{i}.Enable = 'on';
      handles.zStackMin{i}.Enable = 'on';
      handles.zStackMaxText{i}.Enable = 'on';
      handles.zStackMax{i}.Enable = 'on';
      handles.tPointsText{i}.Enable = 'on';
      handles.tPointsMinText{i}.Enable = 'on';
      handles.tPointsMin{i}.Enable = 'on';
      handles.tPointsMaxText{i}.Enable = 'on';
      handles.tPointsMax{i}.Enable = 'on';
    end
  end
  refineModeCB();
  neighbourOptionsCB();
end

function refineModeCB(hObj,event)
  if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotRefineValues),'MMF'),handles.refineMode))
    handles.maxMmfTimeText.Enable = 'on';
    handles.maxMmfTime.Enable = 'on';
    handles.alphaAText.Enable = 'on';
    for iChan=1:3
      if strcmp(mapStrings(handles.refineMode{iChan}.Value,spotRefineValues),'MMF');
        handles.alphaAchText{iChan}.Enable = 'on';
        handles.alphaA{iChan}.Enable = 'on';
      else
        handles.alphaAchText{iChan}.Enable = 'off';
        handles.alphaA{iChan}.Enable = 'off';
      end
    end
  else
    handles.maxMmfTimeText.Enable = 'off';
    handles.maxMmfTime.Enable = 'off';
    handles.alphaAText.Enable = 'off';
    for iChan=1:3
      handles.alphaAchText{iChan}.Enable = 'off';
      handles.alphaA{iChan}.Enable = 'off';
    end
  end
end

function chromaticShiftCB(hObj,event)
  if handles.chromaticShift.Value
    handles.minChrShiftSpotsText.Enable = 'on';
    handles.minChrShiftSpots.Enable = 'on';
    handles.chrShiftPanel.ForegroundColor = [0 0 0];
    handles.chrShiftChanVectText.Enable = 'on';
    handles.chrShiftJobsetText.Enable = 'on';
    handles.chrShiftChanOrderText.Enable = 'on';
    handles.ch1to2Text.Enable = 'on';
    handles.ch1to2Arrow.Enable = 'on';
    handles.ch1to2.Enable = 'on';
    handles.ch1to3.Enable = 'on';
    handles.ch1to3Text.Enable = 'on';
    handles.ch1to3Arrow.Enable = 'on';
    handles.ch2to3Text.Enable = 'on';
    handles.ch2to3Arrow.Enable = 'on';
    handles.ch2to3.Enable = 'on';
    if all(~strcmp(handles.ch1to2.String,{'-','Unknown source'}))
      handles.ch1to2_ch1num.Enable = 'on';
      handles.ch1to2_ch2num.Enable = 'on';
    else
      handles.ch1to2_ch1num.String = '';
      handles.ch1to2_ch2num.String = '';
      handles.ch1to2_ch1num.Enable = 'off';
      handles.ch1to2_ch2num.Enable = 'off';
    end
    if all(~strcmp(handles.ch1to3.String,{'-','Unknown source'}))
      handles.ch1to3_ch1num.Enable = 'on';
      handles.ch1to3_ch3num.Enable = 'on';
    else
      handles.ch1to3_ch1num.String = '';
      handles.ch1to3_ch3num.String = '';
      handles.ch1to3_ch1num.Enable = 'off';
      handles.ch1to3_ch3num.Enable = 'off';
    end
    if all(~strcmp(handles.ch2to3.String,{'-','Unknown source'}))
      handles.ch2to3_ch2num.Enable = 'on';
      handles.ch2to3_ch3num.Enable = 'on';
    else
      handles.ch2to3_ch2num.String = '';
      handles.ch2to3_ch3num.String = '';
      handles.ch2to3_ch2num.Enable = 'off';
      handles.ch2to3_ch3num.Enable = 'off';
    end
    handles.chrShiftFilter.Enable = 'on';
  else
    handles.minChrShiftSpotsText.Enable = 'off';
    handles.minChrShiftSpots.Enable = 'off';
    handles.chrShiftPanel.ForegroundColor = [0.5 0.5 0.5];
    handles.chrShiftChanVectText.Enable = 'off';
    handles.chrShiftJobsetText.Enable = 'off';
    handles.chrShiftChanOrderText.Enable = 'off';
    handles.ch1to2Text.Enable = 'off';
    handles.ch1to2Arrow.Enable = 'off';
    handles.ch1to2.Enable = 'off';
    handles.ch1to3Text.Enable = 'off';
    handles.ch1to3Arrow.Enable = 'off';
    handles.ch1to3.Enable = 'off';
    handles.ch2to3Text.Enable = 'off';
    handles.ch2to3Arrow.Enable = 'off';
    handles.ch2to3.Enable = 'off';
    handles.ch1to2_ch1num.Enable = 'off';
    handles.ch1to2_ch2num.Enable = 'off';
    handles.ch1to3_ch1num.Enable = 'off';
    handles.ch1to3_ch3num.Enable = 'off';
    handles.ch2to3_ch2num.Enable = 'off';
    handles.ch2to3_ch3num.Enable = 'off';
    handles.chrShiftFilter.Value = 0;
    chrShiftFilterCB();
    handles.chrShiftFilter.Enable = 'off';
  end
end

function chrShiftFilterCB(hObj,event)
  if handles.chrShiftFilter.Value
    handles.chrShiftamplitudeText.Enable = 'on';
    handles.chrShiftamplitude.Enable = 'on';
    handles.chrShiftnnDistText.Enable = 'on';
    handles.chrShiftnnDist.Enable = 'on';
  else
    handles.chrShiftamplitudeText.Enable = 'off';
    handles.chrShiftamplitude.Enable = 'off';
    handles.chrShiftnnDistText.Enable = 'off';
    handles.chrShiftnnDist.Enable = 'off';
  end
end

function ch1to2CB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  handles.ch1to2.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{1,2} = file;
  handles.ch1to2_ch1num.Enable = 'on';
  handles.ch1to2_ch1num.String = jobset.options.chrShift.chanOrder{1,2}(1);
  handles.ch1to2_ch2num.Enable = 'on';
  handles.ch1to2_ch2num.String = jobset.options.chrShift.chanOrder{1,2}(2);
end

function ch1to3CB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  handles.ch1to3.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{1,3} = file;
  handles.ch1to3_ch1num.Enable = 'on';
  handles.ch1to3_ch1num.String = jobset.options.chrShift.chanOrder{1,3}(1);
  handles.ch1to3_ch3num.Enable = 'on';
  handles.ch1to3_ch3num.String = jobset.options.chrShift.chanOrder{1,3}(2);
end

function ch2to3CB(hObj,event)
  [file,path] = uigetfile('*.mat','Select chromatic shift jobset file');
  if isequal(file,0)
    return
  end
  handles.ch2to3.String = file;
  file = fullfile(path,file);
  jobset.options.chrShift.jobset{2,3} = file;
  handles.ch2to3_ch2num.Enable = 'on';
  handles.ch2to3_ch2num.String = jobset.options.chrShift.chanOrder{2,3}(1);
  handles.ch2to3_ch3num.Enable = 'on';
  handles.ch2to3_ch3num.String = jobset.options.chrShift.chanOrder{2,3}(2);
end

function neighbourOptionsCB(hObj,event)
  neighChans = cellfun(@(x) strcmp(mapStrings(x.Value,spotDetectValues),'Neighbour'),handles.spotMode);
  switch sum(neighChans)
    case 0
      handles.neighbourMaskShapeText.Enable = 'off';
      handles.neighbourMaskShape.Enable = 'off';
      handles.neighbourMaskRadiusText.Enable = 'off';
      handles.neighbourMaskRadius.Enable = 'off';
      handles.neighbourOrientPanel.ForegroundColor = [0.5 0.5 0.5];
      handles.neighbourChanVectText.Enable = 'off';
      handles.neighbourChanNumText.Enable = 'off';
      handles.neighbourInnerText.Enable = 'off';
      for iChan=1:3
        handles.neighbourOrient{iChan}.Enable = 'off';
      end
      handles.neighbourOuterText.Enable = 'off';
    case 1
      handles.neighbourMaskShapeText.Enable = 'on';
      handles.neighbourMaskShape.Enable = 'on';
      handles.neighbourMaskRadiusText.Enable = 'on';
      handles.neighbourMaskRadius.Enable = 'on';
      handles.neighbourOrientPanel.ForegroundColor = [0 0 0];
      handles.neighbourChanVectText.Enable = 'on';
      handles.neighbourChanNumText.Enable = 'on';
      handles.neighbourInnerText.Enable = 'on';
      handles.neighbourOuterText.Enable = 'on';
      for iChan=1:3
        if ismember(iChan,[handles.coordSysChNum find(neighChans)])
          handles.neighbourOrient{iChan}.Enable = 'on';
        else
          handles.neighbourOrient{iChan}.Enable = 'off';
        end
      end
    case 2
      handles.neighbourMaskShapeText.Enable = 'on';
      handles.neighbourMaskShape.Enable = 'on';
      handles.neighbourMaskRadiusText.Enable = 'on';
      handles.neighbourMaskRadius.Enable = 'on';
      handles.neighbourOrientPanel.ForegroundColor = [0 0 0];
      handles.neighbourChanVectText.Enable = 'on';
      handles.neighbourChanNumText.Enable = 'on';
      handles.neighbourInnerText.Enable = 'on';
      handles.neighbourOuterText.Enable = 'on';
      for iChan = 1:3
        handles.neighbourOrient{iChan}.Enable = 'on';
      end
  end
  if strcmp(mapStrings(handles.jobProc.Value,jobProcessValues),'Chromatic shift') || ...
          strcmp(mapStrings(handles.coordSys.Value,coordSystemValues),'Centre of mass')
    handles.neighbourMaskShapeText.Enable = 'off';
    handles.neighbourMaskShape.Value = 1;
    handles.neighbourMaskShape.Enable = 'off';
  end
end

function intensityOptionsCB(hObj,event)
  
  if any(cellfun(@(x) x.Value,handles.intensityExecute))
      handles.intensityMaskShapeText.Enable = 'on';
      handles.intensityMaskShape.Enable = 'on';
      handles.intensityMaskRadiusText.Enable = 'on';
      handles.intensityMaskRadius.Enable = 'on';
  else
      handles.intensityMaskShapeText.Enable = 'off';
      handles.intensityMaskShape.Enable = 'off';
      handles.intensityMaskRadiusText.Enable = 'off';
      handles.intensityMaskRadius.Enable = 'off';
  end
    
end

function executeCB(hObj,event)
  if ~checkControls()
    return
  end
  updateJobset();
  kitSaveJobset(jobset);
  % % Ask which tasks to run.
  %  taskStrs = {'Spot finding','Plane fitting','Tracking','Sister grouping','Intensity measurement'};
  %  [tasks,ok] = listdlg('ListString',taskStrs,'InitialValue',1:length(taskStrs),'PromptString','Select tasks to execute','Name','Select tasks...');
  %  if ok
    tasks = [1,2,4,5];
    % Map to task numbers and add defaults.
    taskMap = [1 2 3 4 9];
    tasks = [taskMap(tasks) 6:8];
    if any(cellfun(@(x) strcmp(mapStrings(x.Value,spotDetectValues),'Neighbour'),handles.spotMode))
      tasks = [tasks 5];
    end
    
    execmode = 'serial'; %parallel option removed
    progh = waitbar(0,sprintf('Tracking progress (%d/%d)',0,length(jobset.ROI)));
    kitRunJobs(jobset,'callback',@trackProgress,'tasks',tasks,'exec',execmode);
    delete(progh);
    uiresume(gcf);
%  end

  function trackProgress(idx)
    waitbar(idx/length(jobset.ROI),progh,sprintf('Tracking progress (%d/%d)',idx,length(jobset.ROI)));
  end
end

function saveCB(hObj,event)
  if ~checkControls()
    return
  end
  updateJobset();
  kitSaveJobset(jobset);
  uiresume(gcf);
end

function validateCB(hObj,event)
  updateJobset();
  % check that ROIs are already provided
  if ~isfield(jobset,'ROI')
    errorbox('Must provide ROIs prior to validating their metadata.')
  else
    % loop over each movie
    for iMov = 1:length(jobset.ROI)
      [jobset,applyAll] = kitValidateMetadata(jobset,iMov);
      % skip showing remaining movies, just save metadata for all
      if applyAll
        for jMov = iMov+1:length(jobset.ROI)
          jobset.metadata{jMov} = jobset.metadata{iMov};
        end
        handles.validateMetadata.String = 'Re-validate...';
        break
      end
    end
    
  end
  % show that movies have been validated
  handles.validateMetadata.String = 'Re-validate...';
  
end

function populateMovieBox()
  movieDir = handles.movieDirectory.String;
  if isempty(movieDir)
    handles.movies.String = [];
  else
    % Find movie files.
    movieFiles = kitFindFiles(movieDir, kitSupportedFormats(),1,0,1);
    % Strip search directory from filenames.
    for i=1:length(movieFiles)
      movieFiles{i} = strrep(movieFiles{i},[movieDir filesep],'');
    end
    set(handles.movies,'String',movieFiles,'Value',1:length(movieFiles));
  end
end

function populateROIBox()
  handles.ROIs.String=[];
  maxMovLen = 32;
  movieFiles = handles.movies.String;
  handles.ROIs.Value = 0; % Keep MATLAB quiet about invalid selection.
  if ~isempty(movieFiles)
    for i=1:length(jobset.ROI)
      handles.ROIs.String{i} = [strshorten(jobset.ROI(i).movie,maxMovLen) ' [' ...
                          num2str(round(jobset.ROI(i).crop),'%d ') ']'];
    end
  end
  if (handles.ROIs.Value <= 0 && ~isempty(handles.ROIs.String)) || handles.ROIs.Value > length(jobset.ROI)
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
  opts.jobProcess = mapStrings(handles.jobProc.Value,jobProcessValuesJS);
  opts.coordSystem = mapStrings(handles.coordSys.Value,coordSystemValuesJS);
  opts.coordSystemChannel = handles.coordSysChNum;
  for i=1:3
    opts.spotMode{i} = spotDetectValuesJS{handles.spotMode{i}.Value};
    opts.coordMode{i} = spotRefineValuesJS{handles.refineMode{i}.Value};
  end
  opts.minSpotsPerFrame = str2double(handles.minSpotsPerFrame.String);
  opts.maxSpotsPerFrame = str2double(handles.maxSpotsPerFrame.String);
  mmf = opts.mmf;
  mmf.maxMmfTime = str2double(handles.maxMmfTime.String);
  for iChan=1:3
    mmf.alphaA(iChan) = str2double(handles.alphaA{iChan}.String);
  end
  opts.mmf = mmf;

  if handles.chromaticShift.Value
    chrShift = opts.chrShift;
    chrShift.minSpots = str2double(handles.minChrShiftSpots.String);
    chrShift.chanOrder{1,2}(1) = str2double(handles.ch1to2_ch1num.String);
    chrShift.chanOrder{1,2}(2) = str2double(handles.ch1to2_ch2num.String);
    chrShift.chanOrder{1,3}(1) = str2double(handles.ch1to3_ch1num.String);
    chrShift.chanOrder{1,3}(2) = str2double(handles.ch1to3_ch3num.String);
    chrShift.chanOrder{2,3}(1) = str2double(handles.ch2to3_ch2num.String);
    chrShift.chanOrder{2,3}(2) = str2double(handles.ch2to3_ch3num.String);
    if handles.chrShiftFilter.Value
      chrShift.filtering = 1;
      chrShift.intensityFilter = str2double(handles.chrShiftamplitude.String);
      chrShift.neighbourFilter = str2double(handles.chrShiftnnDist.String);
    end
    result = getChromaticShiftResults(chrShift);
    chrShift.result = result;
    opts.chrShift = chrShift;
  end
  neighbourSpots = opts.neighbourSpots;
  neighbourSpots.maskShape = mapStrings(handles.neighbourMaskShape.Value,maskValuesJS);
  neighbourSpots.maskRadius = str2double(handles.neighbourMaskRadius.String);
  for iChan=1:3
    neighbourSpots.channelOrientation(iChan) = str2double(handles.neighbourOrient{iChan}.String);
  end
  opts.neighbourSpots = neighbourSpots;
  intensity = opts.intensity;
  for iChan=1:3
    intensity.execute(iChan) = handles.intensityExecute{iChan}.Value;
  end
  intensity.maskShape = mapStrings(handles.intensityMaskShape.Value,maskValuesJS);
  intensity.maskRadius = str2double(handles.intensityMaskRadius.String);
  opts.intensity = intensity;
  
  jobset.options = opts;
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

function cellResult = getChromaticShiftResults(chrShift)
  cellResult = chrShift.result;
  for i = 1:2
    for j = (i+1):3
      jS = chrShift.jobset{i,j};
	  if i==j || isempty(jS) || strcmp(jS,'Unknown source'); continue; end
      jS = kitLoadJobset(jS);
      neighFilt = chrShift.neighbourFilter;
      intFilt = chrShift.intensityFilter/100;
      mS = kitLoadAllJobs(jS);
      if chrShift.filtering
        mS = chrsFilterSpots(mS,'revert',1, ...
          'neighbourFilter',neighFilt,'intensityFilter',intFilt, ...
          'referenceChan',handles.coordSysChNum);
      end
      [result,~] = chrsCalculateChromaticShift(mS,[i j],...
          'filtered',1);
      cellResult{i,j} = result; cellResult{j,i} = result.*[-1 -1 -1 1 1 1];
    end
  end       
end

end % kitGUI
