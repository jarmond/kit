function kitSpotCategories(jobset,desc)
% KITSPOTCATEGORIES(JOBSET) Displays a GUI to allow handling of categories
% for a given JOBSET.
%
%    KITSPOTCATEGORIES(JOBSET)
%
% Copyright (c) 2019 C. A. Smith
           
% Setup GUI.
if nargin<2 || isempty(desc)
    desc = 0;
end
dockStatus = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal');
handles = createControls(jobset);
handles.jobset = jobset;
handles.fig.Visible = 'on';
uiwait(gcf);
close(gcf);
set(0,'DefaultFigureWindowStyle',dockStatus);

function hs = createControls(jobset)
  
    % Create figure.
    figw = 50;
    figh = 28+8*desc;
    hs.fig = figure('Visible','off','Resize','off','Units','characters','Position',[100 35 figw figh]);
    hs.fig.DockControls = 'off';
    hs.fig.MenuBar = 'none';
    hs.fig.Name = 'Spot categorisation';
    hs.fig.NumberTitle = 'off';
    hs.fig.IntegerHandle = 'off';
    hs.fig.ToolBar = 'none';

    % Define font sizes.
    largefont = 16;
    medfont = 14;
    smallfont = 12;
    tinyfont = 10;

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

    % Get the data.
    job = kitLoadAllJobs(jobset);
    % Get list of categories.
    allCats = fieldnames(job{1}.categories);
    hs.nCats = length(allCats);
    
    % Add a category subtitle.
    x = dx; labw = w;
    t = label(hs.fig,'1. Add a new category',[x y labw h],medfont);
    t.FontWeight = 'bold';
    y = y-h;
    
    % Brief explanation of use.
    if desc
        labw = w; nLines=3;
        y = y-(nLines-1)*h;
        label(hs.fig,'To add a new category, type a new name, then click the ''Add category'' button. Tick the check button to also see previously allocated categories.',[x y labw nLines*h],smallfont);
        y = y-h;
    end
    
    % Add a category.
    x = dx;
    labw = 27;
    label(hs.fig,'Select spots in channel...',[x y labw h]);
    radx = x+(labw+ddx); radw = (figw-(2*dx+labw))/4;
    for i=1:3
        hs.spotSelCh{i} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String',num2str(i),...
            'Position',[radx y radw h],'Callback',@spotSelChCB,'FontSize',tinyfont);
        radx = radx+radw;
    end
    hs.spotSelCh{jobset.options.coordSystemChannel}.Value = 1;
    hs.selChan = 1;
    x = dx; y = y-h;
    labw = figw-(2*dx);
    hs.showExisting = checkbox(hs.fig,'Show existing categories during selection',[x y labw h]);
    btnw = 12; btnh = 2;
    x = figw-(btnw+dx); y = y-h;
    hs.addBtn = button(hs.fig,'Add',[x y-0.25 btnw btnh],@addCB);
    x = dx;
    labw = 17.5;
    label(hs.fig,'New category label:',[x y labw h]);
    x = x+(labw+ddx); labw = w-(btnw+dx+labw+ddx);
    hs.newCategory = editbox(hs.fig,'',[x y labw h]);
    y = y-lh;
    
    % Add a category subtitle.
    x = dx; labw = w;
    t = label(hs.fig,'2. Handle existing categories',[x y labw h],medfont);
    t.FontWeight = 'bold';
    y = y-h;
    
    % Brief explanation of use.
    if desc
        labw = w; nLines=2;
        y = y-(nLines-1)*h;
        label(hs.fig,'Existing categories can be edited, shown or deleted by highlighting them individually in the list below.',[x y labw nLines*h],smallfont);
        y = y-h;
    end
    
    % Categories panel
    labx = dx;
    labw = figw-2*dx;
    hs.categoriesLabel = label(hs.fig,'List of categories:',[labx y labw h]);
    panw = labw; panh = 3*h;
    panx = labx; y = y-panh;
    hs.categories = uicontrol(hs.fig,'Style','listbox','Units','characters','Position',[panx y panw panh],'Max',1,'Min',0,'FontSize',smallfont);
    hs.categories.String = allCats;

    % Run buttons.
    x = dx; y = y-lh;
    btnw = (figw-(2*dx+3*ddx))/4; btnh = 2;
    hs.editBtn = button(hs.fig,'Edit',[x y btnw btnh],@editCB);
    x = x+(btnw+ddx);
	hs.delBtn = button(hs.fig,'Delete',[x y btnw btnh],@delCB);
    x = dx; y = y-h;
    labw = figw-(2*dx);
    hs.showExistingEdit = checkbox(hs.fig,'Show other categories during editing',[x y labw h]);
    btnw = 2*btnw+ddx;
    x = dx; y = y-lh;
    
    % Checks subtitle.
    x = dx; labw = w;
    t = label(hs.fig,'3. Overview',[x y labw h],medfont);
    t.FontWeight = 'bold';
    y = y-lh;
    
    % Show all categories button.
    hs.showBtn = button(hs.fig,'Show all categories',[x y btnw btnh],@showCatsCB);
    x = dx; y = y-h;
    labw = 27;
    label(hs.fig,'Show categories in channel...',[x y labw h]);
    radx = x+(labw+ddx); radw = (figw-(2*dx+labw))/4;
    for i=1:3
        hs.spotShowCh{i} = uicontrol('Parent',hs.fig,'Units','characters','Style','radio','String',num2str(i),...
            'Position',[radx y radw h],'Callback',@spotShowChCB,'FontSize',tinyfont);
        radx = radx+radw;
    end
    hs.spotShowCh{jobset.options.coordSystemChannel}.Value = 1;
    hs.showChan = 1;

    % Close button.
    btnw = 15;
    x = figw-btnw-dx; y = 1;
    hs.closeBtn = button(hs.fig,'Close',[x y btnw btnh],@closeCB);

    movegui(hs.fig,'center');
  
end

%% Callback functions

function spotSelChCB(hObj,event)
  chan = str2double(hObj.String);
  handles.selChan = chan;
  for notChan = setdiff(1:3,chan)
      handles.spotSelCh{notChan}.Value = 0;
  end
end

function spotShowChCB(hObj,event)
  chan = str2double(hObj.String);
  handles.showChan = chan;
  for notChan = setdiff(1:3,chan)
      handles.spotShowCh{notChan}.Value = 0;
  end
end

function addCB(hObj,event)
  hs = handles;
  if hs.nCats == 4
      errorbox('Cannot store more than 4 categories per jobset.');
      return
  end
  newCat = strrep(hs.newCategory.String,' ','_');
  jS = hs.jobset;
  kitNewCategory(jS,'channel',hs.selChan,'category',newCat, ...
      'showExisting',hs.showExisting.Value)
  
  % update handles and store
  hs.nCats = hs.nCats + 1;
  hs.categories.String{end+1} = newCat;
  hs.newCategory.String = '';
  handles = hs;
  
end

function editCB(hObj,event)
  hs = handles;
  iCat = hs.categories.Value;
  thisCat = hs.categories.String{iCat};
  jS = hs.jobset;
  kitNewCategory(jS,'channel',hs.selChan,'category',thisCat, ...
      'showExisting',hs.showExistingEdit.Value)
end

function delCB(hObj,event)
  hs = handles;
  iCat = hs.categories.Value;
  thisCat = hs.categories.String{iCat};
  
  % Ask if sure.
  outcome = selectbox({'Delete','Cancel'},['Delete category ' thisCat '?']);
  if outcome == 2
      return
  end
  
  % Delete data.
  hs.categories.String(iCat) = [];
  jobs = kitLoadAllJobs(hs.jobset);
  for iJob = 1:length(jobs)
      jobs{iJob}.categories = rmfield(jobs{iJob}.categories,thisCat);
      jobs{iJob} = kitSaveJob(jobs{iJob});
  end
  hs.categories.Value = 1;
  handles = hs;
  
end

function showCatsCB(hObj,event)
  hs = handles;
  jS = hs.jobset;
  chan = hs.showChan;
  kitShowCategories(jS,chan);
end

function closeCB(hObj,event)
  uiresume(gcf);
end

function errorbox(msg)
    h=msgbox(msg,'Error','Error','modal');
    uiwait(h);
end

function outcome = selectbox(opts,msg)
    
  nopts = length(opts);  
  
  % Create figure.
  figw = 40;
  figh = 5;
  f = figure('Resize','off','Units','characters','Position',[100 35 figw figh]);
  f.DockControls = 'off';
  f.MenuBar = 'none';
  f.Name = 'Are you sure?';
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
  label(f,msg,[x y labw h],smallfont);
  y = y-lh;
  btnw = (w-2*ddx)/nopts; btnh = 2;
  for i=1:nopts
    button(f,opts{i},[x y btnw btnh],@choiceCB,smallfont);
    x = x+(btnw+ddx);
  end
  
  movegui(f,'center');
  uiwait(f);
    
  function choiceCB(hObj,event)
    outcome = find(cellfun(@(x) strcmp(x,hObj.String),opts));
    uiresume(f);
    close(f);
  end
  
end

end