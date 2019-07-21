function kitShowCategories(jobset,channel)
% KITSHOWCATEGORIES(JOBSET) Displays a GUI to show categorisation of
% spots within a JOBSET.
%
%    KITSHOWCATEGORY(JOBSET)
%
% Copyright (c) 2019 C. A. Smith

% Define colours for spots.
handles.labcol = [ 1 , 1 , 0;
                   0 , 1 , 1;
                  0.5, 1 ,0.5;
                  0.6, 0 , 1];
           
% Get the data.
job = kitLoadAllJobs(jobset);
handles.nMovs = length(job);

% Define channel.
if nargin<2 || isempty(channel)
    channel = jobset.options.coordSystemChannel;
end

% Predefine some handles required during looping
handles.movID = 1;
handles.prevEnable = 'off';
if length(job) > 1
    handles.nextEnable = 'on';
else
    handles.nextEnable = 'off';
end
    
% Start while loop until aborted
handles.stop = 0;
while ~handles.stop
        
    % get this image ID
    iMov = handles.movID;
    
    % check whether there is any data contained within this movie
    if ~isfield(job{iMov},'dataStruct') || ~isfield(job{iMov}.dataStruct{channel},'failed') || job{iMov}.dataStruct{channel}.failed
        handles.movID = handles.movID+1;
        if handles.movID > handles.nMovs
            handles.stop = 1;
        end
        continue
    end
    % get dataStruct
    dS = job{iMov}.dataStruct{channel};
    % get number of spots
    nSpots = length(dS.initCoord.allCoord);
    handles.nSpots = nSpots;
    
    % get categories
    cG = job{iMov}.categories;
    catlabs = fieldnames(cG);
    nCats = length(catlabs);
    cats = struct;
    for iCat = 1:nCats
        cats(iCat).label = catlabs{iCat};
        cats(iCat).list = cG.(cats(iCat).label);
        cats(iCat).colour = handles.labcol(iCat,:);
    end
    
    % give nROIs to the job
    job{iMov}.nROIs = handles.nMovs;

    % show all spots - defined using tracks
    rectDims = griddedSpots(job{iMov},'channel',channel,'rawData',0);
    
    % label figure title
    figtit = sprintf('Categorised spots: Image %i%s, channel %i',job{iMov}.index,handles.nMovs,channel);
    set(gcf,'Name',figtit);
    
    % get image information
    rectPos = rectDims(:,1:2);
    rectWid = rectDims(1,3);
    handles.rectPos = rectPos;
    handles.rectWid = rectWid;
    
    % label existing categories
    hold on
    for iCat = 1:nCats
        labelCategory(gca,cats(iCat),rectDims,iCat);
    end
    
    % Buttons and labels.
    btnw = [12 7]; btnh = 2; h = 1.5;
    figpos = get(gcf,'Position');
    dx = 2.5; ddx = 1;
    % finish, next and previous
    x = figpos(3)-(btnw(1)+dx); y = dx;
    handles.finishBtn = button(gcf,'Finish',[x y btnw(1) btnh],@finishCB);
    x = x-(btnw(2)+ddx);
    handles.nextBtn = button(gcf,'Next',[x y btnw(2) btnh],@nextMovCB);
    handles.nextBtn.Enable = handles.nextEnable;
    x = x-(btnw(2)+ddx/2);
    handles.prevBtn = button(gcf,'Prev',[x y btnw(2) btnh],@prevMovCB);
    handles.prevBtn.Enable = handles.prevEnable;
    % categories key
    labw = 17.5;
    x = dx;
    handles.categories = label(gcf,'Categories key:',[x y labw h],12);
    x = x+labw; y = y-h;
    labw = 20+(3*ddx);
    t = label(gcf,'',[x y labw 2*h],12); t.BackgroundColor = [0 0 0];
    y = y+h; x = x+ddx;
    labw = (labw-3*ddx)/2;
    for iCat = 1:nCats
        if iCat == 2
            x = x+(labw+ddx);
        elseif iCat == 3
            y = y-h;
        elseif iCat == 4
            x = x-(labw+ddx);
        end
        handles.catLab = label(gcf,cats(iCat).label,[x y labw h],10);
        handles.catLab.ForegroundColor = handles.labcol(iCat,:);
        handles.catLab.BackgroundColor = [0 0 0];
        if ismember(iCat,2:3)
            handles.catLab.HorizontalAlignment = 'right';
        end
    end
    
    % GUI now set up, wait for user
    uiwait(gcf);
    % close the figure for the next movie
    close(gcf);

end

%% Callback functions

function prevMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID-1;
  if handles.movID == 1
    handles.prevEnable = 'off';
  end
  handles.nextEnable = 'on';
  handles.nextChan = 0;
  % continue the function
  uiresume(gcf);
end

function nextMovCB(hObj,event)
  % update the handles
  handles.movID = handles.movID+1;
  handles.prevEnable = 'on';
  if handles.movID == handles.nMovs
    handles.nextEnable = 'off';
  end
  handles.nextChan = 0;
  % continue the function
  uiresume(gcf);
end

function finishCB(hObj,event)
  % force stop
  handles.stop = 1;
  % continue the function
  uiresume(gcf);
end

end