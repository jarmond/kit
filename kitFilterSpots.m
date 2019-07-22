function kitFilterSpots(jobset)
% KITFILTERSPOTS(JOBSET) Displays a GUI to allow removal of erroneous
% spots from a JOBSET.
%
% Copyright (c) 2018 C. A. Smith

% Define colours for rectangles.
handles.col = [1 0   0;...
               0 0.75 0];
           
% Get the data.
job = kitLoadAllJobs(jobset);
handles.nMovs = length(job);
handles.chans = find(cellfun(@(x) ~strcmp(x,'none'),jobset.options.spotMode));

% Get basic information from the jobset.
opts = jobset.options;
handles.chanID = opts.coordSystemChannel;

nFrames = job{1}.metadata.nFrames;
if nFrames > 1
  error('This function cannot yet be used for movies.')
end

% Predefine some handles required during looping.
handles.movID = 1;
handles.prevEnable = 'off';
if length(job) > 1
    handles.nextEnable = 'on';
else
    handles.nextEnable = 'off';
end
handles.nextChan = 0;

% Start progress.
prog = kitProgress(0);

% Start while loop until aborted
handles.stop = 0;
while ~handles.stop
    
    iMov = handles.movID;
    iChan = handles.chanID;

    % Check whether there is any data contained within this movie.
    if ~isfield(job{iMov},'dataStruct') || ~isfield(job{iMov}.dataStruct{iChan},'failed') || job{iMov}.dataStruct{iChan}.failed
        if handles.movID > handles.nMovs
            handles.stop = 1;
        end
        continue
    end
    
    dS = job{iMov}.dataStruct{iChan};
    iC = dS.initCoord;
    % Get IDs for all spots not previously filtered.
    nonNaNs = find(~isnan(iC.allCoord(:,1)));
    
    % Back up the full initCoord and spotInt.
    if isfield(dS,'rawData')
      iC = dS.rawData.initCoord;
      if isfield(dS.rawData,'spotInt')
        sI = dS.rawData.spotInt;
      end
    elseif isfield(dS,'spotInt')
        sI = dS.spotInt;
    end
    raw.initCoord = iC;
    if exist('sI','var')
        raw.spotInt = sI;
    end
    dS.rawData = raw;
    job{iMov}.dataStruct{iChan} = dS;
    
    % Get and store number of spots.
    nSpots = length(dS.initCoord.allCoord);
    handles.nSpots = nSpots;

    % Show all spots, numbered by track.
    rectDims = griddedSpots(job{iMov},'channel',handles.chanID);
    
    % Get image information
    rectPos = rectDims(:,1:2);
    rectWid = rectDims(1,3);
    handles.rectPos = rectPos;
    handles.rectWid = rectWid;
    if handles.nextChan
      handles.keep = ismember(handles.keep,nonNaNs);
    else
      handles.keep = ismember(1:nSpots,nonNaNs);
    end
    
    % Reset channel information if necessary.
    handles.nextChan = 0;
    
    % Draw rectangles.
    hold on
    for iSpot = 1:nSpots
        % Get the colour for this spot.
        keep = handles.keep(iSpot);
        icol = handles.col(keep+1,:);
        
        rectangle('Position',[rectPos(iSpot,:)-0.5 rectWid rectWid],...
            'EdgeColor',icol,'LineWidth',3);
    end
    
    % Buttons and labels.
    btnw = [12 7]; btnh = 2; h = 1.5;
    figpos = get(gcf,'Position');
    dx = 2.5; ddx = 1;
    % Label at top left for instructions.
    x = dx; y = figpos(4)-(btnh+ddx);
    labw = 60;
    handles.instructions = label(gcf,'Click on spots to keep (green) or remove (red).',[x y labw h],12);
    % Finish, next and previous.
    x = figpos(3)-(btnw(1)+dx);
    handles.finishBtn = button(gcf,'Finish',[x y btnw(1) btnh],@finishCB);
    x = x-(btnw(2)+ddx);
    handles.nextBtn = button(gcf,'Next',[x y btnw(2) btnh],@nextMovCB);
    handles.nextBtn.Enable = handles.nextEnable;
    x = x-(btnw(2)+ddx/2);
    handles.prevBtn = button(gcf,'Prev',[x y btnw(2) btnh],@prevMovCB);
    handles.prevBtn.Enable = handles.prevEnable;
    % Invert.
    x = figpos(3)-(btnw(1)+dx); y = dx;
%     handles.deselectBtn = button(gcf,'Deselect all',[x y btnw(1) btnh],@deselectAllCB);
    handles.invertBtn = button(gcf,'Invert all',[x y btnw(1) btnh],@invertCB);
    % Next channel.
    x = x-(btnw(1)+dx); y = dx;
    handles.nextChanBtn = button(gcf,'Next chan',[x y btnw(1) btnh],@nextChanCB);
    
    % Set up remove environment.
    set(get(gca,'Children'),'ButtonDownFcn',@rmvCB);
    
    % GUI now set up, wait for user...
    uiwait(gcf);
    
    if handles.nextChan
        % Close the figure for the next channel.
        close(gcf);
        continue
    end
    
    % Get the spots requiring removal, and check for any not in the list.
    rmvList = find(~handles.keep);
    incorrect = setdiff(rmvList,1:nSpots);
    if ~isempty(incorrect)
        warning('The following selected spots do not exist: %s. Will ignore.',num2str(incorrect));
    end

    % Push all removed spots to NaNs.
    if ~isempty(rmvList)
      % initCoord.
      for jChan = handles.chans
        dS = job{iMov}.dataStruct{jChan};
        iC = dS.initCoord;
        iC.allCoord(rmvList,:) = NaN;
        iC.allCoordPix(rmvList,:) = NaN;
        iC.nSpots = sum(~isnan(iC.allCoord(:,1)));
        iC.amp(rmvList,:) = NaN;
        iC.bg(rmvList,:) = NaN;
        
        % Back up results.
        dS.initCoord = iC;
        job{iMov}.dataStruct{jChan} = dS;
      end
      % spotInt.  
      for jChan = find(job{iMov}.options.intensity.execute)
        dS = job{iMov}.dataStruct{jChan};
        sI = dS.spotInt;
        sI.intensity_mean(rmvList,:) = NaN;
        sI.intensity_median(rmvList,:) = NaN;
        sI.intensity_min(rmvList,:) = NaN;
        sI.intensity_max(rmvList,:) = NaN;
        sI.intensity_ratio(rmvList,:) = NaN;
        
        % Back up results.
        dS.spotInt = sI;
        job{iMov}.dataStruct{jChan} = dS;
      end
      
    end
    
    % Save results.
    job{iMov} = kitSaveJob(job{iMov});

    % Update progress.
    prog = kitProgress(iMov/handles.nMovs,prog);
    
    % Close the figure for the next movie.
    close(gcf);

end

kitLog('Manual filtering complete.');

% Re-run plane fitting for jobset with new filtered data.
kitLog('Re-fitting planes to filtered data.');
kitRunJob(jobset,'existing',1,'tasks',[2 6]);

%% Callback functions

function rmvCB(hObj,event)
  
  % get the position of the click
  pos=get(gca,'CurrentPoint');
  xpos = pos(1,1); ypos = pos(1,2);

  % get all positions
  allPos = handles.rectPos;
  
  % get candidates using click's row position
  diffs = xpos-allPos(:,1);
  diffs(diffs<0) = NaN;
  xidx = find(diffs == nanmin(diffs));
  % get candidates using click's column position
  diffs = ypos-allPos(:,2);
  diffs(diffs<0) = NaN;
  yidx = find(diffs == nanmin(diffs));
  % get the common candidate
  idx = intersect(xidx,yidx);

  % if a click is made elsewhere, remind user how to select images
  if isempty(idx)
    handles.instructions.String = 'Click on the images to select/deselect.';
    return
  end
  
  % get the colour for this spot
  keepStat = ~handles.keep(idx);
  icol = handles.col(keepStat+1,:);

  % draw the rectangle
  rectangle('Position',[handles.rectPos(idx,:)-0.5 handles.rectWid handles.rectWid],...
      'EdgeColor',icol,'LineWidth',3);

  handles.keep(idx) = keepStat;
  
end

function deselectAllCB(hObj,event)
  hs = handles;
  % force all stops to be ignored
  handles.keep(1,:) = 0;
  for i = 1:hs.nSpots
    % draw the rectangle
    rectangle('Position',[hs.rectPos(i,:)-0.5 hs.rectWid hs.rectWid],...
        'EdgeColor',hs.col(1,:),'LineWidth',3);
  end
end

function invertCB(hObj,event)
  hs = handles;
  % force all stops to be ignored
  handles.keep = ~handles.keep;
  for i = 1:hs.nSpots
    
    % get the colour for this spot
    keepStat = handles.keep(i);
    jcol = handles.col(keepStat+1,:);
    
    % draw the rectangle
    rectangle('Position',[hs.rectPos(i,:)-0.5 hs.rectWid hs.rectWid],...
        'EdgeColor',jcol,'LineWidth',3);
  end
end

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

function nextChanCB(hObj,event)
  % update the handles
  idx = find(handles.chanID==handles.chans);
  if idx == length(handles.chans)
    handles.chanID = handles.chans(1);
  else
    handles.chanID = handles.chans(idx+1);
  end
  handles.nextChan = 1;
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


    