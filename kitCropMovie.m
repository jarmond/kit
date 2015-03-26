function [crop,cropSize]=kitCropMovie(movieFileName,zPlane)
%KITCROPMOVIE Presents movie projection for user to crop
%
%  SYNOPSIS crop=kitCropMovie(movieFileName)
%
%  INPUT movieFileName: Filename of movie to crop.
%        zPlane: Show single z-plane (optional).
%
%  OUTPUT crop: Vector of crop coordinates as [xmin,ymin,xmax,ymax]
%
% Copyright (c) 2012 Jonathan W. Armond

if nargin<2
  zPlane=[];
end

% Open movie.
[metadata,reader] = kitOpenMovie(movieFileName);

% Load subset of frames and overlay.
nImages = 10;
frameList = unique(round(1:metadata.nFrames/nImages:metadata.nFrames));

% Loop to create max projections.
maxMergeChannels = 3;
rgbImg = zeros([metadata.frameSize(1:2), 3]);
mapChan = [2 1 3];
for c=1:min([maxMergeChannels, metadata.nChannels])
  maxProj = zeros(metadata.frameSize(1:2));
  for f=1:length(frameList)
    % Read stack.
    if isempty(zPlane)
      img = kitReadImageStack(reader, metadata, frameList(f), c);
    else
      img = kitReadImagePlane(reader, metadata, frameList(f), c, zPlane);
    end
    z = size(img,3);
    % Average max projection of each frame.
    maxProj = maxProj + max(img,[],3);
  end
  % Normalize.
  maxProj = maxProj/length(frameList);

  % Merge into RGB image.
  lb = splitModes(maxProj);
  if isempty(lb)
    % Can't identify background. Use sensible default.
    lb = 0.75;
  end
  slim = stretchlim(maxProj,[lb 1]);
  rgbImg(:,:,mapChan(c)) = imadjust(maxProj,slim,[]);
end

% Show image with crop tool.
crop = [];
p = [];
f = figure;
imshow(rgbImg);
title('Draw ROI rectangles. Double click to record each.');
finBtn = uicontrol('Style', 'pushbutton', 'String', 'Finish','FontSize',14,...
        'Position', [20 20 80 30],...
        'Callback', 'close(gcf)');
addBtn = uicontrol('Style', 'pushbutton', 'String', 'Add ROI','FontSize',14,...
        'Position', [120 20 80 30],...
        'Callback', @addROI_cb);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
while ishghandle(f)
  uiwait(gcf);
  if ~isempty(p)
    rectangle('Position',p,'LineWidth',2,'EdgeColor','y');
    crop = [crop; p];
    p = [];
  end
end

% Process ROIs into crop rectangles.
sz = size(rgbImg);
if isempty(crop)
  % If no ROI selected use whole image.
  crop = [1 1 sz(1:2)];
end

for i=1:size(crop,1)
  cropImg = imcrop(rgbImg,crop(i,:));
  sz = size(cropImg);
  cropSize(i,:) = [sz(1:2) z];
end

function addROI_cb(hObj,eventdata,handles)
  set(addBtn,'Enable','off');
  set(finBtn,'Enable','off');
  h = imrect(gca,'PositionConstraintFcn',fcn);
  p = wait(h);
  delete(h);
  set(addBtn,'Enable','on');
  set(finBtn,'Enable','on');
  uiresume(gcf);
end

end