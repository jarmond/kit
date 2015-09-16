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

[rgbImg,zPlanes] = kitMovieProj(movieFileName,zPlane);

% Show image with crop tool.
crop = [];
f = gcf;
p = [];
title('Draw ROI rectangles. Double click to record each.');
finBtn = uicontrol('Style', 'pushbutton', 'String', 'Finish','FontSize',14,...
        'Position', [20 10 80 30],...
        'Callback', 'close(gcf)');
addBtn = uicontrol('Style', 'pushbutton', 'String', 'Add ROI','FontSize',14,...
        'Position', [120 10 80 30],...
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

crop = round(crop);

% Process ROIs into crop rectangles.
sz = size(rgbImg);
if isempty(crop)
  % If no ROI selected use whole image.
  crop = [1 1 sz(1:2)];
end

for i=1:size(crop,1)
  cropImg = imcrop(rgbImg,crop(i,:));
  sz = size(cropImg);
  cropSize(i,:) = [sz(1:2) zPlanes];
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
