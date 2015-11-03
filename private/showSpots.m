function showSpots(img,spots,pixels,transpose)
% Display spots on top of image. Spots can be cell array then multiple colours used in non-pixel mode.

if nargin<3
  pixels=0; % Set to rgb triple for colour.
end
if nargin<4
  transpose=0;
end

if transpose
  tx=2; ty=1;
else
  tx=1; ty=2;
end

h=figure(1);
if ((isscalar(pixels) && pixels) || ~isscalar(pixels)) && ~iscell(spots)
  % If 3D image, max project.
  img = max(img,[],3);

  % Make 3 layers out of original image (normalized).
  img = img/max(img(:));
  img = repmat(img,[1 1 3]);

  % Show pixels.
  for i=1:size(spots,1)
    img(spots(i,1),spots(i,2),:) = [0 0 1];
  end

  % Plot image.
  imshow(img,[]);
  drawnow;
else
  if ~isempty(spots) && ~iscell(spots)
    spots = {spots};
  end
  imshow(max(img,[],3),[]);
  hold on;
  for i=1:length(spots)
        plot(spots{i}(:,tx),spots{i}(:,ty),'x');
  end
  hold off;
end

% Scale up to at least 512 pixels
figpos = getpixelposition(h);
if any(figpos(3:4)<256)
  set(h,'Units','Pixels','Position',[figpos(1:2) (512/min(figpos(3:4)))*figpos(3:4)]);
end
