function showSpots(img,spots,pixels)
% Display spots on top of image.

if nargin<3
  pixels=0;
end

if pixels
  % If 3D image, max project.
  img = max(img,[],3);

  % Make 3 layers out of original image (normalized).
  img = img/max(img(:));
  img = repmat(img,[1 1 3]);

  % Show candidate pixels in red.
  for i=1:size(spots,1)
    img(spots(i,1),spots(i,2),:) = [0 0 1];
  end

  % Plot image.
  imshow(img);
  drawnow;
else
  imshow(max(img,[],3));
  hold on;
  plot(spots(:,2),spots(:,1),'rx');
end
