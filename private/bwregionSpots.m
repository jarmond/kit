function spots=bwregionSpots(movie,options)


nFrames = size(movie,4);
spots = cell(nFrames,1);
for i=1:nFrames
    img = movie(:,:,:,i);

    % Threshold Otsu and BW.
    level = graythresh(img);
    img = img>=level;

    % Erode then dilate.
    img = imerode(img,strel('sphere',1));
    img = imdilate(img,strel('sphere',3));

    % Connected regions.
    s = regionprops(img,'centroid');
    centroid = cat(1, s.Centroid);
    spots{i} = centroid(:,[2 1 3]);
end
