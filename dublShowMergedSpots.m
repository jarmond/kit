function mergedImg = dublShowMergedSpots(job,varargin)
%
%
%
% Copyright (c) 2016 C. A. Smith

% default options
opts.centreChan = 2;
opts.filter = 0;
opts.imageChans = [1 2];
opts.subpixelate = 1;
opts.boxRange = 0.5; % in microns, box will have sides 2*boxRange
% process user's options
opts = processOptions(opts,varargin{:});

% open the movie and get the image data
[~, reader] = kitOpenMovie(fullfile(job.movieDirectory,job.movie));
for iChan = opts.imageChans
  movie{iChan} = kitReadWholeMovie(reader,job.metadata,iChan,job.crop,0,1);
end
% get initCoord, trackList and sisterList for the centre channel
initCoord = job.dataStruct{opts.centreChan}.initCoord;
trackList = job.dataStruct{opts.centreChan}.trackList;
sisterList = job.dataStruct{opts.centreChan}.sisterList;

% get some metadata
nCoords = 2 + job.metadata.is3D;
nSisters = length(sisterList);
nFrames = job.metadata.nFrames;
pixelSize = job.metadata.pixelSize;
chrShift = job.options.chrShift.result;

% calculate box-size
for iCoord = 1:nCoords
    boxRangePix = opts.boxRange/pixelSize(iCoord);
    boxRangePix = round(boxRangePix);%*opts.subpixelate;
    boxSize{iCoord} = -boxRangePix:boxRangePix;
end

% produce a list of trackIDs
trackIDs = sisterList(1).trackPairs(:,1:2);

for iChan = opts.imageChans
  % produce empty image for merges
  if nCoords==3
    mergedImg{iChan} = nan([[length(boxSize{1})...
        length(boxSize{2})...
        length(boxSize{3})]*opts.subpixelate...
        nSisters*2*nFrames]);
  else
    mergedImg{iChan} = nan([[length(boxSize{1})...
        length(boxSize{2})]*opts.subpixelate...
        nSisters*2*nFrames]);
  end
end

% start counter
counter=1;
% start outputting progress
prog = kitProgress(0);

% loop over sister pairs, each sister, then frames
for iSisPair = 1:nSisters
  for iSis = 1:2
    for iFrame = 1:nFrames
      
      % get trackID, then spotID
      trackID = trackIDs(iSisPair,iSis);
      spotID = trackList(trackID).featIndx(iFrame);
      
      % get pixel coordinates for reference channel
      try
      coords = initCoord(iFrame).allCoordPix(spotID,1:3);
      coords = coords([2 1 3]);
      end
      
      % with an interpolation of opts.subpixelate pixels, find nearest
      % sub-pixel location
      centrePixel = round(coords);%*opts.subpixelate);
      
      % get cropped images for each channel
      for iChan = opts.imageChans
        if nCoords==3
%           if opts.subpixelate > 1
%             tempImg = subPixelateImg(movie{iChan}(:,:,:,iFrame),opts.subpixelate);
%           else
%             tempImg = movie{iChan}(:,:,:,iFrame);
%           end
%             try
%           mergedImg{iChan}(:,:,:,counter) = tempImg(...
%             boxSize{1}+centrePixel(1),boxSize{2}+centrePixel(2),boxSize{3}+centrePixel(3));
%             catch
%                 qq=1;
%             end
        try
          tempImg = movie{iChan}(boxSize{1}+centrePixel(1),boxSize{2}+centrePixel(2),...
              boxSize{3}+centrePixel(3),iFrame);
          if opts.subpixelate > 1
            mergedImg{iChan}(:,:,:,counter) = subPixelateImg(tempImg,opts.subpixelate);
          end
        end

        else
          mergedImg{iChan}(:,:,counter) = movie{iChan}(...
            boxSize{1}+centrePixel(1),boxSize{2}+centrePixel(2),iFrame);
        end
      end
      
      % updated progress and increase counter
      prog = kitProgress(counter/(nSisters*2*nFrames),prog);
      counter = counter+1;
      
    end
  end
end
mergedImg = [];
    
end

%% Subfunctions

function newImg = subPixelateImg(img,factor)

imgSize = size(img);
nDims = length(imgSize);

newImgSize = factor*imgSize;
newImg = nan(newImgSize);

switch nDims

    case 3
        for i = 1:imgSize(1)
            for j = 1:imgSize(2)
                for k = 1:imgSize(3)
                    newImg(factor*(i-1)+1:factor*(i-1)+factor, ...
                        factor*(j-1)+1:factor*(j-1)+factor, ...
                        factor*(k-1)+1:factor*(k-1)+factor) = img(i,j,k);
                end
            end
        end
    case 2
        for i = 1:imgSize(1)
            for j = 1:imgSize(2)
                newImg(factor*(i-1)+1:factor*(i-1)+factor, ...
                    factor*(j-1)+1:factor*(j-1)+factor) = img(i,j);
            end
        end
        
end

end