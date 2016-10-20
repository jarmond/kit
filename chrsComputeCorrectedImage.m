function [chan1out,chan2out,removed] = chrsComputeCorrectedImage(chan1in,chan2in,chrShift1to2,varargin)

% check input
if isempty(chan1in) || isempty(chan2in)
    error('Need image data to align.')
elseif size(chan1in) ~= size(chan2in)
    error('Image data structure needs to be consistent within each channel.')
elseif isempty(chrShift1to2)
    error('Need a chromatic shift vector from channel 1 to channel 2.')
end

% default options
opts.coords = [1 2];
opts.image = 0;
opts.interpol = 33; % best to use an odd number here
opts.pixelRes = [0.069384 0.069384 0.2]; % maybe remove the requirement for this, make the user give
opts.transpose = 0; % need to rotate chromatic shift cell if tranposing for a 2D image

% process options
opts = processOptions(opts,varargin{:});

if opts.transpose && length(opts.coords)==3
    warning('Cannot transpose across three coordinates. Will not transpose.')
    opts.transpose = 0;
end

%% Sub-pixelate each channel

chan1out = subPixelateImg(chan1in,opts.interpol);
chan2out = subPixelateImg(chan2in,opts.interpol);

%% Adjust channel 2 by chromatic shift

% convert chrShift to pixels
chrShift1to2 = chrShift1to2(1:3)./opts.pixelRes;
% if transposing required, transpose chrShift
if opts.transpose
    chrShift1to2(opts.coords) = fliplr(chrShift1to2(opts.coords));
end

% adjust chrShift to sub-pixelated image coordinates, and round to nearest
% pixel
chrShift1to2 = chrShift1to2*opts.interpol;
chrShift1to2 = round(chrShift1to2);

% find direction of shift
cSstat = chrShift1to2./abs(chrShift1to2);

% produce removed structure - (coord,[start end],channel)
removed = zeros(3,2,2);

% make adjustment for:
% first coordinate
iCoord = 1;
if iCoord <= length(opts.coords)
    actualCoord = opts.coords(iCoord);
    switch cSstat(actualCoord)
        case 1  % positive shift means add pixels to start of channel 2
            chan2out(chrShift1to2(actualCoord)+1:end,:,:) = chan2out(1:end-chrShift1to2(actualCoord),:,:);
            chan2out(1:chrShift1to2(actualCoord),:,:) = 0;
            removed(iCoord,1,2) = chrShift1to2(actualCoord);
            removed(iCoord,2,2) = -chrShift1to2(actualCoord);
        case -1 % negative shift means remove pixels from start of channel 2
            chan2out(1:abs(chrShift1to2(actualCoord)),:,:) = [];
            chan2out(end+1:end+abs(chrShift1to2(actualCoord)),:,:) = 0;
            removed(iCoord,1,2) = -chrShift1to2(actualCoord);
            removed(iCoord,2,2) = chrShift1to2(actualCoord);
    end
end

% second coordinate
iCoord = 2;
if iCoord <= length(opts.coords)
    actualCoord = opts.coords(iCoord);
    switch cSstat(actualCoord)
        case 1  % positive shift means add pixels to start of channel 2
            chan2out(:,chrShift1to2(actualCoord)+1:end,:) = chan2out(:,1:end-chrShift1to2(actualCoord),:);
            chan2out(:,1:chrShift1to2(actualCoord),:) = 0;
            removed(iCoord,1,2) = chrShift1to2(actualCoord);
            removed(iCoord,2,2) = -chrShift1to2(actualCoord);
        case -1 % negative shift means remove pixels from start of channel 2
            chan2out(:,1:abs(chrShift1to2(actualCoord)),:) = [];
            chan2out(:,end+1:end+abs(chrShift1to2(actualCoord)),:) = 0;
            removed(iCoord,1,2) = -chrShift1to2(actualCoord);
            removed(iCoord,2,2) = chrShift1to2(actualCoord);
    end
end

% third coordinate
iCoord = 3;
if iCoord <= length(opts.coords)
    actualCoord = opts.coords(iCoord);
    switch cSstat(actualCoord)
        case 1  % positive shift means add pixels to start of channel 2
            chan2out(:,:,chrShift1to2(actualCoord)+1:end) = chan2out(:,:,1:end-chrShift1to2(actualCoord));
            chan2out(:,:,1:chrShift1to2(actualCoord)) = 0;
            removed(iCoord,1,2) = chrShift1to2(actualCoord);
            removed(iCoord,2,2) = -chrShift1to2(actualCoord);
        case -1 % negative shift means remove pixels from start of channel 2
            chan2out(:,:,1:abs(chrShift1to2(actualCoord))) = [];
            chan2out(:,:,end+1:end+abs(chrShift1to2(actualCoord))) = 0;
            removed(iCoord,1,2) = -chrShift1to2(actualCoord);
            removed(iCoord,2,2) = chrShift1to2(actualCoord);
    end
end

% remove any duplication of coordinates not required for chromatic shifting
if length(size(chan1out)) == 3
    for iChan = setdiff(1:3,opts.coords)
        switch iChan
            case 1
                chan1out = chan1out(1,:,:);
                chan2out = chan2out(1,:,:);
            case 2
                chan1out = chan1out(:,1,:);
                chan2out = chan2out(:,1,:);
            case 3
                chan1out = chan1out(:,:,1);
                chan2out = chan2out(:,:,1);
        end 
    end
end

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