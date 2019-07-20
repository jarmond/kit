function [chan1out,chan2out,removed] = chrsComputeCorrectedImage(chan1in,chan2in,chrShift1to2,varargin)
%CHRSCOMPUTECORRECTEDIMAGE Produces a sub-pixelated two-colour image to
%    allow sub-pixel chromatic shifted image production.
%
%    CHRSCOMPUTECORRECTEDIMAGE(CHAN1IN,CHAN2IN,CHRSHIFT1TO2,...) Produces
%    a combined image of two channels, CHAN1IN and CHAN2IN, shifting by
%    chromatic shift defined from channel 1 to 2, CHRSHIFT1TO2. Supply
%    options as string/value pairs following JOB.
%
%    Options, defaults in {}:-
%
%    coords: {[1 2]} or two numbers from 1 to 3. Coordinates in which the
%         printed image will be orientated.
%
%    image: {0} or 1. Whether or not to print the resulting image to a
%         figure.
%
%    pixelSize: {[0.0694 0.0694 0.2]}, or 3-coordinate vector. Size of each
%         pixel in the provided images, in microns.
%
%    subpixelate: {9} or other odd number. Factor by which to sub-pixelate
%         images. It is recommended that this number be odd.
%
%
% Copyright (c) 2016 C. A. Smith

% check input
if isempty(chan1in) || isempty(chan2in)
    error('Need image data to correct.')
elseif size(chan1in) ~= size(chan2in)
    error('Image data structure needs to be consistent between channels.')
elseif isempty(chrShift1to2)
    error('Need a chromatic shift vector from channel 1 to channel 2.')
end

% default options
opts.coords = [1 2];
opts.image = 0;
opts.pixelSize = [0.069384 0.069384 0.2]; % maybe remove the requirement for this, make the user give
opts.subpixelate = 9; % best to use an odd number here

% process options
opts = processOptions(opts,varargin{:});

%% Sub-pixelate each channel

chan1out = subPixelateImg(chan1in,opts.subpixelate);
chan2out = subPixelateImg(chan2in,opts.subpixelate);

%% Adjust channel 2 by chromatic shift

% convert chrShift to pixels
chrShift1to2 = chrShift1to2(1:3)./opts.pixelSize;

% adjust chrShift to sub-pixelated image coordinates, and round to nearest
% sub-pixel
chrShift1to2 = chrShift1to2*opts.subpixelate;
chrShift1to2 = round(chrShift1to2);

% find direction of shift
cSstat = sign(chrShift1to2);
chrShift1to2 = abs(chrShift1to2);

% produce removed structure - (coord,[start end],channel)
removed = zeros(3,2,2);

% make adjustment for:
% first coordinate
iCoord = 1;
if iCoord <= length(opts.coords)
    actualCoord = opts.coords(iCoord);
    switch cSstat(actualCoord)
        case +1  % positive shift means correcting channel 2 in the negative direction, i.e. remove pixels from start of channel 2
            chan2out(1:chrShift1to2(actualCoord),:,:) = [];
            chan2out(end+1:end+chrShift1to2(actualCoord),:,:) = 0;
            removed(actualCoord,1,2) = chrShift1to2(actualCoord);
        case -1 % negative shift means add pixels to start of channel 2
            chan2out(chrShift1to2(actualCoord)+1:end,:,:) = chan2out(1:end-chrShift1to2(actualCoord),:,:);
            chan2out(1:chrShift1to2(actualCoord),:,:) = 0;
            removed(actualCoord,2,2) = chrShift1to2(actualCoord);
    end
end

% second coordinate
iCoord = 2;
if iCoord <= length(opts.coords)
    actualCoord = opts.coords(iCoord);
    switch cSstat(actualCoord)
        case +1  % positive shift means correcting channel 2 in the negative direction, i.e. remove pixels from start of channel 2
            chan2out(:,1:chrShift1to2(actualCoord),:) = [];
            chan2out(:,end+1:end+chrShift1to2(actualCoord),:) = 0;
            removed(actualCoord,1,2) = chrShift1to2(actualCoord);
        case -1 % negative shift means add pixels to start of channel 2
            chan2out(:,chrShift1to2(actualCoord)+1:end,:) = chan2out(:,1:end-chrShift1to2(actualCoord),:);
            chan2out(:,1:chrShift1to2(actualCoord),:) = 0;
            removed(actualCoord,2,2) = chrShift1to2(actualCoord);
    end
end

% third coordinate
iCoord = 3;
if iCoord <= length(opts.coords)
    actualCoord = opts.coords(iCoord);
    switch cSstat(actualCoord)
        case +1  % positive shift means correcting channel 2 in the negative direction, i.e. remove pixels from start of channel 2
            chan2out(:,:,1:chrShift1to2(actualCoord)) = [];
            chan2out(:,:,end+1:end+chrShift1to2(actualCoord)) = 0;
            removed(actualCoord,1,2) = chrShift1to2(actualCoord);
        case -1 % negative shift means add pixels to start of channel 2
            chan2out(:,:,chrShift1to2(actualCoord)+1:end) = chan2out(:,:,1:end-chrShift1to2(actualCoord));
            chan2out(:,:,1:chrShift1to2(actualCoord)) = 0;
            removed(actualCoord,2,2) = chrShift1to2(actualCoord);
    end
end

chan1out = squeeze(chan1out);
chan2out = squeeze(chan2out);

end