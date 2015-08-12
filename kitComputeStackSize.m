function stackSize=kitComputeStackSize(cropRect,frameSize)
% KITCOMPUTESTACKSIZE Compute pixel size after cropping

if isempty(cropRect)
    stackSize = frameSize(1:3);
else
  cropRect = cropRect + 0.5; % move to pixel coordinates.
  sy=min(frameSize(2),ceil(cropRect(3)+cropRect(1)))-floor(cropRect(1));
  sx=min(frameSize(1),ceil(cropRect(4)+cropRect(2)))-floor(cropRect(2));
  stackSize = [sx,sy,frameSize(3)];
end
