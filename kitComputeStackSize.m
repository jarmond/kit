function stackSize=kitComputeStackSize(cropRect,frameSize)
% KITCOMPUTESTACKSIZE Compute pixel size after cropping

if isempty(cropRect)
    stackSize = frameSize(1:3);
else
    sy=round(cropRect(3)+cropRect(1))-round(cropRect(1))+1;
    sx=round(cropRect(4)+cropRect(2))-round(cropRect(2))+1;
    stackSize = [sx,sy,frameSize(3)];
end
