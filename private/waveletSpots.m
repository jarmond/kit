function [spots,spotsAmp,ld]=waveletSpots(img,opts,dataProperties)
% Find spots using multiscale product wavelet transform.
%
% See: J. C. Olivo-Marin, Pattern Recognition 35 (2002) 1989-1996 and
%      J.-L. Starck et al., Image Processing and Data Analysis: The
%      Multiscale approach, CUP (1998).
%
% Copyright (c) 2015 Jonathan W. Armond

verbose = opts.debug.showWavelet;
tk = opts.waveletLevelThresh; % threshold scale for local MAD thresholding
levels = opts.waveletNumLevels;  % number of wavelet levels
localmad = opts.waveletLocalMAD; % locally estimated MAD
backsub = opts.waveletBackSub;  % background subtraction
prefilter = opts.waveletPrefilter; % prefilter with Gaussian.

% 3D image.
[sx,sy,sz] = size(img);
W = zeros(sx,sy,sz,levels);

% Normalization and background subtraction.
if backsub
  A = img - imfilter(img,fspecial('gaussian',round(min([sx,sy])/2)),'symmetric');
else
  A = img;
end
A = A/max(A(:));

% Scaling function: B3-spline.
scalingFn = [1,4,6,4,1]/16;
h = cell(levels,1);

% Precompute filters.
for j=1:levels
  if j>1
    % Augment filter. Insert nz zeros between each pair of taps.
    nz = 2^(j-1) - 1;
    h{j} = [scalingFn; zeros(nz,size(scalingFn,2))];
    h{j} = h{j}(1:end-nz);
  else
    h{j} = scalingFn;
  end
end

for L=1:levels
  % Convolve filter.
  A2 = imfilter(A,h{L}'*h{L},'symmetric','conv');

  % Compute wavelet plane.
  W(:,:,:,L) = A - A2;
  A = A2;
end

if verbose
  figure(1);
  n=levels+1;
  fig_n=ceil(sqrt(n));
  fig_m=ceil(levels/fig_n);
  for i=1:levels
    subplot(fig_m,fig_n,i);
    imshow(max(W(:,:,:,i),[],3),[])
    title(['level ' num2str(i)]);
  end
  subplot(fig_m,fig_n,n);
  imshow(max(A2,[],3),[]);
  title('final smooth');
  suptitle('wavelet planes');

  figure(2);
  R = A2 + sum(W,4);
  imshow(max(R,[],3),[]);
  title('reconstruction');
end


% Hard threshold wavelet coefficients.
for L=2:levels
  % Threshold is median + k * median absolute deviation of wavelet coefficients, over 0.67.
  WL = W(:,:,:,L);
  if localmad
    h = fspecial('average',64);
    madest = imfilter(abs(WL-imfilter(WL,h,'symmetric')),h,'symmetric');
    h = fspecial('average',16);
    bkgd=imfilter(WL,h,'symmetric');
    t = bkgd + tk * 1.4826 * madest;
  else
    t = median(WL(:)) + tk * mad(WL(:),1) * 1.4826;
  end
  WL(WL<t) = 0;
  W(:,:,:,L) = WL;
end

% Compute multiscale product on level 2..L.
P = prod(W(:,:,:,2:L),4);
if verbose
  figure(3);
  subplot(1,2,1);
  F = max(abs(P),[],3);
  imshow(log(F),[]);
  title('log multiscale product');
end

% Threshold spots.
Pt = P(P>0 & P<max(P(:))/4);
if length(Pt)>50 && max(Pt)>min(Pt)
  % Estimate histogram mode.
  [f,xi]=ksdensity(Pt);
  [~,i]=max(f);
  ld=xi(i);
  P(P<ld) = 0;
else
  ld = nan;
end

if verbose
  subplot(1,2,2);
  imshow(log(max(abs(P),[],3)),[]);
  title('log thresholded multiscale product');
end

% Locate spots.
CC = imregionalmax(P);
spots = regionprops(CC,P,'Centroid','MeanIntensity');
spotsAmp = vertcat(spots.MeanIntensity);
spots = vertcat(spots.Centroid);

% Convert to image coordinates, for compatibility with later processing.
if ~isempty(spots)
  spots(:,1:2) = fliplr(spots(:,1:2));
end

if verbose
  h = figure(4);
  imshow(max(img,[],3),[]);
  % scale = 3;
  % figpos = get(h,'Position');
  % set(h,'Position',[figpos(1:2) figpos(3:4)*scale]);

  if ~isempty(spots)
    hold on;
    plot(spots(:,2),spots(:,1),'rx');
    hold off;
  end
  title('spots');
end
