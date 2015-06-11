function [spots,spotsAmp]=waveletSpots(img,varargin)
% Find spots using multiscale product wavelet transform.
%
% See: J. C. Olivo-Marin, Pattern Recognition 35 (2002) 1989-1996 and
%      J.-L. Starck et al., Image Processing and Data Analysis: The
%      Multiscale approach, CUP (1998).
%
% Copyright (c) 2015 Jonathan W. Armond

options.verbose = 0;
options.tk=2; % hard threshold scale
options.ld=0.01; % detection threshold
options.opsz=2; % should be about PSF sigma.
options=processOptions(options,varargin{:});

% 3D image.
levels = 3;
[sx,sy,sz] = size(img);
W = zeros(sx,sy,sz,levels);
A = img / max(img(:));

% Scaling function: B3-spline.
scalingFn = [1,4,6,4,1]/16;
%scalingFn2 = scalingFn' * scalingFn; % 2D filter.
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
  % Mirror symmetric cols.
  A = padarray(A,[0,sy],'symmetric');
  A2 = zeros(size(A));

  % Compute separable wavelet, row-by-row, col-by-col, z-col-by-z-col.
  % Rows.
  for k=1:sz
    for i=1:sx
      A2(i,:,k) = conv(A(i,:,k),h{L},'same');
    end
  end
  % Trim edges.
  A = A(:,sy+1:2*sy,:);

  % Mirror symmetric rows.
  A2 = padarray(A2(:,sy+1:2*sy,:),[sx,0],'symmetric');
  % Cols.
  for k=1:sz
    for j=1:sy
      A2(:,j,k) = conv(A2(:,j,k),h{L}','same');
    end
  end
  % Trim edges.
  A2 = A2(sx+1:2*sx,:,:);

  % Mirror in Z.
  A2 = padarray(A2,[0,0,sz],'symmetric');
  % Z-cols.
  A3 = A2;
  A2 = shiftdim(A2,2);
  for j=1:sy
    for i=1:sx
      % A2(i,j,:) = conv(squeeze(A2(i,j,:)),h{L}','same');
      A2(:,i,j) = conv(A2(:,i,j),h{L}','same');
    end
  end
  A2 = shiftdim(A2,1);

  % Trim edges.
  A2 = A2(:,:,sz+1:2*sz);

  % Compute wavelet plane.
  W(:,:,:,L) = A - A2;
  A = A2;
end

if options.verbose
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
for L=1:levels
  % Threshold is background + k * local median absolute deviation of wavelet
  % coefficients, over 0.67.
  WL = W(:,:,:,L);
  bkgd = imgaussfilt3(WL,20);
  madest = imgaussfilt3(abs(WL-bkgd),20);
  t = bkgd + options.tk * 1.4826 * madest;
  pass = WL >= t;
end

% Compute multiscale product on level 2..L.
P = prod(W(:,:,:,2:levels),4);
if options.verbose
  figure(3);
  subplot(1,2,1);
  imshow(log(max(P,[],3)),[]);
  title('log multiscale product');
end

% Threshold spots.
P(abs(P)<options.ld) = 0;
se = strel('disk',options.opsz);
P = imclose(P,se);

if options.verbose
  subplot(1,2,2);
  imshow(log(max(P,[],3)),[]);
  title('log thresholded multiscale product');
end

% Locate spots.
CC = bwconncomp(P);
spots = regionprops(CC,P,'Centroid','MeanIntensity','FilledArea');
% Remove isolated single pixels.
spots(vertcat(spots.FilledArea) == 1) = [];
spotsAmp = vertcat(spots.MeanIntensity);
spots = vertcat(spots.Centroid);

% Convert to image coordinates, for compatibility with later processing.
if ~isempty(spots)
  spots = spots(:,[2 1 3]);
end

if options.verbose
  h = figure(4);
  imshow(max(img,[],3),[]);
  scale = 3;
  figpos = get(h,'Position');
  set(h,'Position',[figpos(1:2) figpos(3:4)*scale]);

  if ~isempty(spots)
    hold on;
    plot(spots(:,2),spots(:,1),'rx');
    hold off;
  end
  title('spots');
end
