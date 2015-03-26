function [spots,spotsAmp]=waveletSpots(img,verbose)
% Find spots using multiscale product wavelet transform.
%
% See: J. C. Olivo-Marin, Pattern Recognition 35 (2002) 1989-1996 and
%      J.-L. Starck et al., Image Processing and Data Analysis: The
%      Multiscale approach, CUP (1998).
%
% Copyright (c) 2015 Jonathan W. Armond

tk=4; % hard threshold scale
ld=1.0; % detection threshold

% 3D image.
levels = 3;
[sx,sy,sz] = size(img);
W = zeros(sx,sy,sz,levels);
A = img / std(img(:));

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
for L=1:levels
  % Threshold is k * median absolute deviation of wavelet coefficients, over 0.67.
  t = tk * mad(reshape(W(:,:,:,L),[numel(W(:,:,:,L)), 1]),1) / 0.67;
  pass = W(:,:,:,L) >= t;
  W(:,:,:,L) = W(:,:,:,L) .* pass;

  % Jeffrey's noninformative prior.
  %t = k * mad(reshape(W(:,:,L),[numel(W(:,:,L)), 1]),1);
  %W(:,:,L) = max((W(:,:,L).^2 - t),0) ./ W(:,:,L);
end

% Compute multiscale product.
P = prod(W,4);
if verbose
  figure(3);
  subplot(1,2,1);
  imshow(log(max(P,[],3)),[]);
  title('log multiscale product');
end

% Threshold spots.
P(abs(P)<ld) = 0;
%se = strel('square',3)
%P = imerode(imdilate(P,se),se);

if verbose
  subplot(1,2,2);
  imshow(log(max(P,[],3)),[]);
  title('log thresholded multiscale product');
end

% Locate spots.
CC = bwconncomp(P,18);
spots = regionprops(CC,P,'Centroid','MeanIntensity','FilledArea');
% Remove isolated single pixels.
spots(vertcat(spots.FilledArea) == 1) = [];
spotsAmp = vertcat(spots.MeanIntensity);
spots = vertcat(spots.Centroid);

if verbose
  h = figure(4);
  imshow(max(img,[],3),[]);
  scale = 3;
  figpos = get(h,'Position');
  set(h,'Position',[figpos(1:2) figpos(3:4)*scale]);

  if ~isempty(spots)
    hold on;
    plot(spots(:,1),spots(:,2),'rx');
    hold off;
  end
  title('spots');
end
