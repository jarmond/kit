function spots=histcutSpots(img,options,dataProperties,verbose)
% Find spots using histogram mode cutoff.
%
% Based on code from MaKi, by (mostly) K. Jaqaman and J. Dorn.
%
% Copyright (c) 2015 Jonathan W. Armond

if nargin<4
  verbose = 0;
end

img = img/max(img(:));

[sx,sy,sz] = size(img);
if sz == 1
  ndims = 2;
else
  ndims = 3;
end

filters = createFilters(ndims,dataProperties);

% Filter image.
imageF = fastGauss3D(img,[],filters.signalP,filters.border,filters.signal);
background = fastGauss3D(imageF,[],filters.backgroundP,filters.border,filters.background);
imageNoise = (img-imageF).^2;
imageNoise = fastGauss3D(imageNoise,[],filters.signalP,filters.border,filters.noise);

% get local maxima from the image
bw = imregionalmax(img);
localMax1DIndx = find(bw);
if ndims == 3
  [x,y,z]=ind2sub([sx sy sz],localMax1DIndx);
  locMax=[x,y,z];
else
  [x,y]=ind2sub([sx sy],localMax1DIndx);
  locMax=[x,y];
end
if any(x>sx) || any(y>sy) || (ndims==3 && any(z>sz))
  error('invalid index');
end

% get signal strength of maxima
amp = imageF(localMax1DIndx) - background(localMax1DIndx);
noise = imageNoise(localMax1DIndx);
dark = amp./sqrt(noise);
poisson = amp./sqrt(noise./max(amp,eps));

% Filter out spots by cutting first histogram mode.
z = zeros(size(amp));
% amplitude
cutoff(1) = splitModes(amp(~isApproxEqual(amp,z,[],'absolute')),[],[],[]);
% amplitude/sqrt(noise) - dark noise
cutoff(2) = splitModes(dark(~isApproxEqual(dark,z,[],'absolute')),[],[],[]);
% amplitude/sqrt(noise/amp) - poisson
cutoff(3) = splitModes(poisson(~isApproxEqual(poisson,z,[],'absolute')),[],[],[]);

if verbose
  figure(1)
  subplot(3,1,1);
  histogram(amp);
  hold on;
  yl = ylim;
  plot(repmat(cutoff(1),[1 2]),yl,'r');
  title('amplitude');

  subplot(3,1,2);
  histogram(dark);
  hold on;
  yl = ylim;
  plot(repmat(cutoff(2),[1 2]),yl,'r');
  title('dark noise');

  subplot(3,1,3);
  histogram(poisson);
  hold on;
  yl = ylim;
  plot(repmat(cutoff(3),[1 2]),yl,'r');
  title('poisson noise');
  hold off
end

% Check if sufficient spots found in descending order of cutoff strictness.
nn(3) = sum(poisson>cutoff(3))/options.minSpotsPerFrame;
nn(2) = sum(dark>cutoff(2))/options.minSpotsPerFrame;
nn(1) = sum(amp>cutoff(1))/options.minSpotsPerFrame;
if nn(3) > 1
  passIdx = poisson>cutoff(3);
elseif nn(2) > 1
  passIdx = dark>cutoff(2);
else
  passIdx = amp>cutoff(1);
end

% keep these spots
spots = locMax(passIdx,:);
%candsBg = background(localMax1DIndx(passIdx));
%candsAmp = imageF(localMax1DIndx(passIdx));

if verbose
  h = figure(2);
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
