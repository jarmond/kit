function kitPlayMovie(filename)
% KITPLAYMOVIE Play movie file

if nargin<1
  [filename,pathname] = uigetfile(...
    kitSupportedFormats(1),'Locate movie to play');
else
  [pathname,filename,ext] = fileparts(filename);
  filename = [filename ext];
end

% Open movie.
[md,reader] = kitOpenMovie(fullfile(pathname,filename));
opts.saturate = [1 99.9];
if size(opts.saturate,1)<md.nChannels
  opts.saturate = repmat(opts.saturate,[md.nChannels,1]);
end
mapChans = [2 1 3]; % Green, red, blue
maxMergeChannels = 3;
dt = md.frameTime(1,2)-md.frameTime(1,1);
len = md.frameTime(1,end)-md.frameTime(1,1);
fprintf('dt = %.2fs  length = %.1f min\n', dt, len/60);

figure;
clf;

for i=1:md.nFrames
  rgbImg = zeros([md.frameSize(1:2), 3]);
  for c=1:min(md.nChannels, maxMergeChannels)
    % Read stack.
    img = kitReadImageStack(reader, md, i, c, [], 0);

    % Max project.
    img = max(img, [], 3);

    % First frame defines contrast stretch.
    if i==1
      irange(c,:)=stretchlim(img,opts.saturate(c,:)/100);
    end

    % Contrast stretch.
    rgbImg(:,:,mapChans(c)) = imadjust(img, irange(c,:), []);
  end

  imshow(rgbImg);
  axis image;
  drawnow;
end