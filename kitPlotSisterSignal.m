function x=kitPlotSisterSignal(job,idx,sigchannel,varargin)
% KITPLOTSISTERSIGNAL Plot sister tracks and signal
%
%    KITPLOTSISTERSIGNAL(JOB,IDX,SIGCHANNEL,...) Plots sister with index IDX
%    against signal intensity in channel SIGCHANNEL. Supply options as
%    string/value pairs following IDX.
%
%    Options, defaults in {}:-
%
%    channel: {1} or number. Channel for which tracks are plotted.
%
% Copyright(c) Jonathan W. Armond 2014

if nargin<3
    error('Must supply JOB, IDX and SIGCHANNEL');
end

% Set defaults
opts.channel = 1;

% Process options
opts = processOptions(opts, varargin{:});

ds = job.dataStruct{opts.channel};
sisterList = ds.sisterList;
trackList = ds.trackList;
trackInt = ds.trackInt;
t = job.metadata.frameTime;
pair = sisterList(1).trackPairs(idx,1:2);

figure;
subplot(2,1,1);
plot(t,sisterList(idx).coords1(:,1),'b-', ...
     t,sisterList(idx).coords2(:,1),'r-')
hold on;
sc=sisterList(idx).coords1(:,1);
sc(trackList(pair(1)).direction~=+1)=nan;
plot(t,sc,'g-','linewidth',2);
sc=sisterList(idx).coords2(:,1);
sc(trackList(pair(2)).direction~=+1)=nan;
plot(t,sc,'g-','linewidth',2);
sc=sisterList(idx).coords1(:,1);
sc(trackList(pair(1)).direction~=-1)=nan;
plot(t,sc,'m-','linewidth',2);
sc=sisterList(idx).coords2(:,1);
sc(trackList(pair(2)).direction~=-1)=nan;
plot(t,sc,'m-','linewidth',2);

title('KT position; green P, magenta AP');
xlim([min(t) max(t)])
subplot(2,1,2);
hold on;
plot(t,trackInt(pair(1)).intensity_max(:,sigchannel),'b-',...
     t,trackInt(pair(2)).intensity_max(:,sigchannel),'r-')
title('intensity');
xlim([min(t) max(t)])

% Return data.
x=[t' sisterList(idx).coords1(:,1) sisterList(idx).coords2(:,1) ...
   trackInt(pair(1)).intensity_max(:,sigchannel) trackInt(pair(2)).intensity_max(:,sigchannel)];