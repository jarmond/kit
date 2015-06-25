function lastUpdate = kitProgress(progress, lastUpdate)
% KITPROGRESS Print formatted progress message
%
% SYNOPSIS: kitProgress(progress, lastUpdate)
%
% INPUT progress: Fraction complete [0,1].
%
%       lastUpdate: Last timestamp progress was called with (optional).
%                   Omit for first call.
%
% Copyright (c) 2012 Jonathan W. Armond

currentClock = clock;

% Add timestamp to message.
fmt = [datestr(now, 'HH:MM:SS') ': % 6.2f%% complete\n'];
msgLength = 27;
carriageReturn = repmat('\b', 1, msgLength);

updateInterval = 1.0; % seconds
if nargin>=2 && etime(currentClock, lastUpdate) > updateInterval
    msg = sprintf(fmt, 100*progress);
    % Simulate carriage-return.
    fprintf(carriageReturn);
    % Print message.
    fprintf('%s',deblank(msg));
    lastUpdate = currentClock;
end

if nargin<2
    fprintf(fmt, 100*progress);
    lastUpdate = currentClock;
end
