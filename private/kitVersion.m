function version=kitVersion(mode)
% KITVERSION Return string containing version number

if nargin < 1
  mode = 1;
end

switch mode
  case 1
    % KiT version.
    version = '1.5.4';
  case 2
    % Jobset structure version.
    version = 6;
  otherwise
    error('Unknown version mode');
end
