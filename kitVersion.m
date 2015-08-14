function version=kitVersion(mode)
% KITVERSION Return string containing version number

if nargin < 1
  mode = 1;
end

switch mode
  case 1
    % KiT version.
    version = '1.4.0';
  case 2
    % Jobset structure version.
    version = 4;
  otherwise
    error('Unknown version mode');
end
