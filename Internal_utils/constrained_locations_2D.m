
function [X,Y] = constrained_locations_2D(numpoints, boundsX, boundsY, varargin)

 
% Input parameters: 
%       numpoints: Total # of points 
%       boundsX: [min max] bounds along X axis
%       boundsY: [min max] bounds along Y axis
%
%       Optional inputs : 
%           thresh_dist: the minimum distance as a function of the mean nnd
%           
%           
% 
% Returned values: 
%       X: X positions of generated points 
%       Y: Y positions of generated points 




p = inputParser;
addParameter(p,'thresh_dist',[], @isnumeric);
parse(p,varargin{:});
params = p.Results;


% threshold dist (10% of the mean nnd estimate)  
if isempty(params.thresh_dist)
    params.thresh_dist = 0.1*mean([boundsX boundsY])/ceil(sqrt(numpoints)); 
end
 
% Get range 
dX = diff(boundsX); 
dY = diff(boundsY); 

% Initialize first point
keeperX = rand.*dX + boundsX(1);
keeperY = rand.*dY + boundsY(1);

% Try dropping down more points 
counter = 2;
for k=2:numpoints 
    minDistance = params.thresh_dist - 1; 
    while minDistance<params.thresh_dist    
        % Get a trial point 
        thisX = rand.*dX + boundsX(1);
        thisY = rand.*dY + boundsY(1); 

      % See how far is is away from existing keeper points
      distances = sqrt((thisX-keeperX).^2 + (thisY - keeperY).^2);
      minDistance = min(distances);
    end
    keeperX(counter) = thisX;
    keeperY(counter) = thisY;
    counter = counter + 1;
end


X = keeperX'; 
Y = keeperY'; 

