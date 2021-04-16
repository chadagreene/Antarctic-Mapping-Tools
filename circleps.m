function h = circleps(lat_or_x,lon_or_y,radius_km,varargin)
% circleps plots circles of given radii on an Antarctic polar stereographic 
% projection map.  Latitude of true scale is 71S. 
% 
%% Syntax 
%
%  circleps(lat,lon,radius_km)
%  circleps(x,y,radius_km)
%  circleps(...,'PropertyName',PropertyValue)
%  circleps(...,'km')
%  circleps(...,'meridian',meridian)
%  h = circleps(...)
% 
%% Description
% 
% circleps(lat,lon,radius_km) plots circle(s) of specified radius in kilometers centered at points 
% given by geo coordinates lat and lon.  
% 
% circleps(x,y,radius_km) lets you input coordinates as polar stereographic meters. Coordinates are 
% automatically determined by the islatlon function. 
% 
% circleps(...,'PropertyName',PropertyValue) specifies patch or fill properties 
% such as 'facecolor' or 'linewidth'. 
%
% circleps(...,'km') plots in polar stereographic kilometers rather than meters. 
% 
% circleps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
%
% h = circleps(...) returns the handle(s) h of the plotted circle(s). 
% 
%% Examples: 
% Given these three places:
% 
%   places = {'mcmurdo station','palmer station','casey station'};
%   [lat,lon] = scarloc(places); 
% 
% Plot circles of 1000 km radius centered on each station:
% 
%   circleps(lat,lon,1000)
% 
% Make them blue filled circles with thick red outlines--and make them 
% have a 500 km radius for McMurdo, 400 km radius for Palmer, and 300 km for Casey: 
% 
%   circleps(lat,lon,[500;400;300],'facecolor','b','edgecolor','r','linewidth',4)
%  
%% Citing Antarctic Mapping Tools
% This function was developed for Antarctic Mapping Tools for Matlab (AMT). If AMT is useful for you,
% please cite our paper: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. 
% http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
% @article{amt,
%   title={{Antarctic Mapping Tools for \textsc{Matlab}}},
%   author={Greene, Chad A and Gwyther, David E and Blankenship, Donald D},
%   journal={Computers \& Geosciences},
%   year={2017},
%   volume={104},
%   pages={151--157},
%   publisher={Elsevier}, 
%   doi={10.1016/j.cageo.2016.08.003}, 
%   url={http://www.sciencedirect.com/science/article/pii/S0098300416302163}
% }
%   
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG). October 2016.  
% http://www.chadagreene.com
% 
% See also plotps, patchps, and pathdistps. 

%% Check inputs: 

narginchk(3,inf) 
assert(isequal(size(lat_or_x),size(lon_or_y))==1,'Input error: Coordinate dimensions do not match.') 
assert(isnumeric(radius_km)==1,'Input error: Make sure radius_km is numeric.') 
if ~isscalar(radius_km) 
   assert(isequal(size(lat_or_x),size(radius_km))==1,'Input error: Circle radius declaration must either be a scalar value or its dimensions must match the dimensions of the circle center coordinates.') 
end

%% Set defaults: 

NOP = 1000; % number of points per circle. 

%% Parse inputs: 


tmp = strcmpi(varargin,'meridian'); 
if any(tmp)
   meridian = varargin{find(tmp)+1}; 
   assert(isscalar(meridian)==1,'Error: meridian must be a scalar longitude.') 
   tmp(find(tmp)+1) = true; 
   varargin = varargin(~tmp); 
else 
   meridian = 0; 
end

% Check input coordinates: 
if islatlon(lat_or_x,lon_or_y)
   [x,y] = ll2ps(lat_or_x,lon_or_y,'meridian',meridian); 
else
   x = lat_or_x; 
   y = lon_or_y; 
end

% Be forgiving if the user enters "color" instead of "facecolor"
tmp = strcmpi(varargin,'color');
if any(tmp)
    varargin{tmp} = 'facecolor'; 
end

% Plot in meters or kilometers? 
tmp = strcmpi(varargin,'km'); 
if any(tmp)
   varargin = varargin(~tmp); 
   x = x/1000; 
   y = y/1000;
   r = radius_km; 
else
   r = radius_km*1000; 
end
   
%% Begin operations:

% Make inputs column vectors: 
x = x(:); 
y = y(:);
r = r(:); 

if isscalar(r)
   r = repmat(r,size(x)); 
end

% Define an independent variable for drawing circle(s):
t = 2*pi/NOP*(1:NOP); 

%% Get initial figure conditions:  

% aspect ratio: 
da = daspect; 
da = [1 1 da(3)]; 

% Query original hold state:
holdState = ishold; 
hold on; 

%% Plot 

% Preallocate object handle: 
h = NaN(size(x)); 

% Plot circles singly: 
for n = 1:numel(x)
    h(n) = fill(x(n)+r(n).*cos(t), y(n)+r(n).*sin(t),'','facecolor','none',varargin{:});
end

%% Clean up: 

if ~holdState
   hold off
end

daspect(da) 

% Delete object handles if not requested by user: 
if nargout==0 
    clear h 
end

end

