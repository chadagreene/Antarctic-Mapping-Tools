function h = surfps(lat,lon,Z,varargin)
% surfps works just like Matlab's surf function, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71 S).
% 
%% Syntax
% 
%  surfps(lat,lon,Z)
%  surfps(...,'PropertyName',PropertyValue,...)
%  surfps(...,'km') 
%  surfps(...,'meridian',meridian)
%  h = surfps(...)
% 
%% Description 
% 
% surfps(lat,lon,Z) plots a surface to represent the data grid Z 
% corresponding to a georeferenced lat,lon grid in South polar stereographic
% eastings and northings. 
% 
% surfps(...,'PropertyName',PropertyValue,...) specifies any number of
% surface properties. 
% 
% surfps(...,'km') plots in polar stereographic kilometers instead of the
% default meters. 
%     
% surfps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
% 
% h = surfps(...) returns a column vector of handles to surface objects.
% 
%% Example
% Load and plot a Bedmap2 DEM of a 100x100 km region around Mt. Vinson: 
% 
% [centerlat,centerlon] = scarloc('mount vinson'); 
% [lat,lon,bed] = bedmap2_data('bed',centerlat,centerlon,50); 
% 
% surfps(lat,lon,bed)
% xlabel 'eastings (m)' 
% ylabel 'northings (m)'
% zlabel 'elevation (m)'
% 
%% Citing Antarctic Mapping Tools
% If AMT is useful for you, please cite the following paper: 
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
% This function was written by Chad Greene of the University of Texas
% Institute for Geophysics (UTIG), November 2015, for inclusion in the
% Antarctic Mapping Tools package. 
%
% See also: pcolorps, plotps, surf, pcolorm, surf, surfm, and ll2ps. 

%% Input checks: 

assert(nargin>2,'The surfps function requires at least three inputs: lat, lon, and Z.') 
assert(islatlon(lat,lon)==1,'I suspect you have entered silly data into surfps because some of the lats or lons fall outside the normal range of geo coordinates.') 

%% Parse inputs

plotkm = false; % by default, plot in meters 
meridian = 0;   % top of the map is Fimbul Ice Shelf

% Has user requested plotting in kilometers? 
if nargin > 2
   tmp = strcmpi(varargin,'km'); 
   if any(tmp)
      plotkm = true; 
      varargin = varargin(~tmp); 
   end
   
   tmp = strcmpi(varargin,'meridian'); 
   if any(tmp)
      meridian = varargin{find(tmp)+1}; 
      assert(isscalar(meridian)==1,'Error: meridian must be a scalar longitude.') 
      tmp(find(tmp)+1) = true; 
      varargin = varargin(~tmp); 
   end
end

%% Convert units and plot: 

[x,y] = ll2ps(lat,lon,'meridian',meridian); 

% Convert to kilometers if user requested:
if plotkm
    x = x/1000; 
    y = y/1000; 
end

h = surf(x,y,Z,varargin{:}); 
shading flat
axis tight
set(gca, 'DataAspectRatio', [diff(get(gca, 'XLim')) diff(get(gca, 'XLim')) diff(get(gca, 'ZLim'))])
hold on; 
grid off

%% Clean up: 

if nargout==0
    clear h
end

end

