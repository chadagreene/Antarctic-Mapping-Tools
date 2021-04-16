function h = patchps(lat,lon,varargin)
% patchps works just like Matlab's patch function, but plots georeferenced
% data in Antarctic polar stereographic coordinates (true latitude 71 S).
%
%  patchps(lat,lon,'b') 
% 
% is equivalent to 
% 
%  [x,y] = ll2ps(lat,lon); 
%  patch(x,y,'b')
% 
%% Syntax
% 
%  patchps(lat,lon,cdata)
%  patchps(...,z,cdata)
%  patchps(...,'PropertyName',PropertyValue,...)
%  patchps(...,'km') 
%  patchps(...,'meridian',meridian)
%  h = patchps(...)
% 
%% Description 
% 
% patchps(lat,lon,cdata) creates a patch object with an outline given by arrays 
% georeferenced coordinates lat and lon. Patch objects are then plotted in South 
% polar stereographic eastings and northings. 
% 
% patchps(...,z,cdata) specifies line or marker style. 
% 
% patchps(...,'PropertyName',PropertyValue,...) specifies any number of
% patch properties. 
% 
% patchps(...,'km') plots in polar stereographic kilometers instead of the
% default meters. 
% 
% patchps(...,'meridian',meridian) specifies a meridian longitude in the 
% polar stereographic coordinate conversion. Default meridian is 0. 
% 
% h = patchps(...) returns a column vector of handles to patch objects. .
% 
%% Example
% 
% Using ASAID grounding line data from the asaid Antarctic Mapping Tools
% plugin, 
% 
% load asaid_gl
% patchps(lat,lon,'b')
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
% This function was written by Chad A. Greene of the University of Texas
% Institute for Geophysics (UTIG), February 2015, for inclusion in the
% Antarctic Mapping Tools package. Updated July 2015 to allow plotting 
% in units of polar stereographic kilometers. 
%
% See also: patch, patchm, plotps, pcolorps, and ll2ps. 

%% Input checks: 

assert(nargin>1,'The patchps function requires at least two input: lat and lon.') 
assert(isnumeric(lat)==1,'patchps requires numeric inputs first.') 
assert(isnumeric(lon)==1,'patchps requires numeric inputs first.') 
assert(max(abs(lat(:)))<=90,'I suspect you have entered silly data into patchps because some of your latitudes have absolute values exceeding 90 degrees.') 

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

%% Get initial conditions 

da = daspect; 
da = [1 1 da(3)]; 
hld = ishold; 
hold on

%% Convert units and plot: 

[x,y] = ll2ps(lat,lon,'meridian',meridian); 

% Convert to kilometers if user requested:
if plotkm
    x = x/1000; 
    y = y/1000; 
end

h = patch(x,y,varargin{:}); 
hold on

%% Put things back the way we found them: 

daspect(da)

if ~hld
   hold off
end

%% Clean up: 

if nargout==0
    clear h
end

end

