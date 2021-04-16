function h = graticuleps(varargin)
% graticuleps places a graticule on a polar stereographic cartesian
% coordinate map of the southern hemisphere. 
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
%% Syntax 
% 
%  graticuleps
%  graticuleps(lats,lons) 
%  graticuleps(...,'km') 
%  graticuleps(...,LineProperty,LineValue) 
%  graticuleps(...,'clipping',true) 
%  h = graticuleps(...)
% 
%% Description  
% 
% graticuleps places a grey graticule on a cartesian-coordinate map. 
% 
% graticuleps(lats,lons) specifies lines of latitude and longitude.  By
% default, lats = -80:10:-30 and lons = -150:30:180. 
%  
% graticuleps(...,'km') plots in polar stereographic kilometers instead of
% the default meters.  
% 
% graticuleps(...,LineProperty,LineValue) formats line properties such as 
% color, linewidth, linestyle, etc. 
% 
% graticuleps(...,'clipping',true) deletes all graticule points outside the extents
% of the current plot. 
% 
% h = graticuleps(...) returns a graphics handle h of the plotted graticule. 
% 
%% Examples: 
% 
% graticuleps 
% graticuleps(-85:5:-60,-135:45:180,'r:','linewidth',2,'km')
%
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Intitute for Geophysics, July 2015. It is included in Antarctic Mapping
% Tools for Matlab. Feel free to drop me a line via http://www.chadagreene.com. 
% Dec 3, 2015: Fixed an issue where axis limits were being changed when graticuleps was called.    
% 
% See also: plotps and antmap. 

%% Set defaults: 

lats = [-85 -80:10:-30]; 
lons = -150:30:180; 
clipping = false; % do not trim data outside plot box 

%% Parse inputs: 

% Has user specified lats and lons? 
if nargin>1
    if isnumeric(varargin{1}) % assume first two inputs are lats and lons.   
        lats = varargin{1}; 
        lons = varargin{2}; 
        tmp = true(size(varargin)); 
        tmp(1:2) = false; 
        varargin = varargin(tmp); 
    end
    
    tmp = strncmpi(varargin,'clipping',4); 
    if any(tmp)
       clipping = varargin{find(tmp)+1}; 
       assert(islogical(clipping)==1,'Error: clipping argument must be logical.') 
       tmp(find(tmp)+1) = 1; 
       varargin = varargin(~tmp); 
    end
    
end

%% Build arrays: 

np = 360; % 360 points per line

% Create lat/lon points of all lines:  
lat = NaN(length(lats)*np + length(lats) + length(lons)*np + length(lons)-1,1); 
lon = lat; 

n = 0; 
% Populate arrays of points for all lines of latitude:
for k = 1:length(lats) 
    lat(n+1:n+np) = lats(k); 
    lon(n+1:n+np) = linspace(-180,180,np); 
    n = n+np+1; % each line segment is separated by a NaN. 
end

% Populate array of points for all lines of longitude: 
for k = 1:length(lons) 
    lon(n+1:n+np) = lons(k); 
    lat(n+1:n+np) = linspace(min(lats),max(lats),np); 
    n = n+np+1; 
end

%% Clip graticule to current axis limits: 

ax = axis; % current axis limits 

if isequal(ax,[0 1 0 1])
    mapwasopen = false; 
else
    mapwasopen = true; 
end

if mapwasopen
    
    % If user requests polar stereographic kilometers, we have to convert current axis limits to km for comparison: 
    if any(strcmpi(varargin,'km')) 
        ax = ax*1000; 
    end
    
    if clipping
       % Indices of x,y inside current axis limits: 
       ind = inpsquad(lat,lon,ax(1:2),ax(3:4),'inclusive');

       % Set everything outside current axis limits to NaN:
       lat(~ind) = NaN;
       lon(~ind) = NaN; 
    end
end

%% Get initial conditions 

da = daspect; 
da = [1 1 da(3)]; 
hld = ishold; 
hold on

%% Plot arrays:

h = plotps(lat,lon,'color',0.5*[1 1 1],varargin{:}); 

% Set or reset axis limits: 
if mapwasopen
    
   if any(strcmpi(varargin,'km')) 
      ax = ax/1000; 
   end
   axis(ax)
end
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