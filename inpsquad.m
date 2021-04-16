function varargout = inpsquad(lat_or_x,lon_or_y,latlim_or_xlim,lonlim_or_ylim,varargin) 
% inpsquad returns true for points in a polar stereographic quadrangle. 
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
%  tf = inpsquad(lat,lon,latlim,lonlim) 
%  tf = inpsquad(lat,lon,xlim,ylim)
%  tf = inpsquad(x,y,latlim,lonlim)
%  tf = inpsquad(x,y,xlim,ylim)
%  tf = inpsquad(...,'inclusive')
%  [rows,cols] = inpsquad(...)
% 
%% Description 
% 
% tf = inpsquad(lat,lon,latlim,lonlim) returns logical matrix tf, which is
% the same size as lat and lon.  tf is true for all points inside the
% polar stereographic extents of all points in latlim, lonlim.  lat and lon
% must be the same size.  latlim and lonlim must be the same size.  With
% this syntax, all points are converted to polar stereographic (re 71°) meters 
% and limits are taken as the limits of the projected values before performing 
% inpolygon calculation. 
% 
% tf = inpsquad(lat,lon,xlim,ylim) as above, but data limits are defined
% by polar stereographic meters. Polar stereographic units are determined based 
% on the ranges of inputs with the islatlon function. 
%
% tf = inpsquad(x,y,latlim,lonlim) Input data points can be polar
% stereographic meters. Polar stereographic units are determined based on the
% ranges of inputs with the islatlon function. 
%
% tf = inpsquad(x,y,xlim,ylim) All inputs can be polar stereographic
% meters. Polar stereographic units are determined based on the
% ranges of inputs with the islatlon function. 
%
% tf = inpsquad(...,'inclusive') includes points on the edge of the
% polygon defined by xlim,ylim or latlim,lonlim. 
% 
% [rows,cols] = inpsquad(...) returns rows and columns of lat,lon or x,y that have 
% *any* points inside the polar stereographic quadrangle.  Note, a row or column
% needs only one point inside the quadrangle to return true for that row or
% column. 
%
%% Example
% Define a grid and some scattered data: 
% 
% [latgrid,longrid] = psgrid('pine island glacier',500,2); 
% scatlat = -75 + randn(15,1); 
% scatlon = -100 + 2*randn(15,1); 
% mapzoom('pine island glacier',1500) 
% plotm(latgrid,longrid,'b.','markersize',2)
% bedmap2 gl
% plotm(scatlat,scatlon,'rh')
%
% Find all scattered points inside the polar stereographic range of the gridded points: 
% 
% in = inpsquad(scatlat,scatlon,latgrid,longrid)
% plotm(scatlat(in),scatlon(in),'mo')
% 
% Or find find all gridded points inside the range of scattered points: 
% 
% in = inpsquad(latgrid,longrid,scatlat,scatlon); 
% plotm(latgrid(in),longrid(in),'k.')
% 
% This function can be used to trim large datasets to a region of interest: 
% [r,c] = inpsquad(latgrid,longrid,scatlat,scatlon);
% trimlat = latgrid(r,c); 
% trimlon = longrid(r,c); 
% 
%% Author Info
% This function was written for Antarctic Mapping Tools for Matlab
% by Chad A. Greene of the University of Texas at Austin's Institute 
% for Geophysics (UTIG), September 2015.  http://www.chadagreene.com. 
% 
% See also: psgrid, ingeoquad, inpolygon, and find. 

%% Error checks: 

narginchk(4,5)
assert(isequal(size(lat_or_x),size(lon_or_y))==1,'Inputs lat_or_x and lon_or_y must be the same size.') 
assert(isequal(size(latlim_or_xlim),size(lonlim_or_ylim))==1,'Inputs latlim_or_xlim and lonlim_or_ylim must be the same size.') 
assert(numel(latlim_or_xlim)>1,'latlim or xlim must have more than one point.')

%% Set defaults: 

inclusive_inpolygon = false; % Return exclusive inpolygon calculation by default. 

%% Parse Inputs: 
% And convert to ps units if necessary. 

% Data points: 
if islatlon(lat_or_x,lon_or_y)
    [x,y] = ll2ps(lat_or_x,lon_or_y); 
else
    x = lat_or_x; 
    y = lon_or_y; 
end

% Requested range of data: 
if islatlon(latlim_or_xlim,lonlim_or_ylim)
    [xi,yi] = ll2ps(latlim_or_xlim,lonlim_or_ylim); 
else
    xi = latlim_or_xlim; 
    yi = lonlim_or_ylim; 
end

% Has user requested inclusive inpolygon? 
if any(strncmpi(varargin,'inclusive',2))
    inclusive_inpolygon = true; 
end

%% Use only edges of limits: 

% Reshape xi and yi: 
xi = xi(:); 
yi = yi(:); 

% Get limits of xi and yi: 
minx = min(xi); 
maxx = max(xi); 
miny = min(yi); 
maxy = max(yi); 

% Build polygon xv, yv: 
xv = [minx maxx maxx minx]; 
yv = [miny miny maxy maxy]; 

%% Find data points inside xv,yv: 

[IN,ON] = inpolygon(x,y,xv,yv); 

% Include points on the edge of the polygon if user requested: 
if inclusive_inpolygon
    IN = any(cat(3,IN,ON),3); 
end

%% Define outputs 

switch nargout
    case {0,1} 
        varargout{1} = IN; 
        
    case 2
        varargout{1} = find(any(IN,2)); % rows
        varargout{2} = find(any(IN,1)); % columns
        
    otherwise
        error('Too many outputs.') 
end

end