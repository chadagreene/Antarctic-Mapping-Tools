function h = geoquadps(latlim,lonlim,varargin)
% geoquadps plots a geographic quadrangle in polar stereographic units.
% 
%% Syntax
% 
%  geoquadps(latlim,lonlim)
%  geoquadps(...,LineProperty,LineValue)
%  geoquadps(...,'meridian',meridian)
%  h = geoquadps(...)
% 
%% Description 
% 
% geoquadps(latlim,lonlim) plots a geographic quadrangle bound by the 
% limits of the two-element vectors latlim,lonlim. The latlim variable
% must be [SouthernLimit NorthernLimit] and lonlim must be [WesternLimit EasternLimit].
% 
% geoquadps(...,LineProperty,LineValue) specifies any line properties. 
% 
% geoquadps(...,'meridian',meridian) specifies a meridian longitude along
% which the polar stereographic projection is centered. Default meridian
% is 0, which puts Fimbul Ice Shelf at the top of the map. To center the 
% map on your quandrangle, try mean(lonlim) as the meridian value. 
% 
% h = geoquadps(...) returns a handle h of the plotted object. 
% 
%% Examples 
% Type 
% 
%   showdemo geoquadps_documentation
% 
% for examples. 
% 
%% Author Info
% This function was written by Chad A. Greene of the University 
% of Texas at Austin, August 2018. 
% 
% See also: plotps, mapzoomps, inset_unproj, and inpsquad. 

%% Error checks: 

narginchk(2,Inf) 
assert(isequal(numel(latlim),numel(lonlim),2)==true,'Error: latlim and lonlim must each be two-element arrays.') 
assert(islatlon(latlim,lonlim)==1,'Error: latlim and lonlim must be geographic coordinates.') 

%% Input parsing: 

% The meridian is the vertical line at the top center of the map: 
tmp = strcmpi(varargin,'meridian'); 
if any(tmp) 
   meridian = varargin{find(tmp)+1}; 
   tmp(find(tmp)+1) = true; 
   varargin = varargin(~tmp); 
else
   meridian = 0; 
end

% Plot in meters by default, unless user wants kilometers: 
tmp = strcmpi(varargin,'km'); 
if any(tmp)
   plotkm = true; 
   varargin = varargin(~tmp); 
else 
   plotkm = false; 
end

%% Build a "quadrangle": 

% Unwrap phase if necessary: 
if diff(lonlim)<0
   lonlim(2) = lonlim(2)+360; 
end

% Estimate the curve along constant latitudes with 300 points per segment, 
% but of course the line along constant longitudes will be straight, so no 
% need to approximate with many points. 

lat = [linspace(latlim(1),latlim(1),200),linspace(latlim(2),latlim(2),200),latlim(1)]; 
lon = [linspace(lonlim(1),lonlim(2),200),linspace(lonlim(2),lonlim(1),200),lonlim(1)];
[x,y] = ll2ps(lat,lon,'meridian',meridian); 

if plotkm
   x = x/1000;
   y = y/1000; 
end

%% Get initial conditions of the plot:

da = daspect; 
da = [1 1 da(3)]; 
hld = ishold; 
hold on

%% Plot the quadrangle: 

h = plot(x,y,varargin{:}); 

%% Put things back the way we found them and clean up: 

daspect(da)

if ~hld
   hold off
end

if nargout==0
    clear h
end

end