function [pathDistance] = pathdist(lat,lon,varargin) 
%PATHDIST uses the distance function to calculate cumulative distance
%traveled along a path given by the arrays lat and lon.  (Requires Mapping
%Toolbox). Always assumes WGS84 ellipsoid. 
%
%% Syntax
% 
%  pathDistance = pathdist(lat,lon)
%  pathDistance = pathdist(...,LengthUnit)
%  pathDistance = pathdist(...,track)
%  pathDistance = pathdist(...,'refpoint',[reflat reflon])
% 
%% Description
% 
% pathDistance = pathdist(lat,lon) returns the cumulative distance
% traveled along the path given by (lat,lon). Distance is in meters
% by default, referenced to the WGS84 ellipsoid. The pathDistance array
% will be the same size as lat and lon. 
%
% pathDistance = pathdist(...,LengthUnit) specifies any valid length unit. 
% The following are a few LengthUnit options. See documentation for 
% validateLengthUnit for a complete list of options.
% 
%     Unit Name     LengthUnit String
%     meter	        'm', 'meter(s)', 'metre(s)' (default)
%     kilometer     'km', 'kilometer(s)', 'kilometre(s)'
%     nautical mile	'nm', 'naut mi', 'nautical mile(s)'
%     foot          'ft', 'international ft','foot', 'international foot', 'feet', 'international feet'
%     inch          'in', 'inch', 'inches'
%     yard          'yd', 'yds', 'yard(s)'
%     mile          'mi', 'mile(s)','international mile(s)'
%
% pathDistance = pathdist(...,track) uses the input string track to specify 
% either a great circle/geodesic or a rhumb line arc. If track equals 'gc' (the default 
% value), then great circle distances are computed on a sphere and geodesic distances are 
% computed on the WGS84 ellipsoid. If track equals 'rh', then rhumb line distances are 
% computed on the WGS84 ellipsoid.
%
% pathDistance = pathdist(...,'refpoint',[reflat reflon]) references the
% path distance to the point along the path nearest to [reflat reflon].
% For this calculation, pathdist finds the point in lat and lon
% which is nearest to [reflat reflon] and assumes this point along
% lat,lon is the zero point. This is only an approximation, and may
% give erroneous results in cases of very sharply-curving, crossing, or 
% otherwise spaghetti-like paths; where [reflat reflon] lies far from any
% point along the path, or where points along the path are spaced far
% apart. 
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
%% Author Info and version history
% 
% Written by Chad A. Greene of the University of Texas at Austin. Institute
% for Geophysics. 
% 
% June 23, 2014: Version 1 written and uploaded to the Mathworks File Exchange site. 
% About an hour later: Version 2 written to allow for NaN inputs. 
% 
% January 2015: Improved input parsing and now allows any valid length unit.  
% 
% See also pathdistps, distance, and validateLengthUnit.

%% Check inputs 

assert(license('test','map_toolbox')==1,'The pathdist function requires Matlab''s Mapping Toolbox. Use pathdistps instead!')
assert(numel(lat)==numel(lon),'Length of lat and lon must match.') 
assert(numel(lat)>1,'lat and lon must have more than one point.') 
assert(isvector(lat),'lat and lon must be vectors.') 


%% Parse inputs 

% Set defaults: 
units = 'm'; 
track = 'gc'; 
referenceto = false; 

% Check if reference point is defined:
tmp = strncmpi(varargin,'ref',3); 
if any(tmp) 
    refpt = varargin{find(tmp)+1}; 
    tmp(find(tmp)+1)=1; 
    varargin = varargin(~tmp); 
    assert(numel(refpt)==2,'Coordinates of reference point can be only a single point given by a latitude/longitude pair in the form [reflat reflon].')
    referenceto = true; 
end 

% In a previous version, units and track needed to be called with name-value pairs. Now pathdist 
% is smart enough to figure out tracks and units without explictly saying
% them. So we'll discard unit and track right now: 
tmp = strncmpi(varargin,'unit',4)+strcmpi(varargin,'track')+strcmpi(varargin,'gc')+strncmpi(varargin,'great',5); 
if any(tmp)
    varargin = varargin(~tmp); 
end

% Check for rhumb line track preference: 
tmp = strncmpi(varargin,'rh',5); 
if any(tmp)
    track = 'rh'; 
    varargin = varargin(~tmp); 
end

% Ensure that if any inputs have made it through the above, that all we have is a length unit: 
if ~isempty(varargin) 
    units = varargin{:}; 
    validateLengthUnit(units);
end


%% 

% Preallocate output: 
if isrow(lat)
    pathDistance = NaN(1,length(lat));
else
    pathDistance = NaN(length(lat),1);
end

% Consider only finite entries in case of gps failure or whatnot:  
la = lat(isfinite(lat)); 
lo = lon(isfinite(lon)); 


distancesBetweenPoints = distance(track,la(1:end-1),lo(1:end-1),... % req's Mapping toolbox. 
    la(2:end),lo(2:end),referenceEllipsoid('wgs 84',units));
    
finitelat = isfinite(lat); 
finitelat(find(finitelat,1,'first'))=0; 
pathDistance(finitelat) = cumsum(distancesBetweenPoints); 

pathDistance(1)=0; 

if referenceto
    d = distance(track,lat,lon,refpt(1),refpt(2),referenceEllipsoid('wgs 84',units));
    [~,mind] = min(d); 
    pathDistance = pathDistance - pathDistance(mind);     
end


end

