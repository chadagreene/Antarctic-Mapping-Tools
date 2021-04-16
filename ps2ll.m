function [lat,lon] = ps2ll(x,y,varargin)
%PS2LL transforms map coordinates to lat/lon data for a polar stereographic 
%system. This is a version of Andy Bliss' polarstereo_inv function, adapted 
%specifically for Antarctica. This function does NOT require the Mapping
%Toolbox. 
% 
%% Syntax
% 
% [lat,lon] = ps2ll(x,y) 
% [lat,lon] = ps2ll(x,y,'TrueLat',ReferenceLatitude) 
% [lat,lon] = ps2ll(x,y,'EarthRadius',RadiusInMeters) 
% [lat,lon] = ps2ll(x,y,'Eccentricity',EarthsMisshapenness) 
% [lat,lon] = ps2ll(x,y,'meridian',MeridianInDegrees) 
% 
%% Description 
% 
% [lat,lon] = ps2ll(x,y) transforms polar stereographic x,y coordinates (re: 
% 71 S) to geographic lat/lon. Inputs x and y  can be scalar, vecotr, or
% matrices of equal size. 
% 
% [lat,lon] = ps2ll(x,y,'TrueLat',ReferenceLatitude) secifies a reference
% latitude of true scale in degrees; also known as the standard parallel.
% Note that although Andy Bliss' polarstereo_inv function used -70 as a
% default, this function uses -71 as the default. NSIDC has been trying to
% standardize this, but take a close look at their reference latitudes for
% any data you put through this function--NSIDC sometimes uses 70 S, and
% sometimes uses 71 S. Again, the default in this function is -71. 
% 
% [lat,lon] = ps2ll(x,y,'EarthRadius',RadiusInMeters) specifies Earth's
% radius in meters. Default is 6378137.0 m, WGS84.
% 
% [lat,lon] = ps2ll(x,y,'Eccentricity',EarthsMisshapenness) specifies
% Earth's eccentricity or misshappenness.  Default values is 0.08181919. 
% 
% [lat,lon] = ps2ll(x,y,'meridian',MeridianInDegrees) specifies the meridian in 
% degrees along the positive Y axis of the map. Default value is 0.
% 
%% Snyder's example: Should return lat = -75 and lon = 150.
% 
% x = -1540033.6;
% y = -560526.4;
% [lat,lon] = ps2ll(x,y,'EarthRadius',6378388.0,'eccentricity',0.0819919,'meridian',-100)
% 
%% Author Info
% 
% This function is a slightly adapted version of Andy Bliss' polarstereo_inv, 
% which can be found here: http://www.mathworks.com/matlabcentral/fileexchange/32907
% You can contact Andy at andybliss@gmail.com. 
% 
% This function was tweaked a bit by Chad A. Greene of the University of Texas 
% at Austin's Institute for Geophysics (UTIG). Changes Chad made include removal
% of deg2rad and rad2deg to remove dependence on Mapping Toolbox, and a change to 
% 71 degrees South as the reference latitude. 
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
%% Futher Reading
%   
% Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987 
% http://kartoweb.itc.nl/geometrics/Publications/Map%20Projections%20-%20A%20Working%20manual%20-%20by%20J.P.%20Snyder.pdf
% See the section on Polar Stereographic, with a south polar aspect and
% known phi_c not at the pole.
%
% WGS84 - radius: 6378137.0 eccentricity: 0.08181919
%   in Matlab: axes2ecc(6378137.0, 6356752.3142)
% Hughes ellipsoid - radius: 6378.273 km eccentricity: 0.081816153
%   Used for SSM/I  http://nsidc.org/data/polar_stereo/ps_grids.html
% International ellipsoid (following Snyder) - radius: 6378388.0 eccentricity: 0.0819919 
% 
% See also: LL2PS, PROJINV, PROJFWD, MINVTRAN, MFWDTRAN, and ROTATEM.

%% Input checks: 

assert(nargin>1,'The ps2ll function requires at least two inputs: mapx and mapy.')

%% Set defaults: 

phi_c = -71;   % standard parallel - this is different from Andy Bliss' function, which uses -70! 
a = 6378137.0; % radius of ellipsoid, WGS84
e = 0.08181919;% eccentricity, WGS84
lambda_0 = 0;  % meridian along positive Y axis

%% Parse user inputs: 

tmp = strncmpi(varargin,'true',4)|strncmpi(varargin,'lat',3)|strncmpi(varargin,'ref',3)|...
    strncmpi(varargin,'earthrad',8)|strncmpi(varargin,'rad',3)|...
    strncmpi(varargin,'ecc',3)|strncmpi(varargin,'PosYLon',4)|...
    strncmpi(varargin,'lon',3)|strncmpi(varargin,'merid',5); 
assert(sum(tmp)==(nargin-2)/2,'There seems to be at least one invalid input string. Are you trying to declare options that do not exist?')

tmp = strncmpi(varargin,'true',4)|strncmpi(varargin,'lat',3)|strncmpi(varargin,'ref',3); 
if any(tmp)
    phi_c = varargin{find(tmp)+1}; 
    assert(isscalar(phi_c)==1,'True lat must be a scalar.')
    if phi_c>0
        disp('I''m assuming you forgot the negative sign for the true latitude, and I am converting your northern hemisphere value to southern hemisphere.')
        phi_c = -phi_c; 
    end
end

tmp = strncmpi(varargin,'earthrad',8)|strncmpi(varargin,'rad',3); 
if any(tmp)
    a = varargin{find(tmp)+1}; 
    assert(isscalar(a)==1,'Earth radius must be a scalar.')
    assert(a > 7e+3,'Earth radius should be something like 6378137 in meters.')
end

tmp = strncmpi(varargin,'ecc',3); 
if any(tmp)
    e = varargin{find(tmp)+1}; 
    assert(isscalar(e)==1,'Earth eccentricity must be a scalar.')
    assert(e>0 & e<1,'Earth eccentricity does not seem like a logical value.')
end

tmp = strncmpi(varargin,'PosYLon',4)|strncmpi(varargin,'lon',3)|...
    strncmpi(varargin,'merid',5); 
if any(tmp)
    lambda_0 = varargin{find(tmp)+1}; 
    assert(isscalar(lambda_0)==1,'PosYLon must be a scalar.')
    assert(lambda_0>=-180 & lambda_0<=360,'PsosYLon does not seem like a logical value.')
end

%% Convert to radians and switch sign because this function is southern-hemisphere-specific: 

phi_c = -phi_c*pi/180;
lambda_0 = -lambda_0*pi/180;
x=-x;
y=-y;

%this is not commented very well. See Snyder for details.
t_c=tan(pi/4 - phi_c/2)/((1-e*sin(phi_c))/(1+e*sin(phi_c)))^(e/2);
m_c=cos(phi_c)/sqrt(1-e^2*(sin(phi_c))^2);
rho=sqrt(x.^2+y.^2); 
t=rho*t_c/(a*m_c);

%find phi with a series instead of iterating.
chi=pi/2 - 2 * atan(t);
lat=chi+(e^2/2 + 5*e^4/24 + e^6/12 + 13*e^8/360)*sin(2*chi)...
    + (7*e^4/48 + 29*e^6/240 + 811*e^8/11520)*sin(4*chi)...
    + (7*e^6/120+81*e^8/1120)*sin(6*chi)...
    + (4279*e^8/161280)*sin(8*chi);

lon=lambda_0 + atan2(x,-y);

%correct the signs and phasing
lat=-lat;
lon=-lon;
lon=mod(lon+pi,2*pi)-pi; %want longitude in the range -pi to pi

%convert back to degrees
lat=lat*180/pi;
lon=lon*180/pi;

%% Make two-column format if user requested no outputs: 

if nargout==0
    lat = [lat(:) lon(:)]; 
end

