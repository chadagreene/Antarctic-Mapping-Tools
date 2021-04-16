function [x,y]=ll2ps(lat,lon,varargin)
% ll2ps (lat/lon to polarstereographic) transforms lat/lon coordinates to map 
% coordinates for a polar stereographic system. This is a version of Andy Bliss' 
% polarstereo_fwd function, adapted specifically for Antarctica. This function 
% does NOT require the Mapping Toolbox. 
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
% [x,y] = ll2ps(lat,lon) 
% [x,y] = ll2ps(lat,lon,'TrueLat',ReferenceLatitude) 
% [x,y] = ll2ps(lat,lon,'EarthRadius',RadiusInMeters) 
% [x,y] = ll2ps(lat,lon,'Eccentricity',EarthsMisshapenness) 
% [x,y] = ll2ps(lat,lon,'meridian',MeridianInDegrees) 
% 
%% Description 
% 
% [x,y] = ll2ps(lat,lon) transforms georeferenced coordinates to
% polarstereographic x,y coordinates referenced to 71 S. Inputs lat and lon
% can be scalar, vecotr, or matrices of equal size. 
% 
% [x,y] = ll2ps(lat,lon,'TrueLat',ReferenceLatitude) secifies a reference
% latitude of true scale in degrees; also known as the standard parallel.
% Note that although Andy Bliss' polarstereo_fwd function used -70 as a
% default, this function uses -71 as the default. NSIDC has been trying to
% standardize this, but take a close look at their reference latitudes for
% any data you put through this function--NSIDC sometimes uses 70 S, and
% sometimes uses 71 S. Again, the default in this function is -71. 
% 
% [x,y] = ll2ps(lat,lon,'EarthRadius',RadiusInMeters) specifies Earth's
% radius in meters. Default is 6378137.0 m, WGS84.
% 
% [x,y] = ll2ps(lat,lon,'Eccentricity',EarthsMisshapenness) specifies
% Earth's eccentricity or misshappenness.  Default values is 0.08181919. 
% 
% [x,y] = ll2ps(lat,lon,'meridian',MeridianInDegrees) specifies the meridian in 
% degrees along the positive Y axis of the map. Default value is 0.
% 
%
%% Snyder's example. Should return x = -1540033.6; y = -560526.4;
% 
% lat = -75; 
% lon = 150;
% 
% [x,y]=ll2ps(lat,lon,'EarthRadius',6378388,'Eccentricity',0.0819919,'meridian',-100)
% 
%% Author Info
% 
% This function is a slightly adapted version of Andy Bliss' polarstereo_fwd, 
% which can be found here: http://www.mathworks.com/matlabcentral/fileexchange/32950
% You can contact Andy at andybliss@gmail.com. 
% 
% This function was tweaked a bit by Chad A. Greene of the University of Texas 
% at Austin's Institute for Geophysics (UTIG). Changes Chad made include removal
% of deg2rad and rad2deg to remove dependence on Mapping Toolbox, and a change to 
% 71 degrees South as the reference latitude. 
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
% See also: PS2LL, PROJINV, PROJFWD, MINVTRAN, MFWDTRAN, and ROTATEM.

%% Input checks: 

assert(nargin>1,'The ll2ps function requires at least two inputs: lat and lon.')
assert(sum(any(lat>0))==0,'Input latitudes must be in southern hemisphere, but you seem to have input positive latitudes.')
assert(sum(any(lon>360))==0,'Input longitudes greater than 360 degrees? This seems wrong.')
assert(sum(any(lon<-360))==0,'Input longitudes less than -360 degrees? This seems wrong.')


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

%% Convert to radians:

lat = lat*pi/180;
phi_c = phi_c*pi/180;
lon = lon*pi/180;
lambda_0 = lambda_0*pi/180;

%this is not commented very well. See Snyder for details.
t=tan(pi/4+lat/2)./((1-e*sin(-lat))./(1+e*sin(-lat))).^(e/2); 

t_c=tan(pi/4 + phi_c/2)./((1-e*sin(-phi_c))./(1+e*sin(-phi_c))).^(e/2);
m_c=cos(-phi_c)./sqrt(1-e^2*(sin(-phi_c)).^2);
rho=a*m_c*t/t_c; %true scale at lat phi_c

x=-rho.*sin(-lon+lambda_0);
y=rho.*cos(-lon + lambda_0);


%% Make two-column format if user requested no outputs: 

if nargout==0
    x = [x(:) y(:)]; 
end


