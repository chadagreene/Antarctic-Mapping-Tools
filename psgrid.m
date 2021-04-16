function [out1,out2] = psgrid(varargin)
% psgrid creates a grid of specified spatial resolution in polar
% stereographic coordinates. 
% 
%% Syntax
% 
%  [lat,lon] = psgrid(CenterLat,CenterLon,width_km,resolution_km)
%  [lat,lon] = psgrid(CenterX,CenterY,width_km,resolution_km)
%  [lat,lon] = psgrid(CenterName,width_km,resolution_km)
%  [x,y] = psgrid(...,'xy')
% 
%% Description 
% 
% [lat,lon] = psgrid(CenterLat,CenterLon,width_km,resolution_km) returns
% a grid of width_km width in kilometers, resolution_km resolution in
% kilometers, centered on the location given by CenterLat, CenterLon.
% If width_km is a two-element vector, the first element is interpreted as 
% width and the second element is height. If resolution_km is a two element vector, 
% the first element is interpreted as horizontal resolution and the second
% element is interpreted as vertical resolution. 
%
% [lat,lon] = psgrid(CenterX,CenterY,width_km,resolution_km) centers the
% grid on polar stereographic coordinate (CenterX,CenterY) in meters. Polar 
% stereographic meters automatically detected if absolute values exceed normal 
% values of geo coordinates. 
% 
% [lat,lon] = psgrid(CenterName,width_km,resolution_km) centers a grid on
% any location or feature found in the SCAR database. 
% 
% [x,y] = psgrid(...,'xy') returns a grid in polar stereographic meters (ps71). 
% 
%% Example 1  
% Get an 600 km wide, 5 km resolution grid of geo coordinates centered at
% (75°S,107°W) and plot on a Landsat Image Mosaic of Antarctica map as blue dots:  
% 
% [lat,lon] = psgrid(-75,-107,600,5);
% lima(-75,-107,800)
% plotm(lat,lon,'b.') 
% 
%% Example 2 
% In a standard polar stereographic projection, Amery Ice Shelf is wider
% than it is tall.  Let's create a 10 km grid over Amery and make it
% 900 km wide by 250 km tall.  Also, we'll include the 'xy' option to
% return polar stereographic cartesian coordinates rather than the default
% geo coordinates. Then plot in cartesian coordinates:  
% 
% [X,Y] = psgrid('amery ice shelf',[900 250],10,'xy'); 
% 
% bedmap2('patchshelves','xy') 
% bedmap2('patchgl','xy') 
% plot(X,Y,'b.')
% axis(1000*[757 3093 -170 1673])
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
% Austin's Institute for Geophysics (UTIG), September 2015. 
% http://www.chadagreene.com 
% 
% See also ps2ll, ll2ps, and meshgrid. 


%% Error checks: 

narginchk(3,5) 
assert(nargout==2,'psgrid requires exactly two outputs--either lat and lon or polar stereographic x and y.') 

%% Parse Inputs: 

% If first input is numeric, assume it's lat or x:  
if isnumeric(varargin{1})
    assert(isnumeric(varargin{2})==1,'If the first input to psgrid is not a location name, the first two inputs must be lat and lon or polar sterographic x and y.')
    assert(isscalar(varargin{1})==1,'Input coordinates must be scalar.') 
    assert(isscalar(varargin{2})==1,'Input coordinates must be scalar.') 
    
    if islatlon(varargin{1},varargin{2})
        [centerx,centery] = ll2ps(varargin{1},varargin{2}); 
    else
        centerx = varargin{1}; 
        centery = varargin{2}; 
    end
    
    % If coordinates are defined, grid width and resolution are inputs 3 and 4: 
    width_km = varargin{3}; 
    resolution_km = varargin{4}; 
   
else % If first input is not numeric, assume it's a place name: 
    [centerx,centery] = scarloc(varargin{1},'xy'); 
    
    % If location name is defined, grid width and resolution are inputs 2 and 3: 
    width_km = varargin{2}; 
    resolution_km = varargin{3};  
end

% Define x and y values for grid width: 
switch numel(width_km)
    case 1
        widthx = width_km*1000; % The *1000 bit converts from km to meters. 
        widthy = width_km*1000; 

    case 2
        widthx = width_km(1)*1000; 
        widthy = width_km(2)*1000; 
        
    otherwise
        error('I must have misinterpreted something. As I understand it, you have requested a grid width with more than two elements. Check inputs and try again.') 
end
        
% Define x and y values for resolution: 
switch numel(resolution_km)
    case 1
        resx = resolution_km*1000; % The *1000 bit converts from km to meters. 
        resy = resolution_km*1000; 

    case 2
        resx = resolution_km(1)*1000; 
        resy = resolution_km(2)*1000; 
        
    otherwise
        error('I must have misinterpreted something. As I understand it, you have requested a grid resolution with more than two elements. Check inputs and try again.') 
end

% Verify that resolution is not greater than width: 
assert(widthx>resx,'It looks like there''s an input error in psgrid because the grid width should be bigger than the grid resolution. Check inputs and try again.') 
assert(widthy>resy,'It looks like there''s an input error in psgrid because the grid width should be bigger than the grid resolution. Check inputs and try again.') 
assert(resx>0,'Grid resolution must be greater than zero.') 
assert(resy>0,'Grid resolution must be greater than zero.') 
assert(widthx>0,'Grid width must be greater than zero.') 
assert(widthy>0,'Grid width must be greater than zero.') 

% Should outputs be polar stereographic? 
if any(strcmpi(varargin,'xy'))
    outputps = true; 
else 
    outputps = false;
end

%% Build grid: 

x = centerx-widthx/2 : resx : centerx+widthx/2; 
y = centery-widthy/2 : resy : centery+widthy/2; 

[X,Y] = meshgrid(x,y); 

%% Convert coordinates if necessary: 

if outputps
    out1 = X; 
    out2 = Y; 
else
    [out1,out2] = ps2ll(X,Y); 
end

end

