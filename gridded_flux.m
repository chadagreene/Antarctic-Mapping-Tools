function [flux,flux_err] = gridded_flux(mask,vx,vy,H,res,vx_err,vy_err,H_err)
% gridded_flux calculates the annual ice flux out of a gridded mask. 
% 
%% Syntax
% 
%  flux = gridded_flux(mask,vx,vy,H,res)
%  [flux,flux_err] = gridded_flux(mask,vx,vy,H,res,vx_err,vy_err,H_err)
% 
%% Description 
%
% flux = gridded_flux(mask,vx,vy,H,res) calculates the ice flux in Gt/yr as 
% it flows out of a binary mask that corresponds to velocity grids vx,vy
% (m/yr) and thickness H (m). The resolution of the grids res (m) must be 
% specified as a scalar. Mass flux calculations assume the density of ice
% is 917 kg/m^3. 
% 
% [flux,flux_err] = gridded_flux(mask,vx,vy,H,res,vx_err,vy_err,H_err)
% estimates flux error using the error grids vx_err, vy_err, and H_err. 
% 
%% Examples
% For examples type 
%
%  amt gridded_flux 
% 
%% Author Info: 
% Written by Chad A. Greene of NASA's Jet Propulsion Laboratory, 
% December 2021. 

%% Input checks: 

narginchk(5,8) 
assert(islogical(mask),'Error: Input mask must be logical.') 
assert(isequal(size(mask),size(vx),size(vy),size(H)),'Error: mask, vx, vy, and H must all be 2D grids of equal dimensions.') 
assert(isscalar(res),'Error: Input resolution res must be a scalar.') 
if nargin>5
   %assert(isequal(size(mask),size(vx_err),size(vy_err),size(H_err)),'Error: The dimensions of vx_err, vy_err, and H_err must match the dimensions of vx, vy, and H.') 
   CalculateError = true; 
   assert(nargout==2,'Flux error estimates (two outputs) require vx_err, vy_err, and H_err as inputs.')
else 
   CalculateError = false; 
   assert(nargout==1,'If error estimates are included as inputs, flux_err must be an output, but you have not requested two outputs.') 
end

%% Define constants

rho_ice = 917; % kg/m^3

%% Calculate pixel displacement after an amount of time dt: 

vx(isnan(vx)) = 0; 
vy(isnan(vy)) = 0; 

x = (1:size(mask,2))*res; 
y = (size(mask,1):-1:1)*res;

[X0,Y0] = meshgrid(x,y); 

% Choose a timestep equivalent to one pixel of the fastest ice along the
% perimeter. (Very tiny timesteps are potentially subject to numerical noise
% for very slow moving ice, and all that noise can add up. Very large timesteps
% can potentially jump a fjord and land on the other side or hit an icy island. 
% Therefore, ~1 pixel timestep is a happy medium: 
perim = bwperim(mask,8); 
dt = res/max(hypot(vx(perim),vy(perim)),[],'all'); 

% Pixel center locations after 1 timestep: 
X1 = X0+vx*dt; 
Y1 = Y0+vy*dt; 

% The mass of ice in each grid cell: 
IceMass = H.*rho_ice*(res^2)*1e-12; 

% Initial total mass of ice: 
IceMass0 = double(mask).*IceMass; 

% Total mass of ice at locations after 1 timestep of advection: 
IceMass1 = interp2(x,y,double(mask),X1,Y1).*IceMass; 

flux = sum(IceMass1-IceMass0,'all','omitnan')/dt; 

if CalculateError 
   
   % Here's the flux error due to thickness error, assuming thickness
   % errors are fully correlated: 
   IceMass_err = (H+H_err).*rho_ice*(res^2)*1e-12; 
   flux_H_err = sum((interp2(x,y,double(mask),X1,Y1)-double(mask)).*IceMass_err,'all','omitnan')/dt; 
   
   % Erroneous pixel locations: 
   X1 = X0+(vx + sign(vx).*vx_err)*dt; 
   Y1 = Y0+(vy + sign(vy).*vy_err)*dt; 
   flux_v_err = sum((interp2(x,y,double(mask),X1,Y1)-double(mask)).*IceMass,'all','omitnan')/dt; 
   
   % Assume thickness errors are independent from velocity errors, 
   % so the total flux error estimate is the root sum square of errors due
   % to thickness error and errors due to velocity error: 
   flux_err = rssq([flux_H_err-flux flux_v_err-flux]); 

end