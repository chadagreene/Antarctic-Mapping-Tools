function m = psdistortion(lat,TrueLat)
% psdistortion approximates the map scale factor for a polar stereographic projection. 
% It is the ratio of distance on a ps projection to distance on a sphere. 
% 
% The equation is: 
% 
%      1 + sin(TrueLat) 
% m = -------------------
%      1 + sin(lat) 
% 
% and is from Eq 3.3 of Numerical Weather and Climate Prediction 
% by Thomas Tomkins Warner (2011). Cambridge Univ. Press.
% 
%% Syntax
% 
%  m = psdistortion(lat)
%  m = psdistortion(lat,TrueLat)
% 
%% Description 
% 
% m = psdistortion(lat) gives the ratio of distance on a polar stereographic 
% projected map to distance on a sphere. Where m is greater than 1, the polar 
% stereographic projection exaggerates (lengthens) the distance compared to 
% reality. Where m is less than 1, projected distances are shrunk compared 
% to reality. 
% 
% m = psdistortion(lat,TrueLat) specifies a latitude of true scale. Default
% is the standard -71 degrees. 
% 
%% Examples 
% For examples, type 
% 
%   showdemo psdistortion_documentation 
% 
%% Author Info
% This function was written by Chad A. Greene. Thanks to Mathieu 
% Morlighem for suggesting the addition of this function. 
% 
% See also: pathdistps. 

%% Input checks: 

narginchk(1,2) 
assert(max(abs(lat(:)))<=90,'Error: inputs must be latitudes.') 

if nargin>1
   assert(isscalar(TrueLat),'Error: TrueLat must be a scalar.') 
   assert(max(abs(TrueLat(:)))<=90,'Error: TrueLat must be in the range -90 to 90.') 
else
   TrueLat = -71; % the default
end

%% Mathematics: 

m = (1+sind(abs(TrueLat)))./(1+sind(abs(lat)));

end