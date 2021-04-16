function T = freeboard2thickness(F,varargin) 
% freeboard2thickness estimates ice thickness from height above sea level, assuming
% hyrostatic equilibrium. 
% 
%% Syntax 
% 
%  T = freeboard2thickness(F) 
%  T = freeboard2thickness(...,'rhoi',iceDensity) 
%  T = freeboard2thickness(...,'rhow',waterDensity) 
%  T = freeboard2thickness(...,'rhos',snowDensity) 
%  T = freeboard2thickness(...,'Ts',snowThickness) 
%  
%% Description 
% 
% T = freeboard2thickness(F) estimates the ice thickness T in meters from height
% above the geoid F in meters. 
% 
% T = freeboard2thickness(...,'rhoi',iceDensity) speccifies ice density in kg/m^3. 
% Default ice density is 917 kg/m^3. 
% 
% T = freeboard2thickness(...,'rhow',waterDensity) specifies water density in kg/m^3. 
% Default water density is 1027 kg/m^3. 
% 
% T = freeboard2thickness(...,'rhos',snowDensity) specifies snow density in kg/m^3.
% Default snow density is 350 kg/m^3, however the default snow thickness is 0 m, so
% the snow density value will only affect calculations if snow thickness is specified.
% 
% T = freeboard2thickness(...,'Ts',snowThickness) specifies snow thickness in meters. 
% Default snow thickness is 0 m. 
%  
%% Example 1: 
% The surface of an iceberg is 80 cm above sea level and it's pure ice. How thick
% is the iceberg? 
% 
%    T = freeboard2thickness(0.8)
%    T =
%           7.47
% 
% And you can go back the other way with the sister function: 
% 
%    F = thickness2freeboard(7.47)
%    F =
%           0.8
%
%% Example 2: 
% The surface of an iceberg is 80 cm above sea level, including 40 cm of snow which 
% has a density of 300 kg/m^3.  How thick is the iceberg? 
% 
%    T = freeboard2thickness(0.8,'Ts',0.4,'rhos',300)
%    T =
%           4.83
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
% This function was written by Chad A. Greene of the University of Texas at Austin
% Institute for Geophysics (UTIG), March 2017. 
% 
% See also thickness2freeboard. 


%% Set defaults

rhoi = 917;  % ice density  (can vary from ~915 to 922 kg/m3 or a wider range)
rhow = 1027; % water density (can also vary from 1024 to 1028 or so) 
rhos = 350;  % snow density (can vary a lot!) 
Ts = 0; % assume no snow

%% Parse inputs 

assert(nargin>0,'Not enough inputs.  In fact, no inputs at all...') 

if nargin>1
   tmp = strcmpi(varargin,'rhoi'); 
   if any(tmp)
      rhoi = varargin{find(tmp)+1}; 
   end
   
   tmp = strcmpi(varargin,'rhow'); 
   if any(tmp)
      rhow = varargin{find(tmp)+1}; 
   end
   
   tmp = strcmpi(varargin,'rhos'); 
   if any(tmp)
      rhos = varargin{find(tmp)+1}; 
   end
   
   tmp = strcmpi(varargin,'Ts'); 
   if any(tmp)
      Ts = varargin{find(tmp)+1}; 
   end
end

%% Perform mathematics: 


T = F.*rhow./(rhow-rhoi) - Ts.*(rhow-rhos)./(rhow-rhoi); 

end