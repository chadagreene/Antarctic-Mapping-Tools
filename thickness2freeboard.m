function F = thickness2freeboard(T,varargin) 
% thickness2freeboard estimates freeboard height above sea level, from ice thickneess
% assuming hyrostatic equilibrium. 
% 
%% Syntax 
% 
%  F = thickness2freeboard(T) 
%  F = thickness2freeboard(...,'rhoi',iceDensity) 
%  F = thickness2freeboard(...,'rhow',waterDensity) 
%  F = thickness2freeboard(...,'rhos',snowDensity) 
%  F = thickness2freeboard(...,'Ts',snowThickness) 
%  
%% Description 
% 
% F = thickness2freeboard(T) estimatesfreeboard height height above sea level F in
% meters above the geoid from ice thickness T in meters. 
% 
% F = thickness2freeboard(...,'rhoi',iceDensity) speccifies ice density in kg/m^3. 
% Default ice density is 917 kg/m^3. 
% 
% F = thickness2freeboard(...,'rhow',waterDensity) specifies water density in kg/m^3. 
% Default water density is 1027 kg/m^3. 
% 
% F = thickness2freeboard(...,'rhos',snowDensity) specifies snow density in kg/m^3.
% Default snow density is 350 kg/m^3, however the default snow thickness is 0 m, so
% the snow density value will only affect calculations if snow thickness is specified.
% 
% F = thickness2freeboard(...,'Ts',snowThickness) specifies snow thickness in meters. 
% Default snow thickness is 0 m. 
%  
%% Example 1: 
% An iceberg is 5 m thick and it's made of pure ice. How high is the iceberg's surface 
% above sea level?
% 
%    F = thickness2freeboard(5)
%    F = 
%         0.54
%
% 
%% Example 2: 
% An iceberg is 4.83 m thick, including 40 cm of snow which has a density of 300 kg/m^3.  
% How high is the iceberg's surface above sea level?
% 
%    F = thickness2freeboard(4.83,'Ts',0.4,'rhos',300)
%    F = 
%         0.80
% 
% And you can go back the other way too: 
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

F = (T + Ts.*(rhow-rhos)./(rhow-rhoi))./(rhow./(rhow-rhoi)); 

end