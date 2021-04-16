function F = base2freeboard(B,varargin) 
% base2freeboard estimates freeboard height above sea level, from ice basal
% elevation, assuming hyrostatic equilibrium. 
% 
%% Syntax 
% 
%  F = base2freeboard(B) 
%  F = base2freeboard(...,'rhoi',iceDensity) 
%  F = base2freeboard(...,'rhow',waterDensity) 
%  F = base2freeboard(...,'rhos',snowDensity) 
%  F = base2freeboard(...,'Ts',snowThickness) 
%  
%% Description 
% 
% F = base2freeboard(B) estimates freeboard height height above sea level F in
% meters above the geoid from basal elevation of ice B in meters. 
% 
% F = base2freeboard(...,'rhoi',iceDensity) speccifies ice density in kg/m^3. 
% Default ice density is 917 kg/m^3. 
% 
% F = base2freeboard(...,'rhow',waterDensity) specifies water density in kg/m^3. 
% Default water density is 1027 kg/m^3. 
% 
% F = base2freeboard(...,'rhos',snowDensity) specifies snow density in kg/m^3.
% Default snow density is 350 kg/m^3, however the default snow thickness is 0 m, so
% the snow density value will only affect calculations if snow thickness is specified.
% 
% F = base2freeboard(...,'Ts',snowThickness) specifies snow thickness in meters. 
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
% What if you know the base of the ice is 4.5 m below sea level. Where's
% surface of the ice? 
% 
%    F = base2freeboard(-4.5)
%    F = 
%         0.54
% 
% What if you know there's a 0.75 m thick layer of fluffy snow on top of 
% the iceberg, and the snow is  200 kg/m^3. Where is the surface? 
% 
%    F = base2freeboard(-4.5,'rhos',200,'Ts',0.75)
%    F = 
%         1.22
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
% See also thickness2freeboard and freeboard2thickness. 

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

F = (B - Ts.*((rhow-rhos)./(rhow-rhoi)))./(1-rhow./(rhow-rhoi)); 

% Assume any base elevations above sea level are error or rock: 
F(F<0) = NaN; 

end