function zs = ice_profile_smoother(d,z,H,varargin) 
% ice_profile_smoother smooths any variable along a glacier flowline, as a 
% function of local ice thickness. The smoothing window is an exponential 
% shape, which is the best approximation of that stress the ice "feels", and 
% unlike the jittery profiles produced by a simple unweighted (boxcar) filter, 
% the exponential window produces smooth profiles. 
% 
% Usage note: This function assumes that postings along d are at least somewhat
% equally spaced. If the gradient(d) varies wildly, the weighting of the moving
% window may not be accurate. Try using pathdistps or pathdistpsn if you want 
% to ensure the spacing along your flowline is equally spaced. 
% 
%% Syntax
% 
%  zs = ice_profile_smoother(d,z,H)
%  zs = ice_profile_smoother(...,'CouplingLength',Nthck)
%  zs = ice_profile_smoother(...,'endpoints','fill')
%  zs = ice_profile_smoother(...,'weights',w)
% 
%% Description 
% 
% zs = ice_profile_smoother(d,z,H) smooths the variable z along a flowline where
% d specifies distance in meters along the flowline and H is the corresponding
% ice thickness. 
%
% zs = ice_profile_smoother(...,'CouplingLength',Nthck) specifies a longitudinal 
% coupling length as a multiple of ice thickesses. This is equivalent to l/H
% in Kamb & Echelmeyer's  paper cited below. Important: Nthick is not the same 
% thing as the total window width. The Kamb & Echelmeyer paper describes it 
% in detail, but the "averaging length" is the full width of a boxcar window 
% and is equal to 4*l. In this function, the default value of Nthck is 2.5, 
% which is equivalent to a moving average window width of 10 ice thicknesses. 
% 
%   For guidance on choosing a value of Nthck, Kamb & Echelmeyer state that
%   "l/H ranges from about 1.5 to 10...for temperate valley glaciers, with f 
%   near 0.5 and with longitudinal strain-rates typically of order 0.01-0.05 /yr, 
%   l/H should be in the range from about 1 to 3, whereas for ice sheets ...
%   the expected l/H is in the range from about 4 to 10, distinctly higher
%   than for valley glaciers."
% 
% zs = ice_profile_smoother(...,'endpoints','fill') similar to the 'endpoints'
% option in the movmean function, the 'fill' option sets the ends of the profile 
% to NaN. This option means you'll lose some data at the end of a profile,
% particularly where ice is thick and when Nthck is large, but it ensures that 
% all resulting data are properly weighted and smoothed. The 'fill' option 
% will also mean you'll lose data anywhere near any NaN values in the input 
% thickness H. 
%
% zs = ice_profile_smoother(...,'weights',w) applies weights to each observation
% within the smoothing window. This is useful if each measurement along a flowline
% has its own error estimate. The variable w must be the same size as d, z, and H, 
% and typically if your measurements z have corresponding 1-sigma error estimates
% z_err, then w=1./z_err.^2. 
% 
%% Example 
% Smooth a surface elevation profile to a coupling length of the default
% 2.5 ice thicknesses:
% 
%  load Kangilerngata_Sermia_flowline 
%  
%  d = pathdistpsn(xfl,yfl); 
%  H = bedmachine_interp('thickness',xfl,yfl,'greenland'); 
%  sfz = bedmachine_interp('surface',xfl,yfl,'greenland'); 
%  
%  sfz_smooth = ice_profile_smoother(d,sfz,H); 
% 
%  plot(d,sfz)
%  hold on
%  plot(d,sfz_smooth,'linewidth',2)
% 
% Try a longer coupling length, and avoid any potential errors at the edges
% of the window by filling the endpoints with NaN (rather than the default, 
% which shrinks the moving window size at the edges of the profile): 
% 
%  sfz_smooth2 = ice_profile_smoother(d,sfz,H,'CouplingLength',5,'endpoints','fill'); 
%  plot(d,sfz_smooth2,'.')
% 
%% Citing this function
% The theory in this function is entirely taken from Kamb and Echelmeyer's 
% 1986 paper. Please cite it. And for repeatability, and to do me a kindness, 
% please also cite my Antarctic Mapping Tools paper: 
% 
% Kamb, B., & Echelmeyer, K. (1986). Stress-Gradient Coupling in Glacier Flow: 
% I. Longitudinal Averaging of the Influence of Ice Thickness and Surface Slope. 
% Journal of Glaciology, 32(111), 267-284. doi:10.3189/S0022143000015604
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. (2017). Antarctic Mapping 
% Tools for MATLAB. Computers & Geosciences, 104, 151-157. 
% https://doi.org/10.1016/j.cageo.2016.08.003

%% Error checks: 

assert(isvector(d),'Error; d must be a vector describing the length in meters along a profile.') 
assert(isequal(size(d),size(z),size(H)),'Dimenions of d, z, and H must all agree.') 
assert(issorted(d(isfinite(d)),'monotonic'),'Spacing of d must be monotonic.') 
if max(d)<1000
   warning('This profile appears to only span less than a kilometer. Is that correct? If the d values you entered were accidentally kilometers instead of meters, multiply by 1000 and try again. If in fact the profile does only go less than a kilometer, then great, you have done everything perfectly and I apologize for the interruption.')
end

%% Parse inputs: 

Nthick = 2.5; % By default l/H of 2.5 corresponds to a boxcar moving window of 10 ice thicknesses. 
FillEndpoints = false; 
w = ones(size(d)); 

tmp = strncmpi(varargin,'CouplingLength',3); 
if any(tmp)
   Nthick = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'endpoints',3); 
if any(tmp)
   if strcmpi(varargin{find(tmp)+1},'fill')
      FillEndpoints = true; 
   else 
      assert(strcmpi(varargin{find(tmp)+1},'shrink'),'Endpoint option can only be fill or shrink.') 
   end
end

tmp = strncmpi(varargin,'weights',3); 
if any(tmp)
   w = varargin{find(tmp)+1}; 
   assert(isequal(size(z),size(w)),'Dimensions of weights must match the dimensions of d, z, and H.')
end

%%

l = Nthick.*H; % works for scalar or vector Nthick

zs = nan(size(z)); 

for k = 1:length(z)
   
   % Generate an exponential weight curve centered on this point: 
   weight = exp(-abs(d - d(k))./l); 
   
   if FillEndpoints
   
      dlo = find(d<d(k) & weight<exp(-2),1,'last'); 
      dhi = find(d>d(k) & weight<exp(-2),1,'first'); 
      if ~isempty(dlo) & ~isempty(dhi)
         if ~any(isnan(l((dlo+1):(dhi-1))))
               zs(k) = wmean(z,weight.*w,'omitnan'); 
         end
      end
   else
      
      zs(k) = wmean(z,weight.*w,'omitnan'); 
   end
   
end

end

%% Subfunction: 

function M = wmean(A,weights,varargin)
% wmean computes the weighted average or weighted mean value. 
% 
%% Syntax
% 
%  M = wmean(A,weights)
%  M = wmean(...,'all')
%  M = wmean(...,'dim',dim)
%  M = wmean(...,'nanflag') 
% 
%% Description
% 
% M = wmean(A,weights) returns the weighted mean of the elements of A along the first
% array dimension whose size does not equal 1. Dimensions of A must match the 
% dimensions of weights. 
%   * If A is a vector, then wmean(A,weights) returns the weighted mean of the elements. 
%   * If A is a matrix, then wmean(A,weights) returns a row vector containing the weighted 
%     mean of each column.
%   * If A is a multidimensional array, then wmean(A,weights) operates along the first 
%     array dimension whose size does not equal 1, treating the elements as vectors.
%     This dimension becomes 1 while the sizes of all other dimensions remain the same.
% 
% M = wmean(...,'all') computes the weighted mean over all elements. (Requires
% Matlab R2018b or later.) 
% 
% M = wmean(...,'dim',dim) returns the mean along dimension dim. For example, if A 
% is a matrix, then wmean(A,weights,'dim',2) is a column vector containing the
% weighted mean of each row.
% 
% M = wmean(...,'nanflag') specifies whether to include or omit NaN values from the 
% calculation for any of the previous syntaxes. wmean(A,weights,'includenan') includes 
% all NaN values in the calculation while wmean(A,weights,'omitnan') ignores them.
% 
%% Examples 
% For examples type 
% 
%   cdt wmean
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas 
% Institute for Geophysics (UTIG). 
% 
% See also: mean and local. 

%% Initial error checks:

narginchk(2,Inf) 
assert(isequal(size(A),size(weights)),'Error: Dimensions of A and weights must agree.') 
assert(any(weights(:)~=0),'Error: all weights are zero and that does not make sense.') 

if (any(weights(:)<0))
    warning('Some of the weights are negative. I''ll do the math on the numbers you gave me, but I wonder if this is what you want to be doing?');
end

%% Parse optional inputs:

% Set defaults: 
dim = max([1 find(size(A)~=1,1,'first')]); % operate along the first non-singleton array (or dim 1 if find turns up empty). 

tmp = strncmpi(varargin,'dimension',3); 
if any(tmp)
   dim = varargin{find(tmp)+1}; 
   tmp(find(tmp)+1) = true; 
   varargin = varargin(~tmp); 
end

% It's possible the user wants the dimension to be 'all': 
tmp = strcmpi(varargin,'all'); 
if any(tmp)
   dim = varargin{tmp}; 
   varargin = varargin(~tmp); 
end
   
%% Perform mathematics: 

% Set weights to zero where they correspond to NaN
if any(strcmpi(varargin,'omitnan'))
   weights(isnan(A)) = NaN; 
end

M = sum(weights.*A,dim,varargin{:})./sum(weights,dim,varargin{:});
end