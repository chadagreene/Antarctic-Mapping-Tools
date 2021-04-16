function psax_m2km(varargin)
% psax_m2km converts polar stereographic axis tick labels from 
% meters to kilometers. It's just a clunky way to fix the ugly scientific notation
% that tends to result from plotting in ps71 cooridinates. 
% 
%% Syntax 
% 
%  psax_m2km
%  psax_m2km('reset') 
% 
%% Description  
% 
% psax_m2km called after a map has already been opened changes axis labels
% from meters to kilometers. This function uses the existing tick locations, 
% and prints the tick values divided by 1000. 
% 
% psax_m2km('reset') is a way to undo calling psax_m2km. 
% 
%% Example 
% 
% antbounds gl
% psax_m2km 
% xlabel 'these units are kilometers'
% ylabel 'these units are also kilometers'
% 
%% Usage note: 
% If you call psax_m2km, then pan or zoom, you'll have to do psax_m2km('reset'), 
% then call psax_m2km again. 
% 
%% Citing Antarctic Mapping Tools
% This function was developed for Antarctic Mapping Tools for Matlab (AMT). If AMT is useful for you,
% please cite our paper: 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. http://dx.doi.org/10.1016/j.cageo.2016.08.003
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
%% 
% 
% See also: plotps, pcolorps, scalebarps. 

%% Parse inputs 

ResetAx = false; % default
if nargin>0
   if strcmpi(varargin,'reset')
      ResetAx = true; 
   end
end

%% Do the work: 

ax = gca; 

if isequal(axis(ax),[0 1 0 1])
   error('The psax_m2km function can only be called after a map is initialized.') 
end

if ResetAx
   set(ax,'XTickLabelMode','auto')
   set(ax,'YTickLabelMode','auto')

else

   set(ax,'XTickLabel',num2str(get(ax,'XTick')'/1000))
   set(ax,'YTickLabel',num2str(get(ax,'YTick')'/1000))

end

end