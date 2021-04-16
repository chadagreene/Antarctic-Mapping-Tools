function [h,corelat,corelon] = corelabel(corename,varargin)
% CORELABEL labels Antarctic ice core locations on a map. Names and 
% locations are given by the 231 ice core locationslisted by iceREADER. 
% (http://www.icereader.org/icereader/listData.jsp)
% 
% Requires Matlab's Mapping Toolbox. 
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
% corelabel('CoreName')
% corelabel('CoreName','TextProperty',PropertyValue)
% corelabel('CoreName','Marker','MarkerStyle') 
% h = corelabel(...)
% [h,corelat,corelon] = corelabel(...)
% 
% 
%% Description 
% 
% corelabel('CoreName') labels an ice core location if a map is current.
% 
% corelabel('CoreName','TextProperty',PropertyValue) formats the label with
% text properties as name-value pairs. 
% 
% corelabel('CoreName','Marker','MarkerStyle') places a marker at the
% center location. 
% 
% h = corelabel(...) returns the handle h of the label. 
% 
% [h,corelat,corelon] = corelabel(...) returns geographic coordinates of
% the ice core. 
%
% 
%% Example: Label a couple of ice core locations:
% 
% load coast
% antmap
% patchm(lat,long,'y')
% 
% corelabel('Vostok')
% corelabel('Dome C','color','b','marker','b*')
% corelabel('byrd','rotation',45,'fontangle','italic')
% 
%
%% Author Info
% Created by Chad A. Greene, October 2013. 
% The University of Texas Institute for Geophysics.
% 
% See also textm, scarclick, scarlabel, coreloc, scarloc, and antmap.
%% Get location: 

[corelat,corelon] = coreloc(corename);

%% Place text label: 

% Apply the text label: 
h = textm(corelat,corelon,corename,'horizontalalignment','center');


%% Place a marker if requested:  

tmp = strcmpi(varargin,'marker'); 
if any(tmp)
    h(2) = plotm(corelat,corelon,varargin{find(tmp)+1});
    varargin{find(tmp)+1}=1; 
    varargin = varargin(~tmp); 
    set(h(1),'VerticalAlignment','bottom'); 
end

%% Format text and marker

% This is brute force, but it works:
for k = 1:length(h) 
    for k2 = 1:length(varargin)
        try
             set(h(k),varargin{k2},varargin{k2+1}); 
        end
    end
end

%% Clean up: 

% Delete output arguments if none are requested:
if nargout==0
    clear h corelat corelon 
end