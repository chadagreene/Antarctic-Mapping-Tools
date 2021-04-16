function [corelat,corelon] = coreloc(corename)
% CORELOC returns the latitude and longitude corresponding
% to any of the 231 Antarctic ice core locations listed by iceREADER.
% (http://www.icereader.org/icereader/index.jsp).  The corename 
% argument must be a string.
% 
% Requires the data file corenames.mat to be in a Matlab-findable
% directory.
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
%% Syntax and Description
% 
% [lat,lon] = coreloc('NameOfAnIceCore') returns the lat/lon location of an
% ice core. 
% 
%% Example: Where is Byrd Ice Core? 
%
% [lat,lon] = coreloc('Byrd')
% lat =
%   -80.0200
% lon =
%   -119.5200
%
%% Author Info
% Created by Chad A. Greene, October 2013 for distribution with
% The Bedmap2 Toolbox for Matlab. Migrated from the Bedmap2 Toolbox 
% to Antarctic Mapping Tools in August 2014. 
% 
% See also scarloc, scarclick, corelabel, scarlabel, and textm.

% The database is in corenames.mat:
try
    load corenames.mat; 
catch err
    error('MATLAB:cantFindCorenames','Can''t find corenames.mat.');
end

% Get lat/lon of feature: 
x=strcmpi(strtrim(corename),names); % logical array of matches corelat = lat(x); 
corelat = lat(x); 
corelon = lon(x);

% If the corename is not found, look for nearby matches and offer some
% helpful advice: 
if sum(x)==0
    fprintf(['"',corename,'" not found.']); fprintf('\n')
    
    nearbynames = strncmpi(strtrim(corename),names,2);
    
    if sum(nearbynames)>0
        fprintf(['Here are some available options beginning with "',corename(1:2),'": \n'])
        disp(names(nearbynames))
    else fprintf('Try typing "load corenames" to explore the available list of features. \n')
    end
    
    return

end

end

