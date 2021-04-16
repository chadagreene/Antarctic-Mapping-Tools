%% |coord| documentation
% |coord| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017). Click <List_of_functions.html here>
% for a complete list of functions in AMT. 
% 
% The |coord| function returns coordinates from mouse clicks.  This function was written 
% for the Antarctic Mapping Toolbox.  
% 
% After calling coord, the current figure title will change to let you know
% you can start clicking.  Here are the commands: 
% 
%     Mouse click:   Places marker a marker and gets coordinates of the marker. 
%     Backspace:     Deletes previous marker like an Undo button. 
%     + or z:        Zooms in, centered on most recent click location (somewhat mediocre performance). 
%     - or x:        Zooms out, centered on most recent click location (somewhat mediocre performance).          
%     Escape key:    Terminates program without outputs. 
%     Return key:    Terminates program with data output.
%
%% Syntax
% 
%  coord
%  [lat,lon] = coord('geo')
%  [x,y] = coord('xy') 
%  [...] = coord(,...'MarkerProperty',MarkerValue)
%  [...] = coord(...,'keep') 
%  [...,h] = coord(...)
%
%% Description
% 
% |coord| without any inputs prints coordinates of mouse clicks to the Command Window
% in a two-column format,| [lat lon]| or |[x y]|.  Coordinates have the units of the 
% current axes unless |'geo'| or |'xy'| are specified.  
% 
% |[lat,lon] = coord('geo')| specifies that you want outputs in geographic coordinates, 
% even if you are clicking around on a polar stereographic x/y map. 
% 
% |[x,y] = coord('xy')| returns coordinates in polar stereographic (ps71) coordinates. 
% 
% |[...] = coord(,...'MarkerProperty',MarkerValue)| specifies marker properties for 
% when you're clickin' around on the map.  
% 
% |[...] = coord(...,'keep')| prevents coord from automatically deleting markers upon
% exiting the user interface.  
% 
% |[...,h] = coord(...)| returns a handle |h| of any plotted markers.   
% 
%% Example 1: With the Matlab's Mapping Toolbox
% Get some coordinates from clicking around in map coordinates. Below 
% I'm initializing a map with |bedmap2| then I'll click once on the
% Antarctic Peninsula and a second time in East Anarctica: 

bedmap2 gl 
[lat,lon] = coord  % (<--click around on the map and press enter when you're satisfied) 

%% 
% By default, the markers disappear after hitting Enter, but you can see from the coordinates 
% that I clicked here: 

plotm(lat,lon,'bo','markersize',16,'markerfacecolor','g')

%% Example 2: Also with Matlab's Mapping Toolbox: 
% This example is a little more involved. Return polar stereographic coordinates, 
% make the makers large red pentagrams, and use the |'keep'| option to retain the markers upon hitting Return:

modismoa 'pine island glacier'
[x,y] = coord('xy','marker','p','markerfacecolor','r','markersize',50,'keep')

%% Example 3: No Mapping Toolbox Required: 
% If you start in polar stereographic meters, coordinates are automatically
% returned in polar stereographic meters unless you specify |'geo'|:  

figure
bedmap2('patch coast','xy') 
xlabel 'eastings (m)' 
ylabel 'northings (m)' 
[x,y] = coord

%% 
% Alternatively, specify |'geo'| and click on the very same map get geo coordinates: 
 
[lat,lon] = coord('geo')

%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the University of 
% Texas at Austin's Institute for Geophysics (UTIG) in November 2015.