%% |scarloc| documentation  
% |scarloc| is part of Antarctic Mapping Tools for Matlab (Greene et al., 2017).  <List_of_functions.html Click here
% for a complete list of functions in AMT.>
% 
% The |scarloc| function returns coordinates of any of the 25,601 features identified by the Scientific Committee 
% on Antarctic Research (<http://www.scar.org/cga SCAR>). 
% 
%% Syntax 
% 
%  [lat,lon] = scarloc(FeatureNames) 
%  latlon = scarloc(FeatureNames) 
%  [x,y] = scarloc(FeatureNames,'xy') 
%  xy = scarloc(FeatureNames,'xy')
%  [...] = scarloc(FeatureNames,'xy','km')
% 
%% Description 
% 
% |[lat,lon] = scarloc(FeatureNames)| returns the geo coordinate location(s) of an
% almost any feature(s) in Antarctica. |FeatureNames| can be string or cell array of 
% multiple locations. 
% 
% |latlon = scarloc('Feature Name')| returns geo coordinate location(s) of an
% almost any feature(s) in Antarctica. If only one output argument is used,
% column 1 corresponds to latitude and column 2 corresponds to longitude. 
% 
% |[x,y] = scarloc(FeatureNames,'xy')| returns polar stereographic (true latitude 
% 71°S) eastings and northings in meters.  
% 
% |[y = scarloc(FeatureNames,'xy')| returns polar stereographic (true latitude 
% 71°S) eastings and northings in meters. If only one output argument is used,
% column 1 corresponds to eastings and column 2 corresponds to northings. 
% 
% |[...] = scarloc(FeatureNames,'xy','km')| returns polar stereographic kilometers 
% instead of the default meters. 
% 
%% Example: Where is McMurdo Station? 

[lat,lon] = scarloc('McMurdo Station')

%% Example 2: Two-column output
% If two outputs are requested as in Example 1, outputs are lat and lon. If
% one output is requested, latitudes populate column 1 and longitudes
% populate column 2: 

latlon = scarloc('mcmurdo station')

%% Example 3: Multiple locations 
% Multiple locations may be entered as a cell array of strings: 

places = {'Byrd Camp';'Casey Station';'Pine Island Glacier'}; 
scarloc(places)

%% Example 4: Polar stereographic coordinates 
% If you want polar stereographic eastings and northings, include the
% |'xy'| tag: 

scarloc(places,'xy')

%% Example 5: Polar stereographic kilometers 
% Sometimes you want polar stereographic kilometers, but dividing meters by
% a thousand feels like far too much of a burden.  Fortunately,
% |scarloc|'ll do it for you: 

scarloc(places,'xy','km')

%% Known Issues
% Unfortunately, scarloc and scarlabel are plagued by problems with special
% characters. This includes letters with accents and umlauts. My sincerest
% apologies to our French, Norwegian, and just about everything aside from
% American friends.  
% 
%% Cite this dataset as: 
% Secretariat SCAR (1992, updated 2015). Composite Gazetteer of Antarctica,
% Scientific Committee on Antarctic Research. GCMD Metadata 
% <http://gcmd.nasa.gov/records/SCAR_Gazetteer.html>
% 
%% Citing AMT
% If this function or any other part of Antarctic Mapping Tools is useful for you, please
% cite the paper that describes AMT.  
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% _Computers & Geosciences_. 104 (2017) pp.151-157. <http://dx.doi.org/10.1016/j.cageo.2016.08.003 doi:10.1016/j.cageo.2016.08.003>.
% 
%% Author Info
% This function and supporting documentation were written by
% <http://www.chadagreene.com Chad Greene> of the University of Texas at
% Austin's Institute for Geophysics (UTIG). Feel free to contact me if you
% have any questions or comments. 