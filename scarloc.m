function varargout = scarloc(featurename,varargin)
% SCARLOC returns coordinates of any of the 25,601 features identified by the Scientific Committee 
% on Antarctic Research (SCAR). http://www.scar.org/cga
% 
% Requires the data file scarnames.mat to be in a Matlab-findable
% directory.
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
% [lat,lon] = scarloc(FeatureNames) returns geographic coordinate location(s) of an
% almost any feature(s) in Antarctica. FeatureNames can be string or cell array of 
% multiple locations. 
% 
% latlon = scarloc('Feature Name') returns the geographic coordinate location(s) of an
% almost any feature(s) in Antarctica. If only one output argument is used,
% column 1 corresponds to latitude and column 2 corresponds to longitude. 
% 
% [x,y] = scarloc(FeatureNames,'xy') returns polar stereographic (true latitude 
% 71°S) eastings and northings in meters.  
% 
% [y = scarloc(FeatureNames,'xy') returns polar stereographic (true latitude 
% 71°S) eastings and northings in meters. If only one output argument is used,
% column 1 corresponds to eastings and column 2 corresponds to northings. 
% 
% [...] = scarloc(FeatureNames,'xy','km') returns polar stereographic kilometers 
% instead of the default meters. 
% 
%% Example: Where is McMurdo Station? 
%
% [lat,lon] = scarloc('McMurdo Station')
% lat =
%   -77.8478
% lon =
%   166.6683
% 
%% Example 2: Same as above, with only one output: 
% 
% latlon = scarloc('mcmurdo station')
% latlon =
%   -77.8478  166.6683
% 
%% Example 3: Multiple locations: 
% 
% scarloc({'Byrd Camp';'Casey Station';'Pine Island Glacier'})
% ans =
%   -80.0833 -119.5333
%   -66.2833  110.5194
%   -75.1667 -100.0000
% 
%% Example 4: Multiple locations in polar stereographic coordinates: 
% 
% scarloc({'Byrd Camp';'Casey Station';'Pine Island Glacier'},'xy')
% ans =
%    1.0e+06 *
%    -0.9397   -0.5324
%     2.4470   -0.9159
%    -1.5958   -0.2814
% 
%% Example 5: As above, but polar stereographic kilometers: 
% 
% scarloc({'Byrd Camp';'Casey Station';'Pine Island Glacier'},'xy','km')
% ans =
%        -939.73       -532.39
%        2447.05       -915.86
%       -1595.76       -281.38
% 
%% Known Issues
% 
% Unfortunately, scarloc and scarlabel are plagued by problems with special
% characters. This includes letters with accents and umlauts. My sincerest
% apologies to our French, Norwegian, and just about everything aside from
% American friends.  
% 
%% Cite this dataset as: 
% 
% Secretariat SCAR (1992, updated 2015). Composite Gazetteer of Antarctica,
% Scientific Committee on Antarctic Research. GCMD Metadata 
% <http://gcmd.nasa.gov/records/SCAR_Gazetteer.html>
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
% Created by Chad A. Greene of the University of Texas Institute 
% for Geophysics, August 2013.
% 
% Updated January 2015 to include new entries in the SCAR database.
% Updated July 2015 to allow polar stereographic units out. 
% 
% See also scarlabel, scarclick, coreloc, corelabel, and textm.

%% Input check:

narginchk(1,4)
assert(isnumeric(featurename)~=1,'Feature name cannot be numeric.')

OfferHelp = true; 
tmp = strncmpi(varargin,'NoWarning',5);
if any(tmp) 
    OfferHelp = false; 
end

pscoords = false; % return geo coordinates by default
if any(strcmpi(varargin,'xy'))
    pscoords = true; 
    
    kmout = false; % if polar stereographic, meters out by default
    if any(strcmpi(varargin,'km'))
        kmout = true; 
    end
end


%% Load data: 

% The database is in scarnames.mat:
try
    load scarnames.mat lat lon names;  
catch err
    error('MATLAB:cantFindScarnames','Can''t find scarnames.mat.');
end

%% Preallocate output: 

if ~ischar(featurename)
    featurelat = NaN(length(featurename),1); 
    featurelon = featurelat; 
else
    featurelat = NaN; 
end


%% Look for each feature name: 

for k = 1:length(featurelat)
    
    % Get lat/lon of feature: 
    if ~ischar(featurename)
        [x,NearbyNames] = strlookup(featurename{k},names); 
    else
        [x,NearbyNames] = strlookup(featurename,names); 

    end

    % If the featurename is not found, look for nearby matches and offer some
    % helpful advice: 

    if isempty(x) && OfferHelp 
        clear featurelat featurelon
        
        if ~ischar(featurename)
        	featurename = featurename{k};
        end
        
        fmsg{1} = ['"',featurename,'" not found.'];
        fmsg{2} = ['Are you sure that "',featurename,'" exists in Antarctica?'];
        fmsg{3} = 'Did a cat walk across your keyboard?';
        fmsg{4} = 'This is the real reason one shouldn''t text and drive. Check your spelling and try again.'; 
        fmsg{5} = 'Now you''re just making things up.'; 
        fmsg{6} = ['SCAR has identified more than 25,000 features in Antarctica, but "',featurename,'" is not one of them.'];
        fmsg{7} = ['Can''t find "',featurename,'."'];
        fmsg{8} = ['"',featurename,'" may exist somewhere in the world, but you won''t find it in Antarctica.'];
        fmsg{9} = ['It is possible that Robert F. Scott named something in Antarctica "',featurename,'," but if he did there are no records of it.']; 
        fmsg{10} = ['You must be thinking of ',featurename,', Kansas, because ',featurename,', Antarctica does not exist.'];
        fmsg{11} = ['Sure, they used to just call it ',featurename,', but not anymore, what with political correctness and all.'];
        fmsg{12} = ['"',featurename,'" is an interesting combination of letters, but I don''t think it''s any place in Antarctica.'];
        fmsg{13} = ['The great Wayne Cochran once sang, "Where oh where can my ',featurename,' be?"  Because it''s not in Antarctica.']; 
        fmsg{14} = ['I''m pretty sure it is in violation of the Antarctic Treaty to refer to any place as "',featurename,'."'];
        fmsg{15} = ['"',featurename,'" does not match any entries in the SCAR database.'];
        fmsg{16} = ['Science is all about formality, so the bigwigs will surely look down their noses at such colloquial jargon as "',featurename,'."']; 
        fmsg{17} = ['My doctor said I need to get my ',featurename,' removed.'];
        fmsg{18} = 'Frostbitten Antarctic researcher mistypes again.';
        fmsg{19} = 'This may be an issue of American English versus British English.'; 
        fmsg{20} = ['Antarctica''s a strange place, but it''s not science fiction. Verify that "',featurename,'" actually exists.'];
        fmsg{21} = ['What''s in a name? I''ll tell you what''s in a name: That which you call "',featurename,'" by any other name may actually exist in Antarctica.'];
        fmsg{22} = ['Did John Carpenter tell you''ll find "',featurename,'" in Antarctica?'];
        fmsg{23} = ['You know, some folks say glaciology is a shrinking field, but I say things are just heating up. In other news, "',featurename,'" does not exist.'];
        fmsg{24} = ['You''re a glaciologist?  Isn''t that a slow-moving field?  Also, I have to tell you, I can''t seem to find any record of "',featurename,'".'];
        fmsg{25} = ['Amazing glaciology, how sweet the sound... "',featurename,'" once was lost, and still has not been found.'];

        
        rngstart = rng; % get initial rng setting before changing it temporarily. 
        rng('shuffle'); 
        fprintf(fmsg{randi(length(fmsg))}); fprintf('\n')
        rng(rngstart); % returns to original rng settings. 


        if ~isempty(NearbyNames)
           disp('Here are the best matches I can find:')
           disp(NearbyNames)
        else fprintf('Try typing "load scarnames" to explore the available list of features. \n')
        end
        return
    end
      
    if isempty(x) 
        featurelat = NaN; 
        featurelon = NaN; 
    else
        featurelat(k) = lat(x);
        featurelon(k) = lon(x);
    end
end

% Convert to polar stereographic coordinates: 
if pscoords
    [out1,out2] = ll2ps(featurelat,featurelon); 
    
    if kmout
        out1 = out1/1000; 
        out2 = out2/1000; 
    end
else
    out1 = featurelat; 
    out2 = featurelon; 
end

% Returning only latitude or only x would not make any sense, so if no outputs are
% requested, or if only one output is requested, return as a lat column and
% lon column or [x y] 
if nargout < 2
    varargout{1} = [out1 out2]; 
else
    varargout{1} = out1; 
    varargout{2} = out2; 
end

end





%% strlookup
% The strlookup function follows: Documentation for strlookup can be found 
% at http://www.mathworks.com/matlabcentral/fileexchange/47577-strlookup/content/strlookup_documentation/html/strlookup_documentation.html
function [ind,CloseNames] = strlookup(string,list,varargin)
% STRLOOKUP uses strcmp or strcmpi to return indices of a list of strings 
% matching an input string. If no matches are found, close matches are
% suggested.
% 
%% Syntax 
% 
% ind = strlookup('string',list)
% ind = strlookup(...,'CaseSensitive')
% ind = strlookup(...,'threshold',ThresholdValue)   
% [ind,CloseNames] = strlookup(...) 
% 
% 
%% Description 
% 
% ind = strlookup('string',list) returns indices ind corresponding to
% cell entries in list matching 'string'. 
% 
% ind = strlookup(...,'CaseSensitive') performs a case-sensitive
% strlookup. 
% 
% ind = strlookup(...,'threshold',ThresholdValue) declares a threshold
% value for close matches. You will rarely (if ever) need to use this.
% The ThresholdValue is a metric of how closely matches should be when
% offering suggestions. Low threshold values limit suggested matches to a
% shorter list whereas high thresholds expand the list size. By default,
% the threshold starts at 1.5, then increases or decreases depending on how
% many close matches are returned.  If fewer than 3 close matches are
% returned, the threshold is increased and it looks for more matches. If
% more than 10 close matches are found, the threshold is tightened
% (reduced) until fewer than 10 matches are found. 
% 
% [ind,CloseNames] = strlookup(...) suppresses command window output if
% no exact match is found, and instead returns an empty matrix ind and a
% cell array of close matches in names. If exact match(es) is/are found, 
% ind will be populated and CloseNames will be empty. 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas 
% at Austin Institute for Geophysics (UTIG), August 2014. 
% 
% Updated January 2015 to include close alphabetical matches and fixed
% an input check based on FEX user Bryan's suggestion. Thanks Bryan. 
% 
% See also strcmp, strcmpi, strfind, regexp, strrep. 

%% Input checks: 

assert(isnumeric(string)==0&&isnumeric(list)==0,'Inputs must be strings.') 

% If user accidentally switches order string and list inputs, fix the order: 
if ischar(list) && iscell(string) 
    tmplist = list; 
    list = string; 
    string = tmplist; 
end

if strcmpi(string,'recursion')
    disp('Did you mean recursion?') 
    return
end

% Allow user to declare match threshold with a name-value pair: 
threshold = []; 
tmp = strncmpi(varargin,'thresh',6); 
if any(tmp)
    threshold = varargin{find(tmp)+1}; 
    assert(isscalar(threshold)==1,'Threshold value must be a scalar.')
    assert(threshold>=0,'Threshold value cannot be negative.')
end

%% Find matches:

% Find case-insensitive matches unless 'CaseSensitive' is requested by user:
if nargin>2 && any(strncmpi(varargin,'case',4))
    TF = strcmp(list,string); 
else
    TF = strcmpi(list,string); 
end


%% Look for close matches if no exact matches are found: 

CloseNames = []; 

if sum(TF)==0
    % Thanks to Cedric Wannaz for writing this bit of code. He came up with a
    % quite clever solution wherein the spectrum of an input string is compared to  
    % spectra of available options in the input list. 
    
    % Define spectrum function: 
    spec = @(name) accumarray(upper(name.')-31, ones(size(name)), [60 1]);
    
    % Get spectrum of input string:
    spec_str = spec(string); 
    
    % Compare spec_str to spectra of all strings available in the list:
    spec_dist = cellfun(@(name) norm(spec(name)-spec_str), list);
    
    % Sort by best matches: 
    [sds,idx] = sort(spec_dist) ;

    % Find list items that closely match input string by spectrum: 
    % If the user has not declared a hard threshold, start with a threshold
    % of 1.5 and then adjust threshold dynamically if there are too many or
    % too few matches:
    if isempty(threshold)
        threshold = 1.5; 
        closeSpectralInd = idx(sds<=threshold);
        
        % If there are more than 10 close matches, try a smaller threshold:
        while length(closeSpectralInd)>10
            threshold = 0.9*threshold;
            closeSpectralInd = idx(sds<=threshold); 
            if threshold<.05
                break
            end
        end
        
        % If there are fewer than 3 close matches, try relaxing the threshold: 
        while length(closeSpectralInd)<3
            threshold = 1.1*threshold;
            closeSpectralInd = idx(sds<=threshold); 
            if threshold>10
                break
            end
        end
        
    else
        % If user declared a hard threshold, stick with it:
        closeSpectralInd = idx(sds<=threshold);
    end
    
    % Check for matches alphabetically by seeing how many first-letter
    % matches there are. If there are more than 4 first-letter matches, 
    % see how many first-two-letter matches there are, and so on:
    for n = 1:length(string)
        closeAlphaInd = find(strncmpi(list,string,n));
        if length(closeAlphaInd)<5
            break
        end
    end
    
    
    % Names of close matches: 
    CloseNames = list(unique([closeSpectralInd;closeAlphaInd]));
    ind = []; 
    
    if isempty(CloseNames)
        disp(['String ''',string,''' not found and I can''t even find a close match. Make like Santa and check your list twice.'])
    else
        if nargout<2
            disp(['String ''',string,''' not found. Did you mean...'])
            disp(CloseNames); 
        end
        if nargout==0
            clear ind
        end
    end
    
else
    ind = find(TF); 
end

end
