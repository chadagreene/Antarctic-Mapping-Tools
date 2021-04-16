function [rowrange,colrange] = find2drange(X,Y,xi,yi,extraIndices) 
% find2drange returns matrix indices which encompass a range of xi,yi
% values. This function is typically used to clip a big 2D or 3D data 
% matrix (e.g., gridded DEM, climate data, etc.) to a smaller geographic
% range. 
% 
%% Syntax 
% 
%  [rowrange,colrange] = find2drange(X,Y,xi,yi)
%  [rowrange,colrange] = find2drange(X,Y,[],yi)
%  [rowrange,colrange] = find2drange(X,Y,xi,yi,extraIndices)
%  [rowrange,colrange] = find2drange(X,Y,xi,yi,[extraRows extraCols])
%  linearInd = find2drange(...) 
% 
%% Description 
% 
% [rowrange,colrange] = find2drange(X,Y,xi,yi) returns range of rows and
% columns of X and Y which fully encompass the ranges of values in xi and
% yi. X and Y must be matrices of equal dimensions. Dimensions of xi and yi
% do not need to match. 
%
% [rowrange,colrange] = find2drange(X,Y,[],yi) allows a full range of X
% values by entering an empty xi input. The same syntax can be used for yi
% to clip X and Y to only a range of xi values. 
% 
% [rowrange,colrange] = find2drange(X,Y,xi,yi,extraIndices) specifies an
% extra number of rows and columns to include in each direction. This may
% be helpful if you'll later do cubic or spline interpolation, which
% require points beyond the immediately-abutting X and Y values. If extra
% rows or columns result in indices less than one or exceeding the
% dimensions of X and Y, rowrange and colrange will be clipped to the range
% 1:size(X). 
%
% [rowrange,colrange] = find2drange(X,Y,xi,yi,[extraRows extraCols])
% specifies different values for extraRows and extraCols. Careful, because
% rows and columns are not x and y, in their units or their order! 
% 
% linearInd = find2drange(...) returns linear indices 
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
% 
% This function was written by Chad A. Greene of the University 
% of Texas at Austin's Institute for Geophysics (UTIG), Jan 2015. 
% http://www.chadagreene.com
% 
% See also: find, sub2ind, ind2sub.  

%% Initial input checks: 

narginchk(4,5) 
assert(isnumeric(X)==1,'X must be numeric.')
assert(ismatrix(X)==1,'This function only works for 1D or 2D X and Y arrays.') 
assert(size(X,1)==size(Y,1),'X and Y must be the same exact size.') 
assert(numel(X)==numel(Y),'X and Y must be the same exact size.') 
assert(isnumeric(xi)==1,'xi must be numeric.') 
assert(isnumeric(yi)==1,'yi must be numeric.') 

%% Parse inputs: 

extrarows = 0; 
extracols = 0; 
if nargin==5
    assert(isnumeric(extraIndices)==1,'extraIndices must be a scalar or two-element numeric array.') 
    extrarows = extraIndices(1); 
    assert(extrarows>=0,'extrarows must be a positive integer.') 
    assert(mod(extrarows,1)==0,'extrarows must be an integer.') 
    if isscalar(extraIndices)
        extracols = extrarows; 
    else
        assert(numel(extraIndices)==2,'extraIndices must be a scalar or two-element numeric array.')
        extracols = extraIndices(2); 
        assert(extracols>=0,'extracols must be a positive integer.') 
        assert(mod(extracols,1)==0,'extracols must be an integer.') 
    end
end

[rowsin,colsin] = size(X); 

% take entire x range if xi is not declared:
if isempty(xi)
    xi = [min(X(:)) max(X(:))]; % assume full X range
else
    xi = [min(xi(:)) max(xi(:))]; 

end

% take entire y range if yi is not declared:
if isempty(yi)
    yi = [min(Y(:)) max(Y(:))]; 
else
    yi = [min(yi(:)) max(yi(:))]; 
end


%% Find range: 

[rowi,coli] = find(X>=xi(1) & X<=xi(2) & Y>=yi(1) & Y<=yi(2)); 

% Range plus extra rows and columns: 
rowrange = (min(rowi)-1-extrarows):(max(rowi)+1+extrarows); 
colrange = (min(coli)-1-extracols):(max(coli)+1+extracols); 

% Clip range to ensure indices do not exceed dimensions of input X and Y:  
rowrange = rowrange(rowrange>=1 & rowrange<=rowsin); 
colrange = colrange(colrange>=1 & colrange<=colsin); 

% Return linear indices instead of subscripts if one output is requested 
if nargout==1 && colsin>1
    rowrange = sub2ind([rowsin colsin],rowrange,colrange); 
end

end