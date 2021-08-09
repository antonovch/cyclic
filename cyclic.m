classdef cyclic
% cyclic object is a substitution for Matlab's grid tensor. For
% multi-dimensional grids, their sizes can grow to large values, while the
% unique values of each represent a vector. This class implements correct
% indexing, while storing this minimal amount of information, along with
% desctiptive metadata, namely, tensor's size and dimension of variation.
% NOTE: under development
    properties
        sz (1,:) {mustBeInteger, mustBeNonnegative} % size vector of the tensor
        dim (1,:) {mustBeInteger, mustBePositive} % number(s) of the dimension(s) of variation
        data (:,1) % values along the dimension(s) dim
    end
    
    methods
        function obj = cyclic(varargin)
            switch nargin
                case 1
                    data = varargin{1};
                    sz = size(data);
                    N = ndims(data);
                    I = [ones(1,N); repmat(cumprod([1 sz(1:end-1)]),max(sz)-1,1)];
                    I = cumsum(I,1);
                    idx = sub2ind(size(I), sz, 1:N);
                    I = min(I, I(idx));
                    dim = find(any(diff(data(I),1)));
                    if isempty(dim) % all identical datapoints, use single scalar
                        dim = find(sz == 1, 1);
                        if isempty(dim)
                            dim = N + 1;
                        end
                        % TODO: write error message, error for empty dim,
                        % check for multiple dimensions of variation,
                        % vecorize stuff, think of cell arrays
                        data = data(1);
                    elseif isequal(dim,1:N)
                        error('cyclic:cyclic', 'Input tensor does not appear to be mesh-like.')
                    else % extract variable data
                        I = I(1:max(sz(dim)),dim);
                        data = data(I);
                        mask = repmat((1:size(data,1)).',1,size(data,2));
                        data(mask > sz(dim)) = NaN;
                    end
                case 3
                    [data, sz, dim] = varargin{:};
                otherwise
                    error('Wrong number of input arguments')
            end
            
%             % if multiple consecutive dimensions are specified, they are
%             % flattened. Equivalent to sz = size(data) = [1 2 3 4 5 6],
%             % dim = 3:5 => data = reshape(data, 1, 2, [], 6).
%             if numel(dim) > 1
%                 assert(all(diff(dim) == 1), 'Dimensions have to be consecutive.')
%                 sz = [sz(1:dim(1)-1), prod(sz(dim)), sz(dim(end)+1)];
%             end
            obj.sz = sz;
            obj.dim = dim;
            obj.data = data;
        end
        
        function varargout = subsref(self, s)
            switch s(1).type
                case '()' % TODO slicing
                    N = length(s(1).subs);
                    if N > self.dim 
                        % subscript for the dimension of interest is given
                        self = self.data(s(1).subs{self.dim});
                    else % needed index is in the flattened dimensions
                        I = s(1).subs{N};
                        if strcmp(I,':')
                            error('Slicing is not supported');
                        else
                            subs = cell(1, numel(self.sz)-N+1); % missing dimensions
                            [subs{:}] = ind2sub(self.sz(N:end), I);
                            self = self.data(subs{self.dim-N+1});
                        end
                    end
                    s = s(2:end);
                case '{}'
                    varargout = subsref(self, struct('type','()','subs',{s(1).subs}));
                    if numel(s) > 1
                        if numel(varargout) == 1
                            varargout = {builtin('subsref', varargout{1}, s(2:end))};
                        else
                            error(['Expected one output from a curly brace ',...
                                'or dot indexing expression, but there were 2 results.']);
                        end
                    end
                    return
            end
            if isempty(s)
                varargout = {self};
            else
                varargout = cell(1,max(1, nargout));
                [varargout{:}] = builtin('subsref', self, s);
            end
        end
        
        function M = full(self)
            self.sz(self.dim) = 1;
            if self.dim == 1
                M = self.data;
            else
                M = reshape(self.data, [ones(1, self.dim-1) numel(self.data)]);
            end
            M = repmat(M, self.sz);
        end
        
        function varargout = size(self, varargin)
            if isempty(varargin)
                if nargout < 2
                    varargout = {self.sz};
                else
                    varargout = arrayfun(@(n)self.sz(n), 1:nargout, 'UniformOutput', false);
                end
            else
                if nargout < 2
                    varargout = {self.sz([varargin{:}])};
                else
                    idx = [varargin{:}]; % to support both {[1 2 3]} and {1,2,3}
                    assert(numel(idx) >= nargout, 'Too many output arguments');
                    varargout = arrayfun(@(n)self.sz(idx(n)), 1:nargout, 'UniformOutput', false);
                end
            end
        end
        
        function tf = isnan(self)
            tf = isnan(self.data);
        end
        
        %%% messes up subsref when asking for a property, but really we
        %%% shouldn't ask for it directly, and numel is needed for
        %%% numerical array compatibility. Use get methods
        function n = numel(self, varargin)
            if isempty(varargin)
                n = prod(self.sz);
            else
                n = numel(subsref(self,struct('type','()','subs',{varargin})));
            end
        end
        
    end
end