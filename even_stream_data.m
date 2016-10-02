function [xy, dist] = even_stream_data(xx, yy, uu, vv, dist_sep, dist_test, varargin)
% [xy, dist] = even_stream_data(xx, yy, uu, vv, dist_sep, dist_test)
% [xy, dist] = even_stream_data(___, Name, Value)
%
% Compute evenly-spaced streamlines with Jobar & Lefer algorithm [1].
% Always returns streamline points. Optionally, return arc length
% (distance along each streamline) and minimum distance to neighboring
% lines.
%
% Required Arguments:
%   xx, yy: Matrices or vectors, x-coord and y-coord. If matrices, the
%       size must match uu and vv. If vectors, xx must match the number of
%       columns in uu and vv, and yy must match the number of rows.
%   uu, vv: Matrices, vector field x-component and y-component
%   dist_sep: Scalar, minimum distance between seed point and streamlines,
%       as a fraction of the minimum domain width (i.e. min( range(xx(:)),
%       range(yy(:)) ) , if set to [], use default = 0.05
%   dist_test: Scalar, minimum distance between streamlines, as a fraction
%       of the minimum domain width, id set to [], use default = 0.02
%
% Optional Parameters (Name - Value):
%   'StepSize': Scalar, streamline step size, see streamline() for details,
%       default = 0.1
%   'MaxNumVertex', Scalar, maximum number of points in each stream line,
%       see streamline() for details, default = 10,000
%   'Verbose': set true to enable verbose messages,
%       default = false
%
% Return:
%   xy: Matrix, [x, y] coordinates for stream line points, each row is a
%       point, individual lines are separated by NaNs.
%   dist: Vector, minimum distance to neighboring stream lines, only
%       computed only if nargout>1.
%
% References:
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
% %

%% parse and check inputs

% handle optional arguments
if isempty(dist_sep); dist_sep = 0.05; end
if isempty(dist_test); dist_test = 0.02; end

% handle optional name-value parameters
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;
parser.addParameter('StepSize', 0.1);
parser.addParameter('MaxNumVertex', 10000);
parser.addParameter('Verbose', false);
parser.parse(varargin{:});
step_size = parser.Results.StepSize;
max_num_vertex = parser.Results.MaxNumVertex;
verbose = parser.Results.Verbose;

% convert coordinate vectors to grids, if indicated
if isvector(xx) && isvector(yy)
    [xx, yy] = meshgrid(xx, yy);
end

% check for sane values
validateattributes(xx, {'numeric'}, {'nonempty'}, mfilename, 'xx');
validateattributes(yy, {'numeric'}, {'size', size(xx)}, mfilename, 'yy');
validateattributes(uu, {'numeric'}, {'size', size(xx)}, mfilename, 'uu');
validateattributes(vv, {'numeric'}, {'size', size(xx)}, mfilename, 'vv');
validateattributes(dist_sep, {'numeric'}, {'scalar', '>' 0, '<', 1}, ...
    mfilename, 'dist_sep');
validateattributes(dist_test, {'numeric'}, {'scalar', '>' 0, '<=', dist_sep}, ...
    mfilename, 'dist_test');
validateattributes(step_size, {'numeric'}, {'scalar', 'positive'}, ...
    mfilename, 'step_size');
validateattributes(max_num_vertex, {'numeric'}, {'scalar', 'positive', 'integer'}, ...
    mfilename, 'max_num_vertex');
validateattributes(verbose, {'numeric', 'logical'}, {'scalar', 'binary'}, ...
    mfilename, 'verbose');

% convert distance parameters from fractions to data units
min_domain_range = min( range(xx(:)), range(yy(:)) ); 
dist_sep = dist_sep*min_domain_range;
dist_test = dist_test*min_domain_range;

%% Compute stream line points

% disable nuisance warning(s)
% ...Delaunay triangulation drops duplicate points from streamlines, OK
warning('off', 'MATLAB:delaunayTriangulation:DupPtsWarnId');

% select random non-NaN grid point for initial seed
while 1
    kk = randi([1, numel(xx)]);
    if ~isnan(uu(kk)) && ~isnan(vv(kk))
        break
    end
end
seed_xy = [xx(kk), yy(kk)];

% add first streamline to triangulation and length list
[stream_xy, ~] = get_streamline(xx, yy, uu, vv, seed_xy, step_size);
stream_tri = delaunayTriangulation(stream_xy);
stream_len = size(stream_tri.Points,1);

% create seed point candidate queue
seed_queue{1} = get_seed_candidates(stream_xy, dist_sep);

% check all seed candidates
while ~isempty(seed_queue)
    
    % pop seed candidates from queue
    seed_xy = seed_queue{1};
    seed_queue(1) = [];
    
    % check each seed candidate in random order
    for ii = randperm(length(seed_xy))
        
        % skip if candidate istoo close to any streamline point
        [~, d_min] = nearestNeighbor(stream_tri, seed_xy(ii,:));
        if d_min < dist_sep
            continue
        end
        
        % create new streamline, skip if empty
        [stream_xy, seed_idx] = get_streamline(xx, yy, uu, vv, seed_xy(ii,:), step_size);
        if size(stream_xy,1) < 2
            continue
        end
        
        % trim new streamline
        for jj = seed_idx:size(stream_xy,1)
            [~, d_min] = nearestNeighbor(stream_tri, stream_xy(jj,:));
            if d_min < dist_test
                break
            end
        end
        for kk = seed_idx:-1:1
            [~, d_min] = nearestNeighbor(stream_tri, stream_xy(kk,:));
            if d_min < dist_test
                break
            end
        end
        stream_xy = stream_xy(kk:jj, :);
        
        % add streamline to triangulation and line length list
        len0 = size(stream_tri.Points, 1);
        stream_tri.Points = [stream_tri.Points; stream_xy];
        stream_len(end+1) = size(stream_tri.Points, 1)-len0; %#ok!
        
        % add seed candidate points to queue
        seed_queue{end+1}  = get_seed_candidates(stream_xy, dist_sep); %#ok!
        
        if verbose
            num_lines = length(stream_len);
            num_seeds = 0;
            for pp = 1:length(seed_queue)
                num_seeds = num_seeds+length(seed_queue{pp});
            end
            fprintf('%s: xy: %d lines, %d seed candidates\n', ...
                mfilename, num_lines, num_seeds);
        end
        
    end
end

% extract stream line points for output as xy
num_pts = size(stream_tri.Points, 1);
num_lines = length(stream_len);
xy = nan(num_pts+num_lines-1, 2);
ii0 = 1;
jj0 = 1;
for kk = 1:num_lines
    ii1 = ii0+stream_len(kk)-1;
    jj1 = jj0+stream_len(kk)-1;
    xy(jj0:jj1,:) = stream_tri.Points(ii0:ii1,:);
    ii0 = ii1+1;
    jj0 = jj1+2;
end

% enable nuisance warning(s)
warning('on', 'MATLAB:delaunayTriangulation:DupPtsWarnId');

%% Compute distance to neighbors

dist = nan(size(xy, 1), 1);
if nargout>1
    kk = 1;
    for ii = 1:num_lines
        if verbose
            fprintf('%s: dist: line %d of %d\n', mfilename, ii, num_lines);
        end
        stream_xy = stream_tri.Points(1:stream_len(ii), :);
        stream_tri.Points(1:stream_len(ii), :) = [];
        [~, stream_dist] = nearestNeighbor(stream_tri, stream_xy);
        dist(kk:kk+stream_len(ii)-1) = stream_dist;
        kk = kk+stream_len(ii)+1;
        stream_tri.Points = [stream_tri.Points; stream_xy];
    end
end

function seed_xy = get_seed_candidates(xy, buf_dist)
%
% Compute the location of stream line seed point candidates that lie at a
% distance 'buf_dist' along a normal vector at each point in 'xy'
%
% Arguments:
%   xy: Matrix, [x,y] coordinates of points along a streamline in rows
%   buf_dist: Scalar, buffer distance between xy and seed point candidates
%   seed_xy, y_seed: Vectors, canditate seed points
% %

% get unit normal vectors at segment midpoints
tangent = diff(xy);
normal = [tangent(:,2), -tangent(:,1)];
normal = bsxfun(@rdivide, normal, sqrt(sum(normal.*normal, 2)));
midpoint = xy(1:end-1, :)+0.5*tangent;

% get candidates offset buf_dist in positive and negative normal direction
seed_xy = [midpoint + buf_dist*normal; midpoint - buf_dist*normal];


function [stream_xy, seed_idx] = get_streamline(xx, yy, uu, vv, seed_xy, step_size)
%
% Compute streamline in both directions starting at seed point
%
% Arguments:
%   xx, yy:
%   uu, vv:
%   seed_xy: Vector, [x, y] coordinates of seed point
%   stream_xy : Matrix, stream line x- and y-coordinates in rows, returns [] if
%       stream line has zero length
%   seed_idx: Scalar, index of seed point in output streamline
% %

fwd = stream2(xx, yy, uu, vv, seed_xy(1), seed_xy(2), step_size);
xy_fwd = fwd{1}(~any(isnan(fwd{1}), 2), :); % drop NaN rows
has_fwd = size(xy_fwd,1) > 1;

rev = stream2(xx, yy, -uu, -vv, seed_xy(1), seed_xy(2), step_size);
xy_rev = rev{1}(~any(isnan(rev{1}), 2), :); % drop NaN rows
has_rev = size(xy_rev,1) > 1;

if has_fwd && has_rev
    stream_xy = [xy_rev(end:-1:2, :); xy_fwd];
    seed_idx = size(xy_rev,1);
elseif has_rev
    stream_xy = xy_rev;
    seed_idx = 1;
elseif has_fwd
    stream_xy = xy_fwd;
    seed_idx = 1;
else
    stream_xy = [];
    seed_idx = [];
end
