function [xs, ys, ls, ds] = ...
    even_stream2(xx, yy, uu, vv, d_sep, d_test, step_size)
%
% Compute evenly-spaced streamlines with Jobar & Lefer algorithm (ref 1).
% Results can be used to plot using ____, ____, ...
%
% Arguments:
%
%   xx, yy:
%
%   uu, vv: 
%
%   d_sep:
%
%   d_test:
%
%   step_size:
%
%   xs, ys: Vectors, x-coord and y-coord for stream line points, individual
%       lines are separated by NaNs
%
%   ls: Vector, length along streamline for stream line points, individual
%       lines are separated by NaNs
%
%   ds: Vector, distance to nearest neighboring streamline for stream line
%       points, individual lines are separated by NaNs
%  
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
% %

%% initialize

% gather some basic coordinate info
x_min = min(xx(:));
x_max = max(xx(:));
x_rng = x_max-x_min; 
y_min = min(yy(:));
y_max = max(yy(:));
y_rng = y_max-y_min;

% get seed point at random (non-NaN) point
u0 = NaN;
v0 = NaN;
while isnan(u0) || isnan(v0)
    x0 = x_min+rand(1)*x_rng;
    y0 = y_min+rand(1)*y_rng;
    u0 = interp2(xx, yy, uu, x0, y0);
    v0 = interp2(xx, yy, vv, x0, y0);
end

% add first stream line to triangulation
[stream_xy, stream_idx] = get_streamline(xx, yy, uu, vv, x0, y0, step_size);
stream_tri = delaunayTriangulation(stream_xy);

% create seed point candidate queue 
seed_queue{1} = get_seed_candidates(stream_xy, d_sep);

keyboard

%% main loop

% TODO: merge x_ and y_ variables
% TODO: dispense with helper functions

d_sep_sq = d_sep*d_sep;
d_test_sq = d_test*d_test;

num_lines = 0;
while ~isempty(x_queue)
    
    % pop seed candidates from queue
    x_seed = x_queue{1}; x_queue(1) = [];
    y_seed = y_queue{1}; y_queue(1) = [];    
    
    % check each seed candidate in random order
    for ii = randperm(length(x_seed))
        
        % compute minimum distance to all streamline points
        dxy = bsxfun(@minus, [x_seed(ii), y_seed(ii)],  [x_line, y_line]);
        d_min_sq = min(sum(dxy.*dxy, 2));
        
        
        % create new streamline
        if d_min_sq >= d_sep_sq
            
            [x_new, y_new, seed_idx] = get_streamline(...
                xx, yy, uu, vv, x_seed(ii), y_seed(ii), step_size);
            
            if ~isempty(x_new) 
            
                % trim new streamline
                for kk = seed_idx:length(x_new)
                    if ~dist_gte(d_test_sq, x_new(kk), y_new(kk), x_line, y_line)
                        x_new(kk:end) = [];
                        y_new(kk:end) = [];
                        break
                    end
                end
                for kk = seed_idx:-1:1
                    if ~dist_gte(d_test_sq, x_new(kk), y_new(kk), x_line, y_line)
                        x_new(1:kk) = [];
                        y_new(1:kk) = [];
                        break
                    end
                end
                
                % add seed candidate points to queue
                [x_queue{end+1}, y_queue{end+1}] = ...
                    get_seed_candidates(x_new, y_new, d_sep); %#ok!
                
                % add streamline to neighbor index
                k_new = xy_to_k(x_new, y_new);
                for kk = unique(k_new)'
                    nbr{kk} = [nbr{kk}; find(k_new==kk)+length(x_line)];
                end
                
                % add trimmed streamline to list
                x_line = [x_line; NaN; x_new]; %#ok!
                y_line = [y_line; NaN; y_new]; %#ok!
                
                % tally and report
                num_lines = num_lines+1;
                fprintf('# streamlines: %d\n', num_lines);                
            end
        end
    end
end

%% prepare outputs

xs = x_line; 
ys = y_line;
ls = [];
ds = [];

%<DEBUG>
keyboard
%</DEBUG>

function [result] = dist_gte(d_min_sq, x_from, y_from, x_to, y_to)
%
% Return TRUE if minimum distance between the point (x_from, y_from) and 
% all points in (x_to, y_to) is greater than or equal to d_min_sq, else 
% return FALSE
%
% Arguments:
%   x_from, y_from: Scalars, single point to compute distance from
%   x_to, y_to: Vectors, many points to compute distance to
%   d_min: Scalar, mininum distance
% %

dxy = bsxfun(@minus, [x_from, y_from],  [x_to, y_to]);
result = min(sum(dxy.*dxy, 2)) >= d_min_sq;

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

function [xy, seed_idx] = get_streamline(xx, yy, uu, vv, x0, y0, step_size)
%
% Compute streamline in both directions starting at x0, y0
%
% Arguments: 
%   See documentation for stream2 for input argument definitions
%   xy : Matrix, stream line x- and y-coordinates in rows, returns [] if
%       stream line has zero length
%   seed_idx: Scalar, index of seed point in output streamline
% %

fwd = stream2(xx, yy, uu, vv, x0, y0, step_size);
xy_fwd = fwd{1}(~any(isnan(fwd{1}), 2), :); % drop NaN rows
has_fwd = size(xy_fwd,1) > 1;

rev = stream2(xx, yy, -uu, -vv, x0, y0, step_size);
xy_rev = rev{1}(~any(isnan(rev{1}), 2), :); % drop NaN rows
has_rev = size(xy_rev,1) > 1;

if has_fwd && has_rev
    xy = [xy_rev(end:-1:2, :); xy_fwd];
    seed_idx = size(xy_rev,1);
elseif has_rev
    xy = xy_rev;
    seed_idx = 1;
elseif has_fwd
    xy = xy_fwd;
    seed_idx = 1;
else
    xy = [];
    seed_idx = [];
end

% %<DEBUG>
% if any(isnan(xs)) || any(isnan(ys))
%     fprintf('NaN in streamline - debugging\n');
%     keyboard
% end
% %</DEBUG>

% %<DEBUG>
% if ~isempty(i0)
%     % i0 = i0+1; % fails, which is good
%     fprintf('x: %g\n', xs(i0)-x0);
%     assert(abs(xs(i0)-x0) < 1e-15, 'seed index is incorrect');
%     fprintf('y: %g\n', ys(i0)-y0);
%     assert(abs(ys(i0)-y0) < 1e-15, 'seed index is incorrect');
% end
% %</DEBUG>