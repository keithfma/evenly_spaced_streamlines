function [x_line, y_line] = ...
    even_stream(xx, yy, uu, vv, d_sep, d_test, step_size) %#ok!
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm (ref 1).
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
%   x_line, y_line: Vectors, x- and y-coordinates of streamline points,
%       individual lines are separated by NaNs
% 
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ’97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43–55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
% %

%% get initial streamline

% get seed point at random (populated) point
x_min = min(xx(:));
x_rng = range(xx(:)); 
y_min = min(yy(:));
y_rng = range(yy(:));
u0 = NaN;
v0 = NaN;
while isnan(u0) || isnan(v0)
    x0 = x_min+rand(1)*x_rng;
    y0 = y_min+rand(1)*y_rng;
    u0 = interp2(xx, yy, uu, x0, y0);
    v0 = interp2(xx, yy, vv, x0, y0);
end

% init stream line and seed candidates
[x_line, y_line] = get_streamline(xx, yy, uu, vv, x0, y0, step_size);
[x_seed, y_seed] = get_seed_candidates(x_line, y_line, d_sep);

%% main loop

% search current seed candidates for valid point
d_sep_sq = d_sep*d_sep;
for ii = randperm(length(x_seed));
    delta_x = x_seed(ii)-x_line;
    delta_y = y_seed(ii)-y_line;
    d_min_sq = min(delta_x.*delta_x+delta_y.*delta_y);
    if d_min_sq >= d_sep_sq
        x_next = x_seed(ii);
        y_next = y_seed(ii);
    end
end


%% debug plot

plot(x_line, y_line, '-k');
hold on
plot(x_seed, y_seed, '.r');
plot(x_next, y_next, 'ob');



function [x_seed, y_seed] = get_seed_candidates(x_line, y_line, d_sep)
%
% Compute the location of stream line seed point candidates that lie at a
% distance d_sep along a normal vector at each point in x_line, y_line 
%
% Arguments:
%   x_line, y_line: Vectors, points along a streamline
%   d_sep: Scalar, desired spacing between streamlines
%   x_seed, y_seed: Vectors, canditate seed points

% get unit normal vectors at segment midpoints
xy = [x_line, y_line];
tangent = diff(xy);
normal = [tangent(:,2), -tangent(:,1)];
normal = bsxfun(@rdivide, normal, sqrt(sum(normal.*normal, 2)));
midpoint = xy(1:end-1, :)+0.5*tangent;

% get candidates offset d_sep in positive and negative normal direction
seed = [midpoint+d_sep*normal; midpoint-d_sep*normal];
x_seed = seed(:,1);
y_seed = seed(:,2);

function [x_line, y_line] = get_streamline(xx, yy, uu, vv, x0, y0, step_size)
%
% Compute streamline in both directions starting at x0, y0
%
% Arguments: 
%   See documentation for stream2 for input argument definitions
%   x_line, y_line :
% %

xy = stream2(xx, yy, uu, vv, x0, y0, step_size);
x_fwd = xy{1}(:,1);
y_fwd = xy{1}(:,2);

xy = stream2(xx, yy, -uu, -vv, x0, y0, step_size);
x_rev = xy{1}(:,1);
y_rev = xy{1}(:,2);

x_line = [x_rev(end:-1:2); x_fwd];
y_line = [y_rev(end:-1:2); y_fwd];    