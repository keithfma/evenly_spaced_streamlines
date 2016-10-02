%% Evenly Spaced Streamlines
% This package plots evenly spaced streamlines for a 2D vector field using
% the algorithm described in Jobar & Lefer, 1997 [1]. Four types of plots
% are included (as in [1]): line, arrow, taper, texture. Streamline
% computation and plotting are split, so that the data can be replotted
% without unnecessarily repeating the computations. In the remainder of
% this document, each type of plot is genereated for an example vector
% field.

%% Example vector field
% The vector field for is the gradient of a simple surface with a single
% local minimum and maximum.

% create vector field
vv = linspace(-2, 2, 20);
hh = vv(2)-vv(1);
[xx, yy] = meshgrid(vv);
zz = xx .* exp(-xx.^2 - yy.^2);
[dzdx, dzdy] = gradient(zz, hh, hh);

% plot surface and vector field
hf = figure;
hf.Name = sprintf('%s: example vector field', mfilename);
imagesc([vv(1), vv(end)], [vv(1), vv(end)], zz);
hold on
quiver(xx, yy, dzdx, dzdy, 'Color', 'k');
title('Example Vector Field');
ax = gca;
ax.XTick = [];
ax.YTick = [];

%% Plot stream lines
% Compute streamlines and plot using a simple line style.

% parameters
dist_sep = 0.05*range(vv);
dist_test = 0.5*dist_sep;

% compute
tic
xy = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', dist_sep, 'DistTest', dist_test);
fprintf('even_stream_data: %.3f s elapsed\n', toc);

% plot
tic
hf = figure;
hf.Name = sprintf('%s: even stream line', mfilename);
even_stream_line(xy, 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-');
title('even\_stream\_line');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_line: %.3f s elapsed\n', toc);

%% Plot stream lines with arrow glyphs 
% Compute streamlines and plot lines with arrow glyphs to indicate flow
% direction. 

% parameters
dist_sep = 0.05*range(vv);
dist_test = 0.5*dist_sep;

% compute
tic
xy = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', dist_sep, 'DistTest', dist_test);
fprintf('even_stream_data: %.3f s elapsed\n', toc);

% plot
tic
hf = figure;
hf.Name = sprintf('%s: even stream arrow', mfilename);
even_stream_arrow(xy, 'LineStyle', '-', 'LineWidth', 0.5, 'Color', 'k', ...
    'ArrowLength', 5, 'ArrowTipAngle', 30, 'ArrowBaseAngle', 10, ...
    'ArrowSpace', 10);
title('even\_stream\_arrow');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_arrow: %.3f s elapsed\n', toc);

%% Plot tapered stream lines
% Plot stream lines with tapering effect, such that line width scales with
% the distance to the nearest neighboring streamline.

% parameters
dist_sep = 0.05*range(vv);
dist_test = 0.5*dist_sep;

% compute
tic
[xy, dist] = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', dist_sep, 'DistTest', dist_test);
fprintf('even_stream_data: %.3f s elapsed\n', toc);

% plot
tic
hf = figure;
hf.Name = sprintf('%s: even stream taper', mfilename);
even_stream_taper(xy, dist, 'LineWidthMin', 0.5, 'LineWidthMax', 2, 'Color', 'k');
title('even\_stream\_taper');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_taper: %.3f s elapsed\n', toc);

%% Plot textured stream lines
% Plot streamlines with texture effect. This is most effective with
% closely-spaced streamlines, in which case it mimics the popular
% line-integral-convilution (LIC) method for visualizing flow fields.

% parameters
dist_sep = 0.003*range(vv);
dist_test = 0.5*dist_sep;

% compute
tic
xy = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', dist_sep, 'DistTest', dist_test);
fprintf('even_stream_data: %.3f s elapsed\n', toc);

% plot
tic
hf = figure;
hf.Name = sprintf('%s: even stream texture', mfilename);
even_stream_texture(xy, 'LineWidth', 1, 'Period', 20);
title('even\_stream\_texture');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_texture: %.3f s elapsed\n', toc);


%% References
% # Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
