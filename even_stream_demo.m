%% Evenly Spaced Streamlines
% This package plots evenly spaced streamlines for a 2D vector in several
% styles described in Jobar & Lefer, 1997 [1]. Four types of plots are
% included: line, arrow, taper, texture. 
% 
% The heavy-lifting of computing evenly-spaced streamlines is handled by a
% modified version of streamslice() function, which uses an algorithm quite
% similar to [1]. The modifications provide more control over streamline
% spacing, and improved settings for high stream line densities.
%
% In the remainder of this document, each type of plot is genereated for an
% example vector field.

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

%% Simple stream lines
% Plot evenly spaced streamlines using a simple line style. This function
% is quite similar to the Mathworks function streamslice(), with added
% control over the maximum streamline density and better results for high
% streamline densities.

tic
hf = figure;
hf.Name = sprintf('%s: even_stream_line', mfilename);
even_stream_line(xx, yy, dzdx, dzdy, 2, 4, 'Color', 'k', 'LineWidth', 1);
title('even\_stream\_line');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_line: %.3f s elapsed\n', toc);

%% Plot stream lines with arrow glyphs 
% Compute streamlines and plot lines with arrow glyphs to indicate flow
% direction. 

tic
hf = figure;
hf.Name = sprintf('%s: even stream arrow', mfilename);
even_stream_arrow(xx, yy, dzdx, dzdy, 2, 4, ...
    'LineStyle', '-', 'LineWidth', 0.5, 'Color', 'k', 'ArrowLength', 4, ...
    'ArrowTipAngle', 30, 'ArrowBaseAngle', 10, 'ArrowSpace', 2);
title('even\_stream\_arrow');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_arrow: %.3f s elapsed\n', toc);

%% Plot tapered stream lines
% Plot stream lines with tapering effect, such that line width scales with
% the distance to the nearest neighboring streamline.

hf = figure;
hf.Name = sprintf('%s: even stream taper', mfilename);
even_stream_taper(xx, yy, dzdx, dzdy, 1, 2, ...
    'LineWidthMin', 0.5, 'LineWidthMax', 3, 'Color', 'k');
title('even\_stream\_taper');
ax = gca;
ax.XTick = [];
ax.YTick = [];
fprintf('even_stream_taper: %.3f s elapsed\n', toc);

%% Plot textured stream lines
% Plot streamlines with texture effect. This is most effective with
% closely-spaced streamlines, in which case it mimics the popular
% line-integral-convilution (LIC) method for visualizing flow fields.

hf = figure;
hf.Name = sprintf('%s: even stream texture', mfilename);
xy = even_stream_texture(xx, yy, dzdx, dzdy, 20, 40, ...
    'LineWidth', 1, 'Period', 1);
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
