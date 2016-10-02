%% Create example vector field

vv = linspace(-2, 2, 20);
hh = vv(2)-vv(1);
[xx, yy] = meshgrid(vv);
zz = xx .* exp(-xx.^2 - yy.^2);
[dzdx, dzdy] = gradient(zz, hh, hh);

%% Plot lines in all styles

% streamline parameters
low_dist_sep = 0.05*range(vv);
low_dist_test = 0.5*low_dist_sep;
high_dist_sep = 0.003*range(vv);
high_dist_test = 0.5*high_dist_sep;
step_size = 0.1;

% create figure
hf = figure;
hf.Name = mfilename;

% normal lines
xyld = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', low_dist_sep, 'DistTest', low_dist_test);

ax = subplot(2,2,1);
even_stream_line(xyld, 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '-');
title('even\_stream\_line');
ax.XTick = [];
ax.YTick = [];

% taper
xyld = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', low_dist_sep, 'DistTest', low_dist_test, 'GetDist', true);

ax = subplot(2,2,2);
even_stream_taper(xyld, 'LineWidthMin', 0.5, 'LineWidthMax', 2, 'Color', 'k');
title('even\_stream\_taper');
ax.XTick = [];
ax.YTick = [];

% arrow
xyld = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', low_dist_sep, 'DistTest', low_dist_test);

ax = subplot(2,2,3);
even_stream_arrow(xyld, 'LineStyle', '-', 'LineWidth', 0.5, 'Color', 'k', ...
    'ArrowLength', 2, 'ArrowTipAngle', 30, 'ArrowBaseAngle', 10, ...
    'ArrowSpace', 10);
title('even\_stream\_arrow');
ax.XTick = [];
ax.YTick = [];

% texture
xyld = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', high_dist_sep, 'DistTest', high_dist_test, 'GetLength', true);

ax = subplot(2,2,4);
even_stream_texture(xyld, 'LineWidth', 1, 'Period', 40);
title('even\_stream\_texture');
ax.XTick = [];
ax.YTick = [];