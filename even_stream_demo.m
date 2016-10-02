%% Create example vector field

vv = linspace(-2, 2, 20);
hh = vv(2)-vv(1);
[xx, yy] = meshgrid(vv);
zz = xx .* exp(-xx.^2 - yy.^2);
[dzdx, dzdy] = gradient(zz, hh, hh);

%% Compute streamline data

dist_sep = 0.05*range(vv);
dist_test = 0.5*dist_sep;
step_size = 0.1;

xyld = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', dist_sep, 'DistTest', dist_test, 'StepSize', step_size, ...
    'GetLength', true, 'GetDist', true, 'Verbose', true);

%% Plot lines in all styles

hf = figure;
hf.Name = mfilename;

% normal lines
ax = subplot(2,2,1);
even_stream_line(xyld, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-');
title('even\_stream\_line');
ax.XTick = [];
ax.YTick = [];

% taper
ax = subplot(2,2,2);
even_stream_taper(xyld, 'LineWidthMin', 0.5, 'LineWidthMax', 3, 'Color', 'k');
title('even\_stream\_taper');
ax.XTick = [];
ax.YTick = [];

% arrow
ax = subplot(2,2,3);
title('even\_stream\_arrow');
ax.XTick = [];
ax.YTick = [];

% texture
ax = subplot(2,2,4);
title('even\_stream\_texture');
ax.XTick = [];
ax.YTick = [];
