%% Create example vector field

vv = linspace(-2, 2, 20);
hh = vv(2)-vv(1);
[xx, yy] = meshgrid(vv);
zz = xx .* exp(-xx.^2 - yy.^2);
[dzdx, dzdy] = gradient(zz, hh, hh);
d_sep = 0.05*range(vv);
d_test = 0.5*d_sep;
step_size = 0.1;

% %<DEBUG> Compute streamline data
% verbose = true;
% xy = get_stream_xy(xx, yy, dzdx, dzdy, d_sep, d_test, step_size, verbose);
% len = get_stream_len(xy, verbose);
% dist = get_stream_dist(xy, verbose);
% %</DEBUG>

%% Example plots

%<DEBUG>
hold off
%</DEBUG>

% hh = even_streamline(xx, yy, dzdx, dzdy, d_sep, d_test, ...
%     'Color', 'r', 'LineWidth', 2, 'Verbose', 1);

% hh = even_streamline_taper(xx, yy, dzdx, dzdy, d_sep, d_test, ...
%     'Color', 'r', 'LineWidthMin', 0.5, 'LineWidthMax', 5, 'Verbose', 1);

[hl, ha] = even_streamline_arrow(xx, yy, dzdx, dzdy, d_sep, d_test, ...
    'Color', 'r', 'LineWidth', 1, 'Verbose', 1);

%<DEBUG>
%</DEBUG>

