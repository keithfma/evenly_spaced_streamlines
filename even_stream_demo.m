%% Create example vector field

vv = linspace(-2, 2, 20);
hh = vv(2)-vv(1);
[xx, yy] = meshgrid(vv);
zz = xx .* exp(-xx.^2 - yy.^2);
[dzdx, dzdy] = gradient(zz, hh, hh);

%% Compute streamline data

d_sep = 0.05*range(vv);
d_test = 0.5*d_sep;
step_size = 0.1;

xyld = even_stream_data(xx, yy, dzdx, dzdy, ...
    'DistSep', d_sep, 'DistTest', d_test, 'StepSize', step_size, ...
    'GetLength', false, 'GetDist', true, 'Verbose', true);

%% Plot lines in all styles

% hold off
% 
% % hh = even_stream_demo(xx, yy, dzdx, dzdy, d_sep, d_test, ...
% %     'Color', 'r', 'LineWidth', 2, 'Verbose', 1);
% 
% hh = even_stream_taper(xx, yy, dzdx, dzdy, d_sep, d_test, ...
%     'Color', 'r', 'LineWidthMin', 0.5, 'LineWidthMax', 5, 'Verbose', 1);
% 
% % [hl, ha] = even_stream_arrow(xx, yy, dzdx, dzdy, d_sep, d_test, ...
% %     'Color', 'r', 'LineWidth', 1, 'Verbose', 1);
% 
% % hh = even_stream_texture(xx, yy, dzdx, dzdy, d_sep, d_test, ...
% %     'Verbose', 1, 'LineWidth', 3);

