
% Run even_stream() for an example vector field and a range of input
% parameter options

%% Intialize example vector field
% see documentation for gradient()

vv = linspace(-2, 2, 20);
hh = vv(2)-vv(1);
[xx, yy] = meshgrid(vv);
zz = xx .* exp(-xx.^2 - yy.^2);
% %<DEBUG> Confirm program can deal with NaNs
% zz(xx>1.5) = NaN;
% %</DEBUG>
[dzdx, dzdy] = gradient(zz, hh, hh);

%% Plot evenly-spaced streamlines

d_sep = 0.05*range(vv);
d_test = 0.5*d_sep;
step_size = 0.1;
verbose = true;

xy = even_stream_xy(xx, yy, dzdx, dzdy, d_sep, d_test, step_size, verbose);
len = even_stream_len(xy, verbose);
dist = even_stream_dist(xy, verbose);

%<DEBUG>
hold off
% imagesc([xx(1), xx(end)], [yy(1), yy(end)], zz, 'AlphaData', ~isnan(zz));
% hold on
plot(xy(:,1), xy(:,2), '-k');
%</DEBUG>

