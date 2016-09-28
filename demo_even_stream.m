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
step_size = 0.1*d_sep;

xyld = even_stream2(xx, yy, dzdx, dzdy, d_sep, d_test, step_size);

%<DEBUG>
figure
imagesc([xx(1), xx(end)], [yy(1), yy(end)], zz, 'AlphaData', ~isnan(zz));
hold on
plot(xyld(:,1), xyld(:,2), '-k');
%</DEBUG>

