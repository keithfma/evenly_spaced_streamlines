function [hl, ha] = even_streamline_arrow(xx, yy, uu, vv, d_sep, d_test, varargin)
% function [hl, ha] = even_streamline_arrow(xx, yy, uu, vv, d_sep, d_test, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] with
% arrow glyphs to indicate the flow direction.
%
% Arguments:
%   xx, yy, uu, vv: Vector field x-coord, y-coord, vector x-component and
%       vector y-component, respectively, sizes must match
%   d_sep: Scalar, minimum distance between seed points and stream lines
%   d_test: Scalar, minimum distance between stream lines
%
% Optional Parameters (Name - Value):
%   'StepSize': stream line step size as in the built-in stream2, default = 0.1
%   'Verbose': set true to enable verbose messages, default = false
%   'LineStyle': line style, as in the plot(), default = '-'
%   'LineWidth': line width, as in the plot(), default = 0.5
%   'Color': line and arrow color, as in plot(), default = 'b'
%   'HeadStyle': arrow head style, as in annotation(), default = 'vback2'
%   'HeadLength': arrow head length, as in annotation(), default = 10
%   'HeadWidth': arrow head width, as in annotation(), default = 10
%   'HeadSpc': arrow head spacing along line in points, default = 10
%   
% Returns:
%   hl = Graphics object for streamlines
%   ha = Vector of graphics objects for arrows
%  
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
% %

% handle inputs
% NOTE: sanity checks are defered to child functions
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

parser.addParameter('stepsize', 0.1);
parser.addParameter('verbose', false);
parser.addParameter('LineStyle', '-');
parser.addParameter('LineWidth', 0.5);
parser.addParameter('Color', 'b');
parser.addParameter('HeadStyle', 'vback2');
parser.addParameter('HeadLength', 10);
parser.addParameter('HeadWidth', 10);
parser.addParameter('HeadSpc', 10);

parser.parse(varargin{:});
step_size = parser.Results.stepsize;
verbose = parser.Results.verbose;
line_style = parser.Results.LineStyle;
line_width = parser.Results.LineWidth;
line_color = parser.Results.Color;
head_style = parser.Results.HeadStyle;
head_length = parser.Results.HeadLength;
head_width = parser.Results.HeadWidth;
head_spc= parser.Results.HeadSpc;

% get streamline data
xy = get_stream_xy(xx, yy, uu, vv, d_sep, d_test, step_size, verbose);

% plot lines
hl = plot(xy(:,1), xy(:,2), ...
    'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_color);

% plot arrows
xy_from = xy(head_spc:head_spc:end, :);
xy_to = xy(1+head_spc:head_spc:end, :);
xy_nan = any(isnan(xy_from), 2) | any(isnan(xy_to), 2);
xy_from = xy_from(~xy_nan, :);
xy_to = xy_to(~xy_nan, :);

arrow(xy_from, xy_to);
keyboard

% % plot arrows 
% % NOTE: uses the annotation() function to get decent arrowhead markers, but
% % this plots in normalized coordinates.
% x_arrow = [xy(head_spc:head_spc:end, 1), xy(1+head_spc:head_spc:end, 1)];
% y_arrow = [xy(head_spc:head_spc:end, 2), xy(1+head_spc:head_spc:end, 2)];
% nan_arrow = any(isnan(x_arrow), 2) | any(isnan(y_arrow), 2);
% x_arrow = x_arrow(~nan_arrow, :);
% y_arrow = y_arrow(~nan_arrow, :);
% keyboard
% ax = gca;
% norm_x_arrow = (x_arrow-ax.XLim(1))/(ax.XLim(2)-ax.XLim(1));
% norm_y_arrow = (y_arrow-ax.YLim(1))/(ax.YLim(2)-ax.YLim(1));
% ha = [];

