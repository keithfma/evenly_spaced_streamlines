function [hh, xy] = pretty_stream_simple(varargin)
% Plot evenly-spaced streamlines in the style of Jobar & Lefer [1], using
% the built-in streamslice() to do the heavy lifting of computing lines.
%
% pretty_stream_simple(xx, yy, uu, vv)
% pretty_stream_simple(xx, yy, uu, vv, density)
% pretty_stream_simple(xy)
% pretty_stream_simple(..., Name, Value)
% [hh, xy] = pretty_stream_simple(...)
%
% Arguments:
%   xx, yy: Matrices or vectors, x-coord and y-coord. If matrices, the
%       size must match uu and vv. If vectors, xx must match the number of
%       columns in uu and vv, and yy must match the number of rows.
%   uu, vv: Matrices, vector field x-component and y-component
%   density: Scalar, coefficient specifying the density (spacing) of
%       streamlines. From the streamslice() documentation: "modifies the
%       automatic spacing of the streamlines. DENSITY must be greater than
%       0. The default value is 1; higher values will produce more
%       streamlines on each plane. For example, 2 will produce
%       approximately twice as many streamlines while 0.5 will produce
%       approximately half as many."
%
% Parameters (Name, Value):
%   'LineStyle': line style as in plot(), default = '-'
%   'LineWidth': line width as in plot(), default = 0.5
%   'Color': line color as in plot(), default = 'b'
%
% Returns:
%   hh = Graphics object for streamlines
%   xy: Matrix, [x, y] coordinates for stream line points, each row is a
%       point, individual lines are separated by NaNs.
%  
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
%
% Example: Plot and replot streamlines
%   [xx, yy] = meshgrid(0:0.2:2, 0:0.2:2);
%   uu = cos(xx).*yy;
%   vv = sin(xx).*yy;
%   subplot(1,2,1)
%   [~, xy] = pretty_stream_simple(xx, yy, uu, vv, 2, 'Color', 'b');
%   subplot(1,2,2)
%   pretty_stream_simple(xy, 'Color', 'r');
% %

% handle inputs
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = true;

parser.addRequired('xx_or_xy');
parser.addOptional('yy', []);
parser.addOptional('vv', []);
parser.addOptional('uu', []);
parser.addOptional('density', []);
parser.addParameter('LineStyle', '-');
parser.addParameter('LineWidth', 0.5);
parser.addParameter('Color', 'b');
parser.addParameter('Verbose', false);

parser.parse(varargin{:});
xx_or_xy = parser.Results.xx_or_xy;
yy = parser.Results.yy;
uu = parser.Results.uu;
vv = parser.Results.vv;
density = parser.Results.density;
line_style = parser.Results.LineStyle;
line_width = parser.Results.LineWidth;
line_color = parser.Results.Color;

% get evenly spaced streamlines
if isempty(yy) && isempty(uu) && isempty(vv) && isempty(density)
    % read xy from input arguments
    xy = xx_or_xy;
else
    % compute stream lines with built-in function and reformat as xy
    xx = xx_or_xy;
    [stream_cell, ~] = streamslice(xx, yy, uu, vv, density);    
    num_lines = length(stream_cell);
    stream_sep = cell(1, num_lines);
    [stream_sep{:}] = deal(nan(1,2));
    stream_cell = reshape([stream_cell; stream_sep], 2*num_lines, 1);
    xy = cell2mat(stream_cell);
end

% create plot
hh = plot(xy(:,1), xy(:,2), ...
    'LineStyle', line_style, 'LineWidth', line_width, 'Color', line_color);