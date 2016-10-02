function [hl, ha] = even_stream_arrow(xyd, varargin)
% function [hl, ha] = even_stream_arrow(xyd, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] with
% arrow glyphs to indicate the flow direction. Uses the 'arrow' package by
% Dr. Erik A. Johnson from the Mathworks File Exchange, which is included
% below.
%
% Arguments:
%   xyd: Matrix with columns [x, y, dist], as produced by
%       even_stream_data. Only the x and y columns are needed.
%
% Optional Parameters (Name - Value):
%   'LineStyle': line style, as in plot(), default = '-'
%   'LineWidth': line width, as in plot(), default = 0.5
%   'Color': line and arrow color, as in plot(), default = 'b'
%   'ArrowLength': arrow head length in pixels, default = 20
%   'ArrowTipAngle': arrow head tip angle in degrees, default = 20
%   'ArrowBaseAngle': arrow head base angle in degrees, default = 10
%   'ArrowSpace': arrow head spacing along line in # points, default = 10
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
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

parser.addParameter('LineStyle', '-');
parser.addParameter('LineWidth', 0.5);
parser.addParameter('Color', 'b');
parser.addParameter('ArrowLength', 7);
parser.addParameter('ArrowTipAngle', 20);
parser.addParameter('ArrowBaseAngle', 45);
parser.addParameter('ArrowSpace', 5);

parser.parse(varargin{:});
line_style = parser.Results.LineStyle;
line_width = parser.Results.LineWidth;
color = parser.Results.Color;
arrow_length = parser.Results.ArrowLength;
arrow_tip_angle = parser.Results.ArrowTipAngle;
arrow_base_angle = parser.Results.ArrowBaseAngle;
arrow_space = parser.Results.ArrowSpace;

% plot lines
hl = plot(xyd(:,1), xyd(:,2), ...
    'LineStyle', line_style, 'LineWidth', line_width, 'Color', color);

% plot arrows
xy_from = xyd(arrow_space:arrow_space:end, 1:2);
xy_to = xyd(1+arrow_space:arrow_space:end, 1:2);
xy_nan = any(isnan(xy_from), 2) | any(isnan(xy_to), 2);
xy_from = xy_from(~xy_nan, :);
xy_to = xy_to(~xy_nan, :);
axis(axis); % recommended in arrow() to fix axes limits 

ha = arrow(xy_from, xy_to, ...
    'Length', arrow_length, ...
    'TipAngle', arrow_tip_angle, ...
    'BaseAngle', arrow_base_angle, ...
    'FaceColor', color, 'LineStyle', 'none');