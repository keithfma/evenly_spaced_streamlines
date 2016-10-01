function [hl, ha] = even_streamline_arrow(xx, yy, uu, vv, d_sep, d_test, varargin)
% function [hl, ha] = even_streamline_arrow(xx, yy, uu, vv, d_sep, d_test, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] with
% arrow glyphs to indicate the flow direction. Uses the 'arrow' package by
% Dr. Erik A. Johnson from the Mathworks FileExchange.
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
parser.addParameter('ArrowLength', 7);
parser.addParameter('ArrowTipAngle', 20);
parser.addParameter('ArrowBaseAngle', 45);
parser.addParameter('ArrowSpacing', 5);

parser.parse(varargin{:});
step_size = parser.Results.stepsize;
verbose = parser.Results.verbose;
line_style = parser.Results.LineStyle;
line_width = parser.Results.LineWidth;
color = parser.Results.Color;
arrow_length = parser.Results.ArrowLength;
arrow_tip_angle = parser.Results.ArrowTipAngle;
arrow_base_angle = parser.Results.ArrowBaseAngle;
arrow_spacing = parser.Results.ArrowSpacing;

% get streamline data
xy = get_stream_xy(xx, yy, uu, vv, d_sep, d_test, step_size, verbose);

% plot lines
hl = plot(xy(:,1), xy(:,2), ...
    'LineStyle', line_style, 'LineWidth', line_width, 'Color', color);

% plot arrows
xy_from = xy(arrow_spacing:arrow_spacing:end, :);
xy_to = xy(1+arrow_spacing:arrow_spacing:end, :);
xy_nan = any(isnan(xy_from), 2) | any(isnan(xy_to), 2);
xy_from = xy_from(~xy_nan, :);
xy_to = xy_to(~xy_nan, :);
axis(axis); % recommended in arrow() to fix axes limits 

ha = arrow(xy_from, xy_to, ...
    'Length', arrow_length, ...
    'TipAngle', arrow_tip_angle, ...
    'BaseAngle', arrow_base_angle, ...
    'FaceColor', color, 'LineStyle', 'none');
