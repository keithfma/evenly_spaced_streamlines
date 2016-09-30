function [] = even_streamline(xx, yy, uu, vv, d_sep, d_test, varargin)
% function [] = even_streamline(xx, yy, uu, vv, d_sep, d_test, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1]
%
% Arguments:
%   xx, yy, uu, vv: Vector field x-coord, y-coord, vector x-component and
%       vector y-component, respectively, sizes must match
%   d_sep: Scalar, minimum distance between seed points and stream lines
%   d_test: Scalar, minimum distance between stream lines
%
% Optional Parameters (Name - Value):
%   'StepSize': stream line step size as in the built-in stream2
%   'Verbose': set true to enable verbose progress messages
%   'LineStyle': line style as in the built-in plot()
%   'LineWidth': line width as in the built-in plot()
%   'Color': line color as in the built-in plot()
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

parser.parse(varargin{:});
step_size = parser.Results.stepsize;
verbose = parser.Results.verbose;
line_style = parser.Results.LineStyle;
line_width = parser.Results.LineWidth;
line_color = parser.Results.Color;

%<DEBUG>
disp(step_size);
disp(verbose);
disp(line_style);
disp(line_width);
disp(line_color);
%</DEBUG>
