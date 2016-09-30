function hh = even_streamline_taper(xx, yy, uu, vv, d_sep, d_test, varargin)
% function hh = even_streamline_taper(xx, yy, uu, vv, d_sep, d_test, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] using the
% line tapering effect.
%
% Arguments:
%   xx, yy, uu, vv: Vector field x-coord, y-coord, vector x-component and
%       vector y-component, respectively, sizes must match
%   d_sep: Scalar, minimum distance between seed points and stream lines
%   d_test: Scalar, minimum distance between stream lines
%
% Optional Parameters (Name - Value):
%   'StepSize': stream line step size as in the built-in stream2
%   'Verbose': set true to enable verbose messages, default = false
%   'LineWidthMin': minimum line width as in the built-in plot(), default = 0.5
%   'LineWidthMax': maximum line width as in the built-in plot(), default = 2
%   'LineStyle': line style as in the built-in plot(), default = '-'
%   'Color': line color as in the built-in plot(), default = 'b'
%
% Returns:
%   hh = Graphics object for streamlines
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
parser.addParameter('LineWidthMin', 0.5);
parser.addParameter('LineWidthMax', 2);
parser.addParameter('LineStyle', '-');
parser.addParameter('Color', 'b');

parser.parse(varargin{:});
step_size = parser.Results.stepsize;
verbose = parser.Results.verbose;
line_width_min = parser.Results.LineWidthMin;
line_width_max = parser.Results.LineWidthMax;
line_style = parser.Results.LineStyle;
line_color = parser.Results.Color;

% create plot
xy = get_stream_xy(xx, yy, uu, vv, d_sep, d_test, step_size, verbose);
dist = get_stream_dist(xy, verbose);

%<DEBUG>
keyboard
%</DEBUG>