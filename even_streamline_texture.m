function hh = even_streamline_texture(xx, yy, uu, vv, d_sep, d_test, varargin)
% function hh = even_streamline_texture(xx, yy, uu, vv, d_sep, d_test, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] using the
% texturing effect to produce results similar to the line integral
% convolution (LIC) technique.
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
%   'LineWidth': line width, as in the plot(), default = 0.5
%   'Period': length of periodic pattern in # points, default = 20;
%   
% Returns:
%   hh = Vector of graphics objects for streamlines
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

parser.addParameter('StepSize', 0.1);
parser.addParameter('Verbose', false);
parser.addParameter('LineWidth', 0.5);
parser.addParameter('Period', 20);

parser.parse(varargin{:});
step_size = parser.Results.StepSize;
verbose = parser.Results.Verbose;
line_width = parser.Results.LineWidth;
period = parser.Results.Period;

validateattributes(period, {'numeric'}, {'scalar', 'positive', 'integer'}, ...
    mfilename, 'period');

% get streamline data
xy = get_stream_xy(xx, yy, uu, vv, d_sep, d_test, step_size, verbose);
len = get_stream_len(xy, verbose);

% get colors and color index for each segment
num_colors = 256;
idx = (1:length(len))';
color_coeff = 0.5*(1+sin(2*pi*idx/period)) + mod(idx, period)/(period-1);
color_coeff = (color_coeff-min(color_coeff))/range(color_coeff);
color_idx = min(num_colors, max(1, round(num_colors*color_coeff)));
colors = flipdud(gray(num_colors));

% create plot
num_segments = size(xy,1)-2*sum(isnan(xy(:,1)))-1;
hh = gobjects(num_segments,1);
current_segment = 0;
hold on;
for ii = 1:length(len)-1           
    % skip line endpoints, nothing to plot
    if any(isnan(xy(ii,:))) || any(isnan(xy(ii+1,:)))
        continue
    end
    % count and report
    current_segment = current_segment+1;
    if verbose 
        fprintf('%s: segment %d of %d\n', mfilename, current_segment, num_segments);
    end 
    % plot segment
    hh(current_segment) = plot(xy(ii:ii+1,1), xy(ii:ii+1,2), ...
        'LineWidth', line_width, 'Color', colors(color_idx(ii), :));
end