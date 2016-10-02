function hh = even_stream_texture(xy, varargin)
% function hh = even_stream_texture(xy, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] using the
% texturing effect to produce results similar to the line integral
% convolution (LIC) technique.
%
% Arguments:
%   xy: Matrix with columns [x, y], containing streamline points as
%       produced by even_stream_data.
%
% Optional Parameters (Name - Value):
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
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

parser.addParameter('LineWidth', 0.5);
parser.addParameter('Period', 20);

parser.parse(varargin{:});
line_width = parser.Results.LineWidth;
period = parser.Results.Period;

validateattributes(period, {'numeric'}, {'scalar', 'positive', 'integer'}, ...
    mfilename, 'period');

% reformat streamlines as segments 
num_segments = size(xy, 1)-2*sum(isnan(xy(:,1)))-1;
x_segment = nan(num_segments, 2);
y_segment = nan(num_segments, 2);
current_segment = 0;
for ii = 1:size(xy,1)-1           
    % skip line endpoints, nothing to plot
    if any(isnan(xy(ii,1:2))) || any(isnan(xy(ii+1,1:2)))
        continue
    end
    current_segment = current_segment+1;
    x_segment(current_segment, :) = xy(ii:ii+1,1);
    y_segment(current_segment, :) = xy(ii:ii+1,2);
end

% get colormap
num_colors = 256;
colors = flipud(gray(num_colors));

% get colormap index for each segment
idx = (1:num_segments)';
color_coeff = 0.5*(1+sin(2*pi*idx/period)) + mod(idx, period)/(period-1);
color_coeff = (color_coeff-min(color_coeff))/range(color_coeff);
c_segment = min(num_colors, max(1, round(num_colors*color_coeff)));

% create plot: for efficieny, plot all segments with the same color at once
unique_colors = unique(c_segment);
num_unique_colors = length(unique_colors); 
hh = gobjects(num_unique_colors,1);
for ii = 1:num_unique_colors
    % get segments of this color
    cc = unique_colors(ii);
    is_cc = find(c_segment==cc);
    num_segment_cc = length(is_cc);
    x_segment_cc = nan(3*num_segment_cc, 1);
    x_segment_cc(1:3:end) = x_segment(is_cc, 1);
    x_segment_cc(2:3:end) = x_segment(is_cc, 2);    
    y_segment_cc = nan(3*num_segment_cc, 1);
    y_segment_cc(1:3:end) = y_segment(is_cc, 1);
    y_segment_cc(2:3:end) = y_segment(is_cc, 2);    
    % plot all segments at once
    hh(ii) = plot(x_segment_cc, y_segment_cc, ...
        'LineWidth', line_width, 'Color', colors(cc, :));
    hold on
end