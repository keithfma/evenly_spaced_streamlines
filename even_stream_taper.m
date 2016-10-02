function hh = even_stream_taper(xyld, varargin)
% function hh = even_stream_taper(xyld, varargin)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1] using the
% line tapering effect.
%
% Arguments:
%   xyld: Matrix with columns [x, y, len, dist], as produced by
%       even_stream_data. Only the x, y, and len columns are needed.
%
% Optional Parameters (Name - Value):
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
parser = inputParser;
parser.CaseSensitive = false;
parser.PartialMatching = false;
parser.KeepUnmatched = false;

parser.addParameter('LineWidthMin', 0.5);
parser.addParameter('LineWidthMax', 2);
parser.addParameter('Color', 'b');

parser.parse(varargin{:});
line_width_min = parser.Results.LineWidthMin;
line_width_max = parser.Results.LineWidthMax;
line_color = parser.Results.Color;

% reformat streamlines as segments 
num_segments = size(xyld,1)-2*sum(isnan(xyld(:,1)))-1;
x_segment = nan(num_segments, 2);
y_segment = nan(num_segments, 2);
d_segment = nan(num_segments, 1);
current_segment = 0;
for ii = 1:size(xyld,1)-1           
    % skip line endpoints, nothing to plot
    if any(isnan(xyld(ii,1:2))) || any(isnan(xyld(ii+1,1:2)))
        continue
    end
    current_segment = current_segment+1;
    x_segment(current_segment, :) = xyld(ii:ii+1,1);
    y_segment(current_segment, :) = xyld(ii:ii+1,2);
    d_segment(current_segment) = mean(xyld(ii:ii+1,4));
end

% get segment width rounded to 0.1 pt
dist_max = max(xyld(:,4));
dist_min = min(xyld(:,4));
w_coef = max(0.001, min(1, (d_segment-dist_min)/(dist_max-dist_min)));
w_segment = line_width_min+w_coef*(line_width_max-line_width_min);
w_segment = round(10*w_segment)/10;

% create plot: for efficieny, plot all segments with the same width at once
unique_widths = unique(w_segment);
num_unique_widths = length(unique_widths); 
hh = gobjects(num_unique_widths,1);
for ii = 1:num_unique_widths
    % get segments of this width
    ww = unique_widths(ii);
    is_ww = find(w_segment==ww);
    num_segment_ww = length(is_ww);
    x_segment_ww = nan(3*num_segment_ww, 1);
    x_segment_ww(1:3:end) = x_segment(is_ww, 1);
    x_segment_ww(2:3:end) = x_segment(is_ww, 2);    
    y_segment_ww = nan(3*num_segment_ww, 1);
    y_segment_ww(1:3:end) = y_segment(is_ww, 1);
    y_segment_ww(2:3:end) = y_segment(is_ww, 2);    
    % plot all segments at once
    hh(ii) = plot(x_segment_ww, y_segment_ww, ... 
        'LineWidth', ww, 'Color', line_color);   
    hold on
end
