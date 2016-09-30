function dist = even_stream_dist(xy)
% function dist = even_stream_dist(xy)
%
% Return minimum distance to neighboring streamlines for streamlines
% created with the even_stream_xy fuction. This step is part of the Jobar &
% Lefer [1] algorithm to plot evenly spaced streamlines with tapered width.
%
% Arguments:
%   xy = Matrix, evenly-spaced streamlines, as created by the
%       even_stream_xy function, with [x,y] points in rows, and lines
%       separated by NaNs
%   dist = Vector, minimum distance to neighboring stream lines for all
%       points in xy
%
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
% %

error('%s is not implemented', mfilename);