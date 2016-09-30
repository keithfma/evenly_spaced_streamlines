function [] = even_streamline(xx, yy, uu, vv, d_sep, d_test, step_size, verbose)
% function [] = even_streamline(xx, yy, uu, vv, d_sep, d_test, step_size, verbose)
%
% Plot evenly-spaced streamlines with Jobar & Lefer algorithm [1]
%
% Arguments:
%   xx, yy, uu, vv: Vector field x-coord, y-coord, vector x-component and
%       vector y-component, respectively, sizes must match
%   d_sep: Scalar, minimum distance between seed points and stream lines
%   d_test: Scalar, minimum distance between stream lines
%   step_size: Scalar, stream line step size as in the built-in stream2
%   verbose: Scalar, set to True to enable verbose progress messages
%  
% References: 
% [1] Jobard, B., & Lefer, W. (1997). Creating Evenly-Spaced Streamlines of
%   Arbitrary Density. In W. Lefer & M. Grave (Eds.), Visualization in
%   Scientific Computing ?97: Proceedings of the Eurographics Workshop in
%   Boulogne-sur-Mer France, April 28--30, 1997 (pp. 43?55). inbook,
%   Vienna: Springer Vienna. http://doi.org/10.1007/978-3-7091-6876-9_5
% %

error('%s not implemented', mfilename);