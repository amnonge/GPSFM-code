% function [trans, points] = normtrans(points)
%
%   Find a transformation to normalize some data points,
%   made of a translation so that the origin is the centroid,
%   and a rescaling so that the average distance to it is sqrt(2).
%
%   Input:
%     * points, the data points, if all the entries of the last row are ones, it is assumed it is homogeneous.
%     * isotropic, boolean indicating if the scaling should be isotropic (default) or not.
%   Output:
%     * trans, the transformation applied.
%     * points, the new data points coordinates.

% This code implements classical computer vision functions.
% It is provided for research purpose only.
% Copyright (C) 2017 Ludovic Magerand <ludovic@magerand.fr>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [trans, points] = normtrans(points, isotropic)

	if nargin < 2
		isotropic = true;
	end

 	[dimension, num_points] = size(points);
	
	if all(points(end,:) == 1)
		homogeneous = true;	
		dimension = dimension - 1;
	else
		homogeneous = false;
	end
	
	centroid = mean(points(1:dimension,:), 2);
	
	diff = points(1:dimension,:) - centroid(:, ones(num_points,1));
	
	if isotropic
		scale =	sqrt(2) ./ mean(sqrt(sum(diff.^2, 1)));

		trans = diag(vertcat(scale*ones(dimension,1), 1));
		trans(1:dimension, end) = -centroid*scale;
	else
		scale = sqrt(2) ./ mean(abs(diff),2);
		
		trans = diag(vertcat(scale, 1));
		trans(1:dimension, end) = -centroid.*scale;
	end
	
	if nargout > 1
		if homogeneous
			points = trans*points;
		else
			points = trans(1:dimension,1:dimension)*points + trans(1:dimension,end*ones(1,num_points));
		end
	end
	
end
