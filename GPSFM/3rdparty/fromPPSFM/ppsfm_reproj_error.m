% function [global_errors, reprojections, global_depths] = ppsfm_reproj_error(cameras, points, valid_cams, valid_pts, ...
%     normalisations, img_measurements, visible)
%
%   Compute reprojection errors and projective depths for the visible and valid points and cameras.
%
%   Input:
%     * cameras: Projective cameras estimation (3Fx4) 
%     * points: Projective points estimation (4xN)
%     * valid_cams: Indexes of cameras that has been estimated
%     * valid_pts: Indexes of points that has been estimated
%     * normalisations: Normalisation transformation for each camera stacked vertically (3Fx3)
%     * img_measurements: Unnormaliazed measurements of the projections (3FxN)
%     * visible: Visibility matrix binary mask (FxN)
%   Output:
%     * global_errors: Reprojections errors (FxN)
%     * reprojections: Reprojections (3FxN)

% This code implements the P2SfM method described in the ICCV'17 paper entitled
% "Practical Projective Structure from Motion (P2SfM)" by Ludovic Magerand and
% Alessio Del Bue.
% It is provided for research purpose only and appropriate citation to this
% paper should be made if you use it.
% The paper can be found at https://bitbucket.org/lmagerand/ppsfm/wiki/ or
% at the IEEE explore website (ieeexplore.ieee.org) without the supplementary
% materials.
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

function [global_errors, global_reprojs, global_depths] = ppsfm_reproj_error(cameras, points, valid_cams, valid_pts, ...
		normalisations, img_measurements, visible)

	valid_cams3 = kron(valid_cams, [3 3 3]) - kron(ones(1,length(valid_cams)), [2 1 0]);

	% Denormalize cameras
	denorm_cams = nan(length(valid_cams3), 4);
	for j = 1:3:length(valid_cams3);
		denorm_cams(j:j+2,:) = normalisations(valid_cams3(j):valid_cams3(j+2),:)...
			\ cameras(valid_cams3(j):valid_cams3(j+2),:);
	end


	% Projective depths and reprojection error
	scaled_measurements = denorm_cams * points(:,valid_pts);
	
	depths = scaled_measurements(3:3:end,:);
	
	reprojections = scaled_measurements ./ depths(kron(1:length(valid_cams), [1 1 1]), :);
	
	errors = reprojections .* visible(kron(valid_cams, [1 1 1]), valid_pts) - img_measurements(valid_cams3, valid_pts);
	errors = sqrt(errors(1:3:end,:).^2 + errors(2:3:end,:).^2);
	
	global_errors = nan(size(cameras,1)/3, size(points,2));
	global_depths = nan(size(cameras,1)/3, size(points,2));
	global_reprojs = nan(size(cameras,1), size(points,2));
	
	global_errors(valid_cams, valid_pts) = errors;
	global_depths(valid_cams, valid_pts) = depths;
	global_reprojs(valid_cams3, valid_pts) = reprojections;
end

