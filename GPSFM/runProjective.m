function [ repError,elapsedTime,sizes ] = runProjective( dataset,selfcalib)
%RUNPROJECTIVE Summary of this function goes here
%   Detailed explanation goes here


load(sprintf('../DataSet Proj/%s.mat', dataset));

% Generating triplet graph
originaltic=tic;
[finalTriplets,v] = buildTripletsAviodCollinearFast(pointMatchesInliers,FN);


% Normalize the measurment matrix
normMat = getNormalizationMatBigCondition(1,M);
FN = normMat'*FN*normMat;
[ FN ] = normalizeForbineusNorm( FN );
FN(isnan(FN)) = 0;


% Perform the triplets optimization
[Xs,Y] = optimizeTripletsIRLS( FN,finalTriplets(1:end,1:3),true,false ); 
nodesNum = length(unique(finalTriplets));

% Extracting camera matrices from a collection of consistent triplets
Pss=getProjectiveCameraMatrices(Xs,v,finalTriplets);
[Fnew,Anew]=getBigFfromCameras(Pss,finalTriplets,nodesNum);
Fnew = inv(normMat')*Fnew*inv(normMat);
[ ~,~,~,Ps,~ ] = getCameraMatrices(Fnew);


% Taking subset of points for bundle adjustment
[newM,~,ignoredPoints] = dilutePoint(M);


% Running two rounds of bundle adjustment
[ Psb , ~,~ ] = improveSolutionBACeres( Ps,newM );
[ Psb , Points3DBA ] = improveSolutionBACeres( Psb,newM );


elapsedTime = toc(originaltic);  
% Recostructing the points we ommited before the BA (not included in
% timing)
[Points3Dtemp,~,~] = get3Dpoint(Psb,M(:,ignoredPoints),M);
Points3D = zeros(4,length(ignoredPoints));
Points3D(:,~ignoredPoints) = Points3DBA;
Points3D(:,ignoredPoints) = Points3Dtemp;
[ PsRep ,visible,MRep] = prepareDataForRep(Psb,0,M);


% Compute reporjection error and RMS
 [global_errors, global_reprojs, global_depths] = ppsfm_reproj_error(PsRep, Points3D, 1:length(PsRep)/3, 1:length(Points3D), ...
		repmat(eye(3),length(PsRep),1),MRep,visible);
 display(['Reprojection Error is' ] )
repError = nanmean(nanmean(global_errors(global_errors~=0)));
% RMS = sqrt(nanmean(nanmean(global_errors(global_errors~=0).^2)));


%Print the final results
display('reprojection error| RMS|TotalTime')
[repError,elapsedTime]
%number of cameras/ number of points
sizes=[size(M,1)/2,size(M,2)]



if selfcalib
[Psb_m,Xs_m,R_output1_N,T_output1_N]= runSelfCalibration( Psb,Points3D,[1000 0 500;0 1000 500;0 0 1],repmat( [288.0000  360.0000],size(FN,1)/3,1));
colors=Xs_m(1:3,:)';
for i=1:3
colors(:,i)=colors(:,i)-min(colors(:,i));
colors(:,i)=colors(:,i)/max(colors(:,i));
end

if ~exist('reconstructions', 'dir')
	mkdir('reconstructions')
end
pcwrite(pointCloud(Xs_m(1:3,:)','Color',uint8(255*colors)),sprintf('reconstructions/%s.ply',dataset));
savefig(sprintf('reconstructions/%s.fig',dataset));

end
end

