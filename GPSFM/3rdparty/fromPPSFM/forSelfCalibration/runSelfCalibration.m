function [Psb_m,Xs_m,R_output1,T_output1]=runSelfCalibration( Psb,Xs,normalization,centers)
Ps=zeros(length(Psb)*3,4);
normmMat=zeros(length(Psb)*3,3);
for i=1:length(Psb)
   
normmMat(3*i-2:3*i,:)= inv(normalization);%normalization(3*i-2:3*i,3*i-2:3*i);%
end

for i=1:length(Psb)
Psb{i}=Psb{i}/norm(Psb{i},'fro');
Psb{i}=Psb{i}/sign(det(Psb{i}(1:3,1:3)));
    Ps(3*i-2:3*i,:)=normmMat(3*i-2:3*i,:)*Psb{i};
Ps(3*i-2:3*i,:)=Ps(3*i-2:3*i,:)/norm(Ps(3*i-2:3*i,:),'fro');

end
options=struct;
options.sedumi_path='';
options.gloptipoly_path='';
options.debug=1;
options.verbose=1;
options.upgrade_threshold=15;
[metric_cameras, metric_points, metric_upgrade] = metr(Ps, Xs, normmMat, centers', options);


Psb_m=metric_cameras;
Xs_m=metric_points;
figure, pcshow(metric_points(1:3,:)');
hold on;
axis off
R_output1=cell(length(Psb),1);
T_output1=cell(length(Psb),1);
[~,~,t1]=vgg_KR_from_P(Psb_m(3*1-2:3*1,:));
minDist=inf;
for i=1:length(Psb)
    [K,RT,t]=vgg_KR_from_P(Psb_m(3*i-2:3*i,:));
    if i>1
        curdist=norm(t-t1);
        if curdist<minDist
            minDist=curdist;
        end
    end
    RT=RT/det(RT(1:3,1:3));
    R_output1{i}=RT(1:3,1:3)';
    
    T_output1{i}=t;
    
end
scaleCam=mean(std(metric_points(1:3,:)'));
for i=1:length(Psb)
    plotCamera('Location',T_output1{i},'Orientation', R_output1{i}','Opacity',0,'Size',scaleCam/3,'Color',[0 1 0 ]);hold on;
end

end
 

function [metric_cameras, metric_points, metric_upgrade] = metr(cameras, points, normalisations, centers, options)

	addpath(options.sedumi_path);
	addpath(options.gloptipoly_path);
	
	num_cams = size(cameras,1)/3;
	
	% Cleaning gloptipoly global var
	mset('clear');
	mset('verbose', options.debug);

	% Move the principal point to the origin
	for j = 1:num_cams
		norm_center = normalisations(3*j-2:3*j,:)*vertcat(centers(:,j),1);
		cameras(3*j-2:3*j,:) = [1 0 -norm_center(1); 0 1 -norm_center(2); 0 0 1] * cameras(3*j-2:3*j,:);
	end


	% Unknowns
	mpol('O',10,1);
    

	Om = [O(1) O(2) O(3) O(4);
		  O(2) O(5) O(6) O(7);
		  O(3) O(6) O(8) O(9);
		  O(4) O(7) O(9) O(10)]; % The dual absolute quadric


	% Constraints

	% Fixing the scale of the dual absolute quadric
	scale_daq = supcon(sum(Om(:).^2), 1, 'eq');

	% Rank deficiency
	rank_def = supcon(det(Om), 0, 'eq');

	% Positive SemiDefiniteness
	psd = supcon();
	psd = psd(ones(1,14));
	psd(1) = supcon(Om(1,1), 0, 'ge');
	psd(2) = supcon(Om(2,2), 0, 'ge');
	psd(3) = supcon(Om(3,3), 0, 'ge');
	psd(4) = supcon(Om(4,4), 0, 'ge');
	psd(5) = supcon(det(Om([1 2],[1 2])), 0, 'ge');
	psd(6) = supcon(det(Om([1 3],[1 3])), 0, 'ge');
	psd(7) = supcon(det(Om([1 4],[1 4])), 0, 'ge');
	psd(8) = supcon(det(Om([2 3],[2 3])), 0, 'ge');
	psd(9) = supcon(det(Om([2 4],[2 4])), 0, 'ge');
	psd(10) = supcon(det(Om([3 4],[3 4])), 0, 'ge');
	psd(11) = supcon(det(Om([2 3 4],[2 3 4])), 0, 'ge');
	psd(12) = supcon(det(Om([1 3 4],[1 3 4])), 0, 'ge');
	psd(13) = supcon(det(Om([1 2 4],[1 2 4])), 0, 'ge');
	psd(14) = supcon(det(Om([1 2 3],[1 2 3])), 0, 'ge');




	% Objective function
	obj = mpol();
	for j = 1:num_cams
		ppx = (cameras(3*j-2,:)*Om*cameras(3*j,:)')^2;
		ppy = (cameras(3*j-1,:)*Om*cameras(3*j,:)')^2;
		skew = (cameras(3*j-2,:)*Om*cameras(3*j-1,:)')^2;
		ratio = (cameras(3*j-2,:)*Om*cameras(3*j-2,:)' - cameras(3*j-1,:)*Om*cameras(3*j-1,:)')^2;
		obj = obj + (skew + ppx + ppy + ratio)/4;
	end
	obj = momcon(obj/num_cams, 'min');


	% Solving

	prob = msdp(obj, psd, rank_def, scale_daq); % chir, 
	status = msol(prob);
	
	if status == 0
		if options.verbose
			fprintf('  Increasing test tolerances (results might be innacurate)\n');
		end
 		mset('testol', options.upgrade_threshold);
		status = msol(prob);
	end
	
	if status == 1
		% Extract solution(s)
		MQ = double(Om);
		if size(MQ,3) > 1
			if options.verbose
				fprintf('  Multiple solutions found (%d), taking first one\n', size(MQ,3));
			end
			MQ = squeeze(MQ(:,:,1));
		end
		
		% Enforce the correct rank exactly
		[U,S,V] = svd(MQ);
		S(4,4) = 0;
		MQ = U*S*V';

		% Enforce exact symetry
		MQ(2,1) = (MQ(1,2)+MQ(2,1))/2; MQ(1,2) = MQ(2,1);
		MQ(3,1) = (MQ(1,3)+MQ(3,1))/2; MQ(1,3) = MQ(3,1);
		MQ(4,1) = (MQ(1,4)+MQ(4,1))/2; MQ(1,4) = MQ(4,1);
		MQ(3,2) = (MQ(2,3)+MQ(3,2))/2; MQ(2,3) = MQ(3,2);
		MQ(4,2) = (MQ(2,4)+MQ(4,2))/2; MQ(2,4) = MQ(4,2);
		MQ(4,3) = (MQ(3,4)+MQ(4,3))/2; MQ(3,4) = MQ(4,3);

		% Computing the ambiguity
		[U,S] = svd(MQ);
		metric_upgrade = U*sqrt(S);
		metric_upgrade(:,4) = [0;0;0;1];
% 		metric_upgrade(4,4) = 1;

		% Upgrade the reconstruction
		metric_cameras = cameras*metric_upgrade;
		for j = 1:num_cams
			KKT = metric_cameras(3*j-2:3*j,1:3)*metric_cameras(3*j-2:3*j,1:3)';
			metric_cameras(3*j-2:3*j,:) = metric_cameras(3*j-2:3*j,:) / sqrt(KKT(end,end));
		end
		
		metric_points = metric_upgrade \ points;
		metric_points = metric_points ./ metric_points([4 4 4 4], :);
		
	else
		warning('Estimation of the dual absolute quadric failed (status = %d)\n', status);
		metric_cameras = [];
		metric_points = [];
		metric_upgrade = [];
	end

end
