clear all; close all;
addpath('../');
addpath('../robotics3D/');


NumOfFeatures = 4;
NumOfPoses = 2;

PoseGraphMatrix = ones(NumOfFeatures, NumOfPoses);

maxIter = 30;
whole_iter = 1;
% for whole_iter = 1:maxIter

%% 3D-WORLD creation


% To ensure f.o.v = 60 deg.
fov = 60*pi/180;
Rfeat = 1; Rcam = Rfeat/sin(fov/2);


% Random pose generator
W_T_C1 = [];
AbsolutePoses_true = zeros(3,4,NumOfPoses);
for k = 1:NumOfPoses
    tt = 2*pi*rand();
    WRC = [-sin(tt) 0 -cos(tt); cos(tt) 0 -sin(tt); 0 -1 0];
    WpC = Rcam*[cos(tt);sin(tt);0];
    if k == 1
        W_T_C1 = [WRC WpC];
        AbsolutePoses_true(:,:,k) = [eye(3) zeros(3,1)];
    else
        C_T_W = [WRC' -WRC'*WpC];
        AbsolutePoses_true(:,:,k) = C_T_W*[W_T_C1;zeros(1,3) 1];
    end
end


% Random feature generator
FeatureBag_true = zeros(4,NumOfFeatures);
%figure(2); hold on;
for k = 1:NumOfFeatures
    tt = 2*pi*rand();
    phi = pi*rand() - pi/2;
    FeatureBag_true(:,k) = [Rfeat*[cos(phi)*cos(tt);cos(phi)*sin(tt);sin(phi)];1]; % <--- w.r.t the world
    FeatureBag_true(1:3,k) = InversePose(W_T_C1)*[FeatureBag_true(1:3,k);1];
end
figure(1);
set(gcf,'visible','off');
plot3(FeatureBag_true(1,:),FeatureBag_true(2,:),FeatureBag_true(3,:),'*');title('3D points');
% saveas(gcf,['/home/jahid/Fei/summer/5Point_new/test/iter',num2str(whole_iter),' 3D points.jpg']);

% Create feature projections as measurement
featureExtracted_true = cell(NumOfPoses, 1);
for k = 1:NumOfPoses    
    numFeatPosek = find(PoseGraphMatrix(:,k));
    featureExtracted_true{k} = zeros(2, length(numFeatPosek));
    for l = 1:length(numFeatPosek)
        % assign correct feature index for pose graph matrix
        PoseGraphMatrix(numFeatPosek(l),k) = l; 
        
        % Pose Transformation
        Ck_T_W = AbsolutePoses_true(:,:,k);
        Cl_T_W = AbsolutePoses_true(:,:,FeatureBag_true(4,numFeatPosek(l)));
        Ck_T_Cl = Ck_T_W*[InversePose(Cl_T_W);zeros(1,3) 1];
        
        % Image Projection
        Ck_p_fl = Ck_T_Cl*[FeatureBag_true(1:3,numFeatPosek(l));1];
        featureExtracted_true{k}(:,l) = Ck_p_fl(1:2)./Ck_p_fl(3);
    end
    featureExtracted_true{k} = [featureExtracted_true{k};ones(1,size(featureExtracted_true{k},2))];
    figure(k+1);
    set(gcf,'visible','off');
    plot(featureExtracted_true{k}(1,:),featureExtracted_true{k}(2,:),'*');title(['2D points in camera',num2str(k)]);
%     saveas(gcf,['/home/jahid/Fei/summer/5Point_new/test/iter',num2str(whole_iter),' 2D points camera',num2str(k),'.jpg']);
end

% [R_true,t_true] = P3P_Solver(featureExtracted_true{1}(:,1:4),featureExtracted_true{2}(:,1:4),eye(3));
[C,p] = P3P_Solver_2(FeatureBag_true(1:3,1:3),featureExtracted_true{2}(:,1:3));
[R_est,t_est,idx]=TruePose(C,p,FeatureBag_true(1:3,4),featureExtracted_true{2}(:,4),eye(3));


%% display
figure;
% set(gcf,'visible','off');
Draw3DWorld([R_est,t_est], FeatureBag_true, W_T_C1);hold on
figure;
Draw3DWorld(AbsolutePoses_true(:,:,2), FeatureBag_true, W_T_C1);
% saveas(gcf,['/home/jahid/Fei/summer/P3P/test/iter',num2str(whole_iter),' camera pose.jpg']);
% close all
disp(['Iter: ',num2str(whole_iter)]);
% end





