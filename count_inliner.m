function [num_inliners,error] = count_inliner(R,t,FeatureBag,featureExtracted,idx_3D,idx_2D,thresh)

num_inliners = 0;
error = zeros(length(idx_3D),1);
for i = 1:length(idx_3D)
    b = featureExtracted(:,idx_2D(i));
    b = b/norm(b);
    p = FeatureBag(1:3,idx_3D(i));
    b_est = R*p + t;
    b_est = b_est/norm(b_est);
    error(i) = norm(b_est-b);
    if error(i) < thresh
        num_inliners = num_inliners + 1;
    end
end
