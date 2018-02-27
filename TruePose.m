function [R_true,t_true,idx]=TruePose(R,t,p,b,K)

%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%
% p: 3D features, w.r.t global frame
% b: 2D features, 3X1 vector, w.r.t camera frame
% K: intrinsic matrix


b = b/norm(b);
for i=1:size(R,3)
    b_est = R(:,:,i)*p + t(:,i);
    b_est = b_est/norm(b_est);
    error(i) = norm(b_est-b);
end
[~,idx] = min(error);
R_true = R(:,:,idx);
t_true = t(:,idx);