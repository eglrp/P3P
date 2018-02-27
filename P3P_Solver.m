function [R_true,t_true] = P3P_Solver(P,b,K)
%%%%%%%% input %%%%%%%%
% Pi: 3D feature refer to global frame
% bi: 2D features refer to camera frame

%%%%%%%% output %%%%%%%%
% R : rotation
% t : tanslation

%%%%%%%% Based on Tong Ke's paper %%%%%%%%%%%%

p1 = P(:,1);
p2 = P(:,2);
p3 = P(:,3);
p4 = P(:,4);
b1 = b(:,1);
b2 = b(:,2);
b3 = b(:,3);
b4 = b(:,4);

% nomalize bi to be unit vector
b1 = b1/norm(b1);
b2 = b2/norm(b2);
b3 = b3/norm(b3);


% C = C(k1,theta1)C(k2,theta2)C(k3,theta3)
% equation(6) 
k1 = (p1-p2)/norm(p1-p2); %%%%% ask about the sign of k1!!!!!!
k3 = cross(b1,b2);k3 = k3/norm(k3);

% equation(11) 
u1 = p1 - p3;
u2 = p2 - p3;
v1 = cross(b1,b3);
v2 = cross(b2,b3);

% equation(22)(23)
sigma = norm(cross(u1,k1));
k3_prim_prim = cross(u1,k1)/sigma;

% equation(29)-(36)
f11 = sigma*transpose(k3)*b3;
f21 = sigma*(transpose(b1)*b2)*(transpose(k3)*b3);
f22 = sigma*(transpose(k3)*b3)*norm(cross(b1,b2));
f13 = sigma*transpose(v1)*k3;
f23 = sigma*transpose(v2)*k3;
f24 = (transpose(u2)*k1)*(transpose(k3)*b3)*norm(cross(b1,b2));
f15 = -(transpose(u1)*k1)*(transpose(k3)*b3);
f25 = -(transpose(u2)*k1)*(transpose(b1)*b2)*(transpose(k3)*b3);

% equation(39)-(50)
g1 = f13*f22;
g2 = f13*f25 - f15*f23;
g3 = f11*f23 - f13*f21;
g4 = -f13*f24;
g5 = f11*f22;
g6 = f11*f25 - f15*f21;
g7 = -f15*f24;
a0 = g7^2 - g2^2 - g4^2;
a1 = 2*(g6*g7 - g1*g2 - g3*g4);
a2 = g6^2 + 2*g5*g7 + g2^2 + g4^2 - g1^2 - g3^2;
a3 = 2*(g5*g6 + g1*g2 + g3*g4);
a4 = g5^2 + g1^2 + g3^2;

% equation(38)
[x1,x2,x3,x4] = Quartic_function(a4,a3,a2,a1,a0);
cos_theta_1_candidate = [x1;x2;x3;x4];
Inx = find(isreal(cos_theta_1_candidate));
cos_theta_1 = cos_theta_1_candidate(Inx,:);
cos_theta_1
sin_theta_1 = zeros(size(cos_theta_1));
cos_sin_theta_3 = zeros(2,length(cos_theta_1));
cos_theta_3 = zeros(size(cos_theta_1));
sin_theta_3 = zeros(size(cos_theta_1));
R = zeros(3,3,length(cos_theta_1));
t = zeros(3,1,length(cos_theta_1));
for i = 1:length(cos_theta_1)
    sin_theta_1(i) = sign(transpose(k3)*b3)*(1-cos_theta_1(i))^0.5;
    cos_sin_theta_3(:,i) = sin_theta_1(i)/(g5*cos_theta_1(i)^2+g6*cos_theta_1(i)+g7)*[g1*cos_theta_1(i)+g2;g3*cos_theta_1(i)+g4];
    cos_theta_3(i) = cos_sin_theta_3(1,i);
    sin_theta_3(i) = cos_sin_theta_3(2,i);
    R_prim = [k1 k3_prim_prim cross(k1,k3_prim_prim)];
    R_prim_prim = [b1 k3 cross(b1,k3)]';
    e1 = [1;0;0];
    R_e1_theta1 = cos_theta_1(i)*eye(3) - sin_theta_1(i)*skew(e1) + (1-cos_theta_1(i))*e1*e1';
    e2 = [0;1;0];
    R_e2_theta3 = cos_theta_3(i)*eye(3) - sin_theta_3(i)*skew(e2) + (1-cos_theta_3(i))*e2*e2';
    R(:,:,i) = R_prim*R_e1_theta1*R_e2_theta3*R_prim_prim;
    t(:,:,i) = p3 - sigma*sin_theta_1(i)/transpose(k3)*b3*R(:,:,i)*b3;
end

[R_true,t_true,idx]=TruePose(R,t,p4,b4,K);
