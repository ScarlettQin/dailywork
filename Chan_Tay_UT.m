function [x_Chan_Tay,x_Chan_Tay_var] = Chan_Tay_UT(r_i1,M,x,cov_R,sigma)
% r_i1 TDOA量测
% M 基站个数
% x 基站位置
% cov_R 噪声方差

%% 最小二乘计算目标位置
for i = 1:M-1
    G1(i,:) = -[(x(:,i+1)-x(:,1))' r_i1(i)];
    h1(i,:) = 0.5*(r_i1(i)^2-x(:,i+1)'*x(:,i+1)+x(:,1)'*x(:,1))+0.5*sigma^2;
end

%%  第一步计算
x_chan1 = inv(G1'*inv(cov_R)*G1)*G1'*inv(cov_R)*h1;

for k = 1:3
    for i = 1:M
        r(i) = sqrt(sum((x(:,i)-x_chan1(1:2)).^2));
    end
    B1 = diag(r(2:M));
    Phi1 = UT_Conversion(M,cov_R,r);
    
    for i = 1:M-1
        G1_T(i,:) = [((x(:,1)-x_chan1(1:2))./r(1)-(x(:,i+1)-x_chan1(1:2))./r(i+1))' 0];
        h1_T(i,:) = r_i1(i) - r(i+1) + r(1) + G1_T(i,1:2)*x_chan1(1:2);
    end
    
    G1_aug = [G1;G1_T];
    h1_aug = [h1;h1_T];
    
    Phi1_aug1 = blkdiag(Phi1,cov_R);
    Phi1_aug = UT_fir_cov(M,r,cov_R);
    x_chan2 = pinv(G1_aug'*pinv(Phi1_aug)*G1_aug)*G1_aug'*pinv(Phi1_aug)*h1_aug;
    cov_x_chan2 = pinv(G1_aug'*pinv(Phi1_aug)*G1_aug);
    x_chan1 = x_chan2;
end

%% 第二步计算
G2 = [1 0; 0 1;1 1];
h2 = [(x_chan2(1,1)-x(1,1))^2; (x_chan2(2,1)-x(2,1))^2; (x_chan2(3,1))^2];
B2 = diag([x_chan2(1)-x(1,1) x_chan2(2,1)-x(2,1) x_chan2(3,1)]);
Phi2 = B2*cov_x_chan2*B2;
x_chan3 = inv(G2'*inv(Phi2)*G2)*G2'*inv(Phi2)*h2;
x_chan31 = sqrt(x_chan3)+x(:,1);


%% 另一种第二步
G3 = [1 0; 0 1; x(:,1)'];
h3 = [x_chan2(1);x_chan2(2);(x_chan2(1:2)'*x_chan2(1:2)+x(:,1)'*x(:,1)-x_chan2(3)^2)*0.5];
B3 = [1 0 0; 0 1 0; x_chan2(1:2)' x_chan2(3)];
Phi3 = B3*cov_x_chan2*B3';
x_chan4 = pinv(G3'*pinv(Phi3)*G3)*G3'*pinv(Phi3)*h3;

x_Chan_Tay = x_chan31;
x_Chan_Tay_var = x_chan4;

end



