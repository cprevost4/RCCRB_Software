clear all
close all
clc
rng(107)
load('SRF_S2')
load('random_vector')

%% Generate model

dim1 = [6 6 30]; dim2 = [24 24 6]; 
L = 2; R = 2;

A2 = [kron(ones(1,R),eye(L)); randn(dim2(1)-L,L*R)]; 
B2 = randn(dim2(2),L*R);
C1 = randn(dim1(3),R); C1(1,:) = 1;

q = 3; phi = gauss_kernel(q); phi = phi/norm(phi); %Cond=1 if q=d
H = toeplitz([phi(1), zeros(1,dim2(1)-1)], [phi zeros(1,dim2(1)-q)]);
d = 4; S = eye(d*dim1(1)); S = S(1:d:end,:);
P1 = S*H; P2 = P1;
Pm = Pm(3:8,25:54);

A1 = P1*A2; B1 = P2*B2; C2 = Pm*C1;
S1 = pw_vec2(A1,B1,R); S2 = pw_vec2(A2,B2,R);

X1 = zeros(dim1);
for r=1:R
    X1 = X1 + outprod(A1(:,(r-1)*L+1:r*L)*B1(:,(r-1)*L+1:r*L)',C1(:,r));
end
X2 = zeros(dim2);
for r=1:R
    X2 = X2 + outprod(A2(:,(r-1)*L+1:r*L)*B2(:,(r-1)*L+1:r*L)',C2(:,r));
end
Y = zeros([dim2(1) dim2(2) dim1(3)]); % SRI
for r=1:R
    Y = Y + outprod(A2(:,(r-1)*L+1:r*L)*B2(:,(r-1)*L+1:r*L)',C1(:,r));
end

SNR1 = 5:5:60; %Noise on first tensor
for s=1:length(SNR1)
    sigma_n1(s) = 10^(-SNR1(s)/10);
end
SNR2 = 20; %Noise on first tensor
for s=1:length(SNR2)
    sigma_n2(s) = 10^(-SNR2(s)/10);
end

%% Conditional CCRB 1

%vector1 = normrnd(1,0.3,[1 size(Pm,2)]);

[~,J_13,J_23] = vec_unfold([dim2(1) dim2(2) dim1(3)]);
Sa = J_13*kron(pw_kr(C1,B2,R),eye(dim2(1)));
Sb = J_23*kron(pw_kr(C1,A2,R),eye(dim2(2)));
Sc = kron(pw_vec2(A2,B2,R),eye(dim1(3)));
tmp1 = [A2(:); B2(:); C1(:)];  ind = sort([find(tmp1==1); find(tmp1==0)]);
M = eye((dim2(1)+dim2(2))*L*R + dim1(3)*R); M(ind,:) = [];
H = [Sa Sb Sc]*M';

for s=1:length(SNR1)
    s
    [~,CCRB] = ccrb_allatonce(A1,B1,C1,A2,B2,C2,P1,P2,vector1.*Pm,sigma_n1(s),sigma_n2);
    CCRB_rec = H*CCRB*H';
    
    CCRB1_A2(s) = sum(diag(CCRB(1:(dim2(1)-L)*L*R,1:(dim2(1)-L)*L*R)));
    CCRB1_B2(s) = sum(diag(CCRB((dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R,(dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R)));
    CCRB1_C1(s) = sum(diag(CCRB((dim2(1)+dim2(2)-L)*L*R+1:end,(dim2(1)+dim2(2)-L)*L*R+1:end)));
    CCRB1_Y(s) = sum(diag(CCRB_rec));
end

%% Conditional CCRB 2

%vector2 = normrnd(1,0.3,[1 size(Pm,2)]);

[~,J_13,J_23] = vec_unfold([dim2(1) dim2(2) dim1(3)]);
Sa = J_13*kron(pw_kr(C1,B2,R),eye(dim2(1)));
Sb = J_23*kron(pw_kr(C1,A2,R),eye(dim2(2)));
Sc = kron(pw_vec2(A2,B2,R),eye(dim1(3)));
tmp1 = [A2(:); B2(:); C1(:)];  ind = sort([find(tmp1==1); find(tmp1==0)]);
M = eye((dim2(1)+dim2(2))*L*R + dim1(3)*R); M(ind,:) = [];
H = [Sa Sb Sc]*M';

for s=1:length(SNR1)
    s
    [~,CCRB] = ccrb_allatonce(A1,B1,C1,A2,B2,C2,P1,P2,vector2.*Pm,sigma_n1(s),sigma_n2);
    CCRB_rec = H*CCRB*H';
    
    CCRB2_A2(s) = sum(diag(CCRB(1:(dim2(1)-L)*L*R,1:(dim2(1)-L)*L*R)));
    CCRB2_B2(s) = sum(diag(CCRB((dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R,(dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R)));
    CCRB2_C1(s) = sum(diag(CCRB((dim2(1)+dim2(2)-L)*L*R+1:end,(dim2(1)+dim2(2)-L)*L*R+1:end)));
    CCRB2_Y(s) = sum(diag(CCRB_rec));
end

%% RCCRB

Nreal = 500;

for s=1:length(SNR1)
    
    CCRB = 0;
    for n=1:Nreal
        [s n]
        scalar = normrnd(1,0.1,[1 size(Pm,2)]);
        
        [~,mat] = ccrb_allatonce(A1,B1,C1,A2,B2,C2,P1,P2,scalar.*Pm,sigma_n1(s),sigma_n2);
        CCRB = CCRB + mat;
    end
    
    RCCRB(:,:,s) = CCRB/Nreal;
end

% For the parameters
for s=1:length(SNR1)
    s
    RCCRB_A2(s) = sum(diag(RCCRB(1:(dim2(1)-L)*L*R,1:(dim2(1)-L)*L*R,s)));
    RCCRB_B2(s) = sum(diag(RCCRB((dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R,(dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R,s)));
    RCCRB_C1(s) = sum(diag(RCCRB((dim2(1)+dim2(2)-L)*L*R+1:end,(dim2(1)+dim2(2)-L)*L*R+1:end,s)));
    RCCRB_omega(s) = sum(diag(RCCRB(:,:,s)));
    
    mat = H*RCCRB(:,:,s)*H';
    RCCRB_Y(s) = sum(diag(mat));
end

%% Figures

load('mse_ll1')

figure(1)
subplot(2,2,1); semilogy(SNR1,CCRB1_C1,'k--'); hold on;
semilogy(SNR1,CCRB2_C1,'r--'); hold on; semilogy(SNR1,mse_C1_c,'bo');
subplot(2,2,2); semilogy(SNR1,CCRB1_A2,'k--'); hold on;
semilogy(SNR1,CCRB2_A2,'r--'); hold on; semilogy(SNR1,mse_A2_c,'bo');
subplot(2,2,3); semilogy(SNR1,CCRB1_B2,'k--'); hold on;
semilogy(SNR1,CCRB2_B2,'r--'); hold on; semilogy(SNR1,mse_B2_c,'bo');
subplot(2,2,4); semilogy(SNR1,CCRB1_Y,'k--'); hold on;
semilogy(SNR1,CCRB2_Y,'r--'); hold on; semilogy(SNR1,mse_Y_c,'bo');
saveas(gcf,'figures/fig2.fig')

figure(2)
subplot(1,2,1); semilogy(SNR1,RCCRB_omega,'k--'); hold on;
semilogy(SNR1,mse_A2_c+mse_B2_c+mse_C1_c,'bo');
subplot(1,2,2); semilogy(SNR1,RCCRB_Y,'k--'); hold on;
semilogy(SNR1,mse_Y_c,'bo');
saveas(gcf,'figures/fig3.fig')
save('RCCRB','RCCRB_omega','RCCRB_Y')