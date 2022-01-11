clear all
close all
clc
rng(107)
load('SRF_S2')

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

%% CCRB

[~,J_13,J_23] = vec_unfold([dim2(1) dim2(2) dim1(3)]);
Sa = J_13*kron(pw_kr(C1,B2,R),eye(dim2(1)));
Sb = J_23*kron(pw_kr(C1,A2,R),eye(dim2(2)));
Sc = kron(pw_vec2(A2,B2,R),eye(dim1(3)));
tmp1 = [A2(:); B2(:); C1(:)];  ind = sort([find(tmp1==1); find(tmp1==0)]);
M = eye((dim2(1)+dim2(2))*L*R + dim1(3)*R); M(ind,:) = [];
H = [Sa Sb Sc]*M';

for s=1:length(SNR1)
    s
    [~,CCRB] = ccrb_allatonce(A1,B1,C1,A2,B2,C2,P1,P2,Pm,sigma_n1(s),sigma_n2);
    CCRB_rec = H*CCRB*H';
    
    CCRB_A2(s) = sum(diag(CCRB(1:(dim2(1)-L)*L*R,1:(dim2(1)-L)*L*R)));
    CCRB_B2(s) = sum(diag(CCRB((dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R,(dim2(1)-L)*L*R+1:(dim2(1)+dim2(2)-L)*L*R)));
    CCRB_C1(s) = sum(diag(CCRB((dim2(1)+dim2(2)-L)*L*R+1:end,(dim2(1)+dim2(2)-L)*L*R+1:end)));
    CCRB_Y(s) = sum(diag(CCRB_rec));
end

%% Pre-allocation

Nreal = 200; Ninit = 10; Niter = 500;

se_C1_c = zeros((dim1(3)-1)*R,Nreal,length(SNR1));
se_A2_c = zeros((dim2(1)-L)*L*R,Nreal,length(SNR1));
se_B2_c = zeros(dim2(2)*L*R,Nreal,length(SNR1));

mse_C1_c = zeros(1,length(SNR1));
mse_A2_c = zeros(1,length(SNR1));
mse_B2_c = zeros(1,length(SNR1));

%% Simulations

options.OutputFormat = 'cpd'; options.MaxIter = 200;
P = Nreal*Ninit*length(sigma_n1); p=0; %Progress indicator 

for s=1:length(sigma_n1)
    for n=1:Nreal
    
        %Generate noisy tensors
        Y1 = X1+sigma_n1(s)*randn(dim1);
        Y2 = X2+sigma_n2*randn(dim2);

        %---------------------------    

        %Initial objective values
        obj1 = 10^50; obj2 = 10^50; obj3 = 10^50; obj4 = 10^50;
        for i=1:Ninit

            %Init
            A10 = [kron(ones(1,R),eye(L)); randn(dim1(1)-L,L*R)]; 
            B10 = randn(size(B1));
            C10 = randn(size(C1)); C10(1,:) = 1;
            A20 = [kron(ones(1,R),eye(L)); randn(dim2(1)-L,L*R)]; 
            B20 = randn(size(B2));
            C20 = randn(size(C2)); C20(1,:) = 1;


            % Coupled LL1 ALS algo
            [A2_cu,B2_cu,C1_cu] = coupled_ll1_als(Y1,Y2,B20,C10,P1,P2,Pm,(sigma_n2/sigma_n1(s))^2,Niter);
            cost3 = (sigma_n2/sigma_n1(s))^2*frob(X1-ll1gen({P1*A2_cu,P2*B2_cu,C1_cu},L*ones(1,R)),'squared') + ...
                frob(X2-ll1gen({A2_cu,B2_cu,Pm*C1_cu},L*ones(1,R)),'squared');
            
            if cost3<obj3
                obj3 = cost3;
                A2_hat_c = A2_cu; B2_hat_c = B2_cu; C1_hat_c = C1_cu;
            end
            

            p=p+1; clc
            fprintf('Progress %g %%',(p/P)*100)

        end

        %---------------------------    

        % Rescale
        
        %------

        D = B2_hat_c'; %tmp = [];
        D = pw_kr(eye(R),A2_hat_c(1:L,:),R)*D;
        B2_hat_c = D';
        B2_hat_c = B2_hat_c.*kron(repmat(C1_hat_c(1,:),dim2(2),1),ones(1,L));
        C1_hat_c = C1_hat_c.*repmat(1./C1_hat_c(1,:),dim1(3),1);
        A2_hat_c = A2_hat_c*inv(pw_kr(eye(R),A2_hat_c(1:L,:),R));

  

        %Permute
        %------
        
        ind = sort_columns(C1,C1_hat_c)';
        C1_hat_c = C1_hat_c(:,ind);
        tmpA = []; tmpB = [];
        for r=1:R
            tmpA(:,(r-1)*L+1:r*L) = A2_hat_c(:,(ind(r)-1)*L+1:ind(r)*L);
            tmpB(:,(r-1)*L+1:r*L) = B2_hat_c(:,(ind(r)-1)*L+1:ind(r)*L);
        end
        A2_hat_c = tmpA; B2_hat_c = tmpB;
        Y_hat_c = zeros([dim2(1) dim2(2) dim1(3)]);
        for r=1:R
            Y_hat_c = Y_hat_c + outprod(A2_hat_c(:,(r-1)*L+1:r*L)*B2_hat_c(:,(r-1)*L+1:r*L)',C1_hat_c(:,r));
        end

        %---------------------------    

        %Squared errors
        
        err = A2(L+1:end,:) - A2_hat_c(L+1:end,:); se_A2_c(:,n,s) = (err(:)).^2;
        se_B2_c(:,n,s) = (B2(:)-B2_hat_c(:)).^2;
        err = C1(2:end,:) - C1_hat_c(2:end,:); se_C1_c(:,n,s) = (err(:)).^2; 

        se_Y_c(:,n,s) = (Y(:)-Y_hat_c(:)).^2;
        
    end
    
    %MSE 
    mse_A2_c(s) = sum(mean(se_A2_c(:,:,s),2));
    mse_B2_c(s) = sum(mean(se_B2_c(:,:,s),2));
    mse_C1_c(s) = sum(mean(se_C1_c(:,:,s),2));
    
    mse_Y_c(s) = sum(mean(se_Y_c(:,:,s),2));
end


%% Figures

figure(1)
subplot(1,2,1)
semilogy(SNR1,CCRB_C1+CCRB_A2+CCRB_B2,'k--'); hold on
semilogy(SNR1,mse_C1_c+mse_A2_c+mse_B2_c,'bo'); 
subplot(1,2,2)
semilogy(CCRB_Y,'k--'); hold on; semilogy(mse_Y_c,'bo')
saveas(gcf,'figures/fig1.fig')

save('mse_ll1','mse_Y_c','mse_C1_c','mse_A2_c','mse_B2_c')
