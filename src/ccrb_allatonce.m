function [FIM,CCRB] = ccrb_allatonce(A1,B1,C1,A2,B2,C2,P1,P2,Pm,sigma_n1,sigma_n2)


%Extract dimensions
dim1 = [size(A1,1) size(B1,1) size(C1,1)]; 
dim2 = [size(A2,1) size(B2,1) size(C2,1)];
R = size(C2,2); L = size(A2,2)/R;
%Permutation matrices
[~,J_13_1,J_23_1] = vec_unfold(dim1);
[~,J_13_2,J_23_2] = vec_unfold(dim2);
%Covariance matrix
Sigma = blkdiag((1/sigma_n1^2)*eye(prod(dim1)),(1/sigma_n2^2)*eye(prod(dim2)));

%Selection matrix
tmp1 = [A2(:); B2(:); C1(:);]; ind = sort([find(tmp1==1); find(tmp1==0)]);
M = eye((dim2(1)+dim2(2))*L*R + dim1(3)*R);  M(ind,:) = [];

%Tensor unfoldings

Sa1 = J_13_1*kron(pw_kr(C1,P2*B2,R),P1);
Sb1 = J_23_1*kron(pw_kr(C1,P1*A2,R),P2);
Sc1 = kron(pw_vec2(P1*A2,P2*B2,R),eye(dim1(3))); 
Sa2 = J_13_2*kron(pw_kr(Pm*C1,B2,R),eye(dim2(1)));
Sb2 = J_23_2*kron(pw_kr(Pm*C1,A2,R),eye(dim2(2)));
Sc2 = kron(pw_vec2(A2,B2,R),Pm); 


Mu = [Sa1 Sb1 Sc1;
      Sa2 Sb2 Sc2];

FIM = M*Mu'*Sigma*Mu*M';
CCRB = inv(FIM);

end

