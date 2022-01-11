function [FIM, CRB] = crb(dim,A,B,C,sigma_n1)

R = size(C,2); L = size(A,2)/R;
[~,J_13,J_23] = vec_unfold(dim);

tmp = [A(:); B(:); C(:)];
ind = sort([find(tmp==1); find(tmp==0)]);
M = eye(dim(3)*R + (dim(1)+dim(2))*L*R); M(ind,:) = [];

Sa = J_13*kron(pw_kr(C,B,R),eye(dim(1)));
Sb = J_23*kron(pw_kr(C,A,R),eye(dim(2)));
Sc = kron(pw_vec2(A,B,R),eye(dim(3)));
mat = [Sa Sb Sc];
Sigma = (1/sigma_n1^2)*eye(prod(dim));

FIM = M*mat'*Sigma*mat*M';
CRB = pinv(FIM);


end

