function [A_hat,B_hat,C_hat] = coupled_ll1_als(Y1,Y2,B_hat,C_hat,P1,P2,Pm,lamda,Niter)

R = size(C_hat,2);
L = size(B_hat,2)/R;

cost(1) = Inf; diff_cost(1) = Inf; 
n=1;

while (n<Niter) && (diff_cost(n) > eps) && (cost(n) > eps)
    
    n=n+1;
    
    % Update A
    temp1h=pw_kr(C_hat,P2*B_hat,R);
    temp1m=pw_kr(Pm*C_hat,B_hat,R);
    Kp=pw_kr(Pm*C_hat,B_hat,R)'*pw_kr(Pm*C_hat,B_hat,R);
    K= pw_kr(C_hat,P2*B_hat,R)'*pw_kr(C_hat,P2*B_hat,R);
    inv_K=pinv(K);
    As=lamda*(P1'*P1);
    Bs=Kp*inv_K;
    Cs=(lamda*P1'*tens2mat(Y1,[],1)'*temp1h+tens2mat(Y2,[],1)'*temp1m)*inv_K;
    A_hat=sylvester(full(As),Bs,Cs);
    
    % Update B
    temp1h=pw_kr(C_hat,P1*A_hat,R);
    temp1m=pw_kr(Pm*C_hat,A_hat,R);
    Kp=pw_kr(Pm*C_hat,A_hat,R)'*pw_kr(Pm*C_hat,A_hat,R);
    K= pw_kr(C_hat,P1*A_hat,R)'*pw_kr(C_hat,P1*A_hat,R);
    inv_K=pinv(K);
    As=lamda*(P2'*P2);
    Bs=Kp*inv_K;
    Cs=(lamda*P2'*tens2mat(Y1,[],2)'*temp1h+tens2mat(Y2,[],2)'*temp1m)*inv_K;
    B_hat=sylvester(full(As),Bs,Cs);
    
    % Update C
    temp3h=pw_vec2(P1*A_hat,P2*B_hat,R);
    temp3m=pw_vec2(A_hat,B_hat,R);
    Kp=lamda*pw_vec2(P1*A_hat,P2*B_hat,R)'*pw_vec2(P1*A_hat,P2*B_hat,R);
    K= pw_vec2(A_hat,B_hat,R)'*pw_vec2(A_hat,B_hat,R);
    inv_K=pinv(K);
    As=Pm'*Pm;
    Bs=Kp*inv_K;
    Cs=(lamda*tens2mat(Y1,[],3)'*temp3h+Pm'*tens2mat(Y2,[],3)'*temp3m)*inv_K;
    C_hat=sylvester(full(As),Bs,Cs);
    for r=1:R
        C_hat(:,r) = C_hat(:,r)/norm(C_hat(:,r));
    end
    
    cost(n) = lamda*frob(Y1-ll1gen({P1*A_hat,P2*B_hat,C_hat},L*ones(1,R)),'squared') + ...
                frob(Y2-ll1gen({A_hat,B_hat,Pm*C_hat},L*ones(1,R)),'squared');
    diff_cost(n) = cost(n-1)-cost(n);
    
end
    
end

