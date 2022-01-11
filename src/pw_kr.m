function X = pw_kr(A,B,R)

La = size(A,2)/R;
Lb = size(B,2)/R;
X = [];

for r=1:R
    X = [X kron(A(:,(r-1)*La+1:r*La),B(:,(r-1)*Lb+1:r*Lb))];
end

end

