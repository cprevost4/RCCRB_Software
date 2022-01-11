function X = pw_vec2(A,B,R)

L = size(A,2)/R;
X = [];

for r=1:R
    tmp = A(:,(r-1)*L+1:r*L)*B(:,(r-1)*L+1:r*L)';
    X = [X tmp(:)];
end

end

