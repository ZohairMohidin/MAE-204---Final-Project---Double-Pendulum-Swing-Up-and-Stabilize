function X=RHS(X,E,A,B,R,Q)
% keyboard

X = E'\(-A'*X*E - E'*X*A + E'*X*B*R^-1*B'*X*E+Q)\E;
for i=size(A,3):-1:2
    h=.01;
        f1 = RHS(X(:, :, i), E(:, :, i), A(:, :, i), B, R, Q);
        f2 = RHS(X(:, :, i)-f1/2*h, E(:, :, i), A(:, :, i), B, R, Q);
        f3 = RHS(X(:, :, i)-f2/2*h, E(:, :, i), A(:, :, i), B, R, Q);
        f4 = RHS(X(:, :, i)-f3*h, E(:, :, i), A(:, :, i), B, R, Q);
        X(:, :, i-1) = X(:, :, i) - h*(f1/6+(f2+f3)/3+f4/6);
end

end
