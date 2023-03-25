function dX=LHS_ARE(P,E,A,B,Q2,Q1);
% keyboard
dX = E'\(A'*P*E + E'*P*A - E'*P*B*Q2^-1*B'*P*E+Q1)\E;
end

for i=1:size(A,3)-1
    h=.01;
        f1 = LHS_ARE(X(:, :, i), E(:, :, i), A(:, :, i), C', Q2, Q1);
        f2 = LHS_ARE(X(:, :, i)+f1/2*h, E(:, :, i), A(:, :, i), C', Q2, Q1);
        f3 = LHS_ARE(X(:, :, i)+f2/2*h, E(:, :, i), A(:, :, i), C', Q2, Q1);
        f4 = LHS_ARE(X(:, :, i)+f3*h, E(:, :, i), A(:, :, i), C', Q2, Q1);
        P(:, :, i+1) = P(:, :, i) + h*(f1/6+(f2+f3)/3+f4/6);
end