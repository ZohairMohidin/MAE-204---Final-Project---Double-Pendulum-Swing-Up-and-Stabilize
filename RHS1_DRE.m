function dX=RHS1_DRE_control(X,E,A,B,R,Q);
% keyboard
% dX = - (E')^-1 \ A' X - X * A * (E)^-1 + X * B * R^-1 * B' * X - (E')^-1 * Q * E^-1;
dX = E'\(-A'*X*E - E'*X*A + E'*X*B*R^-1*B'*X*E+Q)\E;
end