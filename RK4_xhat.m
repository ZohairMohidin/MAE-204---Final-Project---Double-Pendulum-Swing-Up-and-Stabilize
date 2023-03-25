function [x_hat,x_hat_sim]=RK4_xhat(A,B,C,E,x0,x_k,u_k,K,L)

% x_k is a row vector of 6 not 9


x_hat=zeros(size(x0));
  h=.01;
for i=1:size(x_k,2)-1
        f1 = CS(x_hat(:, i),E(:, :, i),A(:, :, i), B,C ,u_k(i),K(:,:,i),L(:,:,i),x_k(:,i));
        f2 = CS(x_hat(:, i)+f1/2*h,E(:, :, i),A(:, :, i), B,C ,u_k(i),K(:,:,i),L(:,:,i),x_k(:,i));
        f3 = CS(x_hat(:, i)+f2/2*h,E(:, :, i),A(:, :, i), B,C ,u_k(i),K(:,:,i),L(:,:,i),x_k(:,i));
        f4 = CS(x_hat(:, i)+f3*h,E(:, :, i),A(:, :, i), B,C ,u_k(i),K(:,:,i),L(:,:,i),x_k(:,i));
        x_hat(:,i+1) = x_hat(:, i) + h*(f1/6+(f2+f3)/3+f4/6);
        x_hat_sim(:,i) = x_k(:,i)+h*x_hat(:,i);
end
x_hat_sim(:,length(x_k))= x_k(:,end)+h*x_hat(:,end);
    function dx=CS(x_hat,E,A,B,C,u_k,K,L,x_k)
dx = (E\(A+L*C))*x_hat + (2*L*C - E\B*K)*x_k +(L*C + E\B*K)*x_hat + E\B*u_k;
    end
end
