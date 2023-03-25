function x=RK4_x_LQR(A,B,K,T,x0)
% T =  total time to run RK4
t = 0:.01:T;
x=zeros(6,length(t));
x(:,1)=x0;
for i=1:length(t)
    h=.01;
        f1 = CS(x(:,i), A, B,K);
        f2 = CS(x(:,i)+f1/2*h,A,B,K);
        f3 = CS(x(:,i)+f2/2*h, A, B,K);
        f4 = CS(x(:,i)+f3*h, A,B,K);
        x(:,i+1) = x(:,i) + h*(f1/6+(f2+f3)/3+f4/6);
end
    function dx = CS(x,A,B,K)
        dx = (A+B*K)*x;
    end
end
