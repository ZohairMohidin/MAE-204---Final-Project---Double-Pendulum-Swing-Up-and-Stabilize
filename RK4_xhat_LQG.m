function xhat=RK4_xhat_LQG(A,B,C,K,L,T,x)
% T =  total time to run RK4
t = 0:.01:T;
xhat=zeros(6,length(t));
for i=1:length(t)
    h=.01;
        f1 = CS(xhat(:,i), A,B,C,K,L,x(:,i));
        f2 = CS(xhat(:,i)+f1/2*h,A,B,C,K,L,x(:,i));
        f3 = CS(xhat(:,i)+f2/2*h, A, B,C,K,L,x(:,1));
        f4 = CS(xhat(:,i)+f3*h, A,B,C,K,L,x(:,i));
        xhat(:,i) = xhat(:,i) + h*(f1/6+(f2+f3)/3+f4/6);
end
    function dx = CS(xhat,A,B,C,K,L,x)
        dx = (A-L'*C)*xhat + B*K*xhat + L*C*x;
    end
end