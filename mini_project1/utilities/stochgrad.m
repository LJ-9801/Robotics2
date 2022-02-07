function J=stochgrad(x,fun)

n=length(x);
N=10*n;

e=norm(x)*1e-6;
y0=feval(fun,x);
dx=e*randn(length(x),N);

for i=1:10*n
    y(:,i)=feval(fun,x+dx(:,i));
end

J=(y-y0)*pinv(dx);

end


