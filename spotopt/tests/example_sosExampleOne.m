% Build a homogeneous polynomial, then test the sphere.
n = 4;
x = msspoly('x',n);

d = 10;
mn = monomials(x,d);

p = rand(length(mn),1)'*mn;


pr = spotsosprog;
[pr,f] = pr.newFree(1);
[pr,q] = pr.newFreePoly(monomials(x,0:d-2));
pr = pr.withSOS(p + f + q*(1-x'*x));

solver = spotsolversedumi;

sol = pr.minimize(solver,f);

N = 1000;
X = randn(n,N);
X = X./(repmat(sqrt(sum(X.^2,1)),n,1));

[min(msubs(p,x,X)) -sol.eval(f)]