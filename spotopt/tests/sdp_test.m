pr = spotsdp;


n = 10;
m = 4;
A = randn(n,m);

[pr,p] = pr.newFree(1);
[pr] = pr.withPos(p);
[pr] = pr.withPSD([ p*eye(n) A ; A' p*eye(m)]);

solver = spotsolversedumi(struct('fid',1));

sol = solver.minimize(pr,p);
sol = solver.minimizeDualForm(pr,p);

solver = spotsolversdpnal();
sol = solver.minimizeDualForm(pr,p);

solver = spotsolversdpt3_4();
sol = solver.minimizeDualForm(pr,p);