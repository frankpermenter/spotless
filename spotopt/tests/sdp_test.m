pr = spotsdp;


n = 10;
m = 4;
A = randn(n,m);

[pr,p] = pr.newPos(1);
[pr,X] = pr.withPSD([ p*eye(n) A ; A' p*eye(m)]);

solver = spotsolversedumi(struct('fid',1));

sol = solver.minimize(pr,p);