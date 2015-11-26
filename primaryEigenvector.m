function  [primary maxEig] = primaryEigenvector(rhof);
[a b] = eig(rhof);
pos = find(sum(b)==max(sum(b)));
primary = a(:,pos);
maxEig = max(diag(b));