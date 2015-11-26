function A = sic3()

d=3;

Ap=[-0.9307 + 0.0447i   0.1637 + 0.1717i  -0.2665 - 0.0668i
   0.5512 - 0.0486i   0.0292 - 0.6296i  -0.3378 - 0.4272i
   0.0822 - 0.6981i   0.2644 + 0.2239i   0.1597 + 0.6003i
   0.4229 + 0.1086i   0.4179 - 0.4264i  -0.5166 + 0.4313i
   0.3356 - 0.1249i  -0.7845 + 0.4961i   0.0276 - 0.0972i
  -0.4670 - 0.0360i   0.2255 - 0.1700i  -0.1919 + 0.8149i
  -0.1277 - 0.3374i  -0.4405 - 0.4881i   0.1628 - 0.6412i
  -0.6485 + 0.3144i  -0.6551 - 0.1341i   0.0862 + 0.1617i
   0.2033 - 0.3198i   0.2759 + 0.4163i   0.7432 - 0.2337i];


a=0;
for i = 1:d^2
    for ii = 1:d^2
        a=a+1;
        proj = kron(Ap(i,:),Ap(ii,:));
        A(a,1:d^4) = kron(proj,conj(proj));
    end
end
