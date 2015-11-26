
%A is the matrix of Pauli measurements
%P is the corresponding count vector with noise that is simulated well.
%A full experiment is simulated here.
%Each joint projector is "measured" and noise is added to these. Then the
%Pauli outcomes are built. It was hard and painful. Do not try to reproduce
%this code. Just write it again from scratch if you need to.

function [AA P] = make3QubitHVADRLandPauli(rho)


d=2^3;
x = reshape(rho,[d^2 1]);
%%%%%%%%% PAULi %%%%%%%%%%%%%
s{1} = [1/sqrt(2) 1/sqrt(2); 
    1/sqrt(2) -1/sqrt(2)];

s{2} = [ 1/sqrt(2) 1i/sqrt(2); 
    1/sqrt(2) -1i/sqrt(2)];

s{3} = [1 0; 
    0 1];

p123 = zeros(3,3,3);
p1 = zeros(3,3,3);
p2=zeros(3,3,3);
p3=zeros(3,3,3);
p12=zeros(3,3,3);
p13=zeros(3,3,3);
p23=zeros(3,3,3);
p0 = zeros(3,3,3);

A=[];n=0;
itNumber= 0;
id = 0;
for i=1:3
    mub1 = s{i};
    for ii=1:3
        mub2 = s{ii};
        for iii=1:3
            mub3 = s{iii};
            for k=1:2
                v1 = mub1(k,:);
                for kk=1:2
                    v2 = mub2(kk,:);
                    for kkk=1:2
                        itNumber = itNumber + 1;
                        v3 = mub3(kkk,:);
                        v = kron(kron(v1,v2),v3);
                        vv = kron(v,conj(v));
                        
                        p = abs(vv*x);
                        pExact(itNumber) = p;
                         pNoisy = poissrnd(p);
%                         pNoisy = p; %no noise
                        A = [A;v];
                
                        ka=k;
                        kka=kk;
                        kkka=kkk;
                        
                        p123(i,ii,iii) = p123(i,ii,iii)+pNoisy*(ka-3/2)*2^3*(kka-3/2)*(kkka-3/2);
                        
                        p12(i,ii,iii) = p12(i,ii,iii)+pNoisy*(ka-3/2)*2^2*(kka-3/2);
                        p13(i,ii,iii) = p13(i,ii,iii)+pNoisy*(ka-3/2)*2^2*(kkka-3/2);
                        p23(i,ii,iii) = p23(i,ii,iii)+pNoisy*(kka-3/2)*2^2*(kkka-3/2);
                        
                        p1(i,ii,iii) = p1(i,ii,iii)+pNoisy*(ka-3/2)*2;
                        p2(i,ii,iii) = p2(i,ii,iii)+pNoisy*(kka-3/2)*2;
                        p3(i,ii,iii) = p3(i,ii,iii)+pNoisy*(kkka-3/2)*2;
                        
                        p0(i,ii,iii) =  p0(i,ii,iii) + pNoisy;
                    end
                end
            end
        end
    end
end


q0 = sum(sum(sum(p0)))/3^3;
q1 = sum(sum(p1,2),3)/3^2;
q2 = sum(sum(p2,1),3)/3^2;
q3 = sum(sum(p3,1),2)/3^2;

q1 = reshape(q1,[3 1 1]);
q2 = reshape(q2,[3 1 1]);
q3 = reshape(q3,[3 1 1]);

q12 = sum(p12,3)/3;
q13 = sum(p13,2)/3;
q23 = sum(p23,1)/3;

s= size(q12);
ell= prod(s);

q12 = reshape(q12,[ell 1 1]);
q13 = reshape(q13,[ell 1 1]);
q23 = reshape(q23,[ell 1 1]);

s= size(p123);
ell= prod(s);

p123 = permute(p123,[3 2 1]);
q123 = reshape(p123,[ell 1 1]);

sig{1} = eye(2);
sig{2} = [0 1; 1 0];
sig{3} = [0 1i; -1i 0];
sig{4} = [1 0; 0 -1];

P = [q0;-q3;-q2;-q1;q23;q13;q12;-q123];

a=0;
AA = reshape(eye(d),[1 d^2]);
for pos=3:-1:1
    t = [1 1 1];
    for i=1:3
        t(pos) = i+1;
        a=a+1;
        Atemp = kron(kron(sig{t(1)},sig{t(2)}),sig{t(3)});
       vvv = reshape(Atemp,[1 d^2]);
         AA = [AA;vvv];
    end
end

a=0;
pos = [3 2; 3 1; 2 1];
for po=1:3
    t = [1 1 1];
    for i=1:3
        t(pos(po,1)) = t(pos(po,1))+1;
        t(pos(po,2)) =1;
        for ii=1:3
            t(pos(po,2)) = t(pos(po,2))+1;
            a=a+1;
            Atemp = kron(kron(sig{t(1)},sig{t(2)}),sig{t(3)});
            vvv = reshape(Atemp,[1 d^2]);
            AA = [AA;vvv];
        end
    end
end

a=0;
for i=2:4
    for ii=2:4
        for iii=2:4
            M = kron(kron(sig{i},sig{ii}),sig{iii});
            vvv = reshape(M,[1 d^2]);
            AA = [AA; vvv];
            a=a+1;
%             P(a) = p123(i-1,ii-1,iii-1);
        end
    end
end

 
Pkron = [ real(AA*x) P]




