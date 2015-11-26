function [A ]= MUB3GD()

d=3;

Ap =  [ -0.4393 + 0.0125i   0.0001 + 0.1178i  -0.3397 - 0.8231i
   0.3778 + 0.2595i  -0.0748 - 0.8264i   0.0025 - 0.3183i
  -0.4275 - 0.6433i   0.1993 - 0.5080i  -0.2232 + 0.2361i
  -0.5821 - 0.4731i   0.2602 + 0.0114i   0.5578 + 0.2417i
   0.1762 + 0.6204i   0.2975 + 0.3685i   0.3513 + 0.4861i
   0.1094 - 0.0991i  -0.6064 + 0.5829i   0.4287 - 0.2948i
   0.3163 - 0.0741i   0.7351 + 0.2995i  -0.0903 + 0.5062i
  -0.2152 + 0.2357i   0.0473 + 0.4899i   0.8089 - 0.0403i
  -0.2796 + 0.8451i   0.1452 - 0.3252i  -0.1200 + 0.2576i
  -0.3658 - 0.2897i   0.5970 - 0.4100i   0.0483 - 0.5054i
  -0.7163 + 0.2585i  -0.1722 + 0.5543i   0.1049 - 0.2685i
  -0.2310 - 0.3862i   0.2442 + 0.2829i  -0.7174 + 0.3785i ];



for i=1:(d+1)
   M{i}  = Ap((1:d)+d*(i-1) , 1:d);
end

a=0;




    for i = 1:(d+1)
        for k=1:d
            for ii = 1:(d+1)
                for kk=1:d
                    a=a+1;
                    proj = kron(M{i}(k,:),(M{ii}(kk,:)));
                    A(a,1:d^4) = kron(proj,conj(proj));
                    
                end
            end
        end
    end





