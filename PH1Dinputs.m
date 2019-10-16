function [effortsp,effortsq] = PH1Dinputs(Y,mu,time,X,X1,X2,N,Xdual,option)
if option==1
    effortsq= zeros(2 * N ,1);
    effortsp= zeros(3 * N ,1);
    %--------------------------------------------effortsp
    effortsp(1)= -2*pi*sin(2*pi*(X(1)+time));
 
    k=1;
    for index1=1: N
        effortsp(3*index1-1)=  -2*pi*sin(2*pi*(Xdual(k)+time));
        k=k+1;
    end
    k=2;
    for index1=1:N-1
        effortsp(3*index1)= -2*pi*sin(2*pi*(X(k)+time));
        effortsp(3*index1+1)= -2*pi*sin(2*pi*(X(k)+time));

        k=k+1;
    end
    effortsp(3*N)= -2*pi*sin(2*pi*(X(N+1)+time));
    % -----------------------------------------------------effortsq
    effortsq(1)= -2*pi*sin(2*pi*(X(1)+time));
    
    
    k=2;
    for index1=2:2:(2*N - 1)
        effortsq(index1)=  -2*pi*sin(2*pi*(X(k)+time));
        effortsq(index1+1)= -2*pi*sin(2*pi*(X(k)+time));
        
        
        k=k+1;
        
    end
    effortsq(2*N)= -2*pi*sin(2*pi*(X(N+1)+time));
%--------------------------------------------------------Case 2
else if option==2
        effortsq= zeros(3 * N ,1);
        effortsp= zeros(4 * N ,1);
         %--------------------------------------------effortsp
        effortsp(1)=-2*pi*sin(2*pi*(X(1)+time));

        k=1;
        for index1=1 : N
            effortsp(4*index1-2)=-2*pi*sin(2*pi*(X1(k)+time));
            effortsp(4*index1-1)=-2*pi*sin(2*pi*(X2(k)+time));
            k=k+1;
        end
        k=2;
        for index1=1:N-1
            effortsp(4*index1)= -2*pi*sin(2*pi*(X(k)+time));
            effortsp(4*index1+1)= -2*pi*sin(2*pi*(X(k)+time));
            k=k+1;
        end
        effortsp(4 * N)=-2*pi*sin(2*pi*(X(N+1)+time));
        %-----------------------------------------------------effortsq
        effortsq(1)=-2*pi*sin(2*pi*(X(1)+time));
        k=1;
        for index1=1: N
            effortsq(3*index1-1)=-2*pi*sin(2*pi*(Xdual(k)+time));
            k=k+1;
        end
        k=2;
        for index1=1:N-1
            effortsq(3*index1)= -2*pi*sin(2*pi*(X(k)+time));
            effortsq(3*index1+1)=-2*pi*sin(2*pi*(X(k)+time));
            k=k+1;
        end
        effortsq(3 * N)=-2*pi*sin(2*pi*(X(N+1)+time));
    end
end
end