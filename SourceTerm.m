function [Sourcetermp,Sourcetermq] = SourceTerm(time,X,Xdual,X1,X2,MMp,Sp,Sq,MMq,N,option)
if option==1
    %     --------------------------------------------SourceTermp
    fp= zeros(3 * N ,1);
    Sourcetermp= zeros(3 * N ,1);
    fp(1)= - 4 * pi * pi * cos(2*pi*(X(1)+time));
    k1=1;
    for index1=1:N
        fp(3*index1-1)= - 4 * pi * pi * cos(2*pi*(Xdual(k1)+time));
        k1=k1+1;
    end
    k2=2;
    for index1=1:N-1
        fp(3*index1)=  -4 * pi * pi * cos(2*pi*(X(k2)+time));
        fp(3*index1+1)=  -4 * pi * pi * cos(2*pi*(X(k2)+time));
        k2=k2+1;
    end
    fp(3 * N)=  -4 * pi * pi * cos(2*pi*(X(N+1)+time));
    
    fq=zeros(2*N,1);
    fq(1)=  -2 * pi* sin(2*pi*(X(1)+time));
    k=2;
    for index1=2:2:(2*N - 1)
        fq(index1)=  -2 * pi* sin(2*pi*(X(k)+time));
        fq(index1+1)= -2 * pi* sin(2*pi*(X(k)+time));
        k=k+1;
    end
    fq(2*N)= -2 * pi* sin(2*pi*(X(N+1)+time));
    for i=1:N
        k=1+3*(i-1);
        k1=1+2*(i-1);
        Sourcetermp(k:k+2)=Sourcetermp(k:k+2)+ Sq' * fq(k1:k1+1);
    end
    Sourcetermp = Sourcetermp + MMp * fp;
    %         --------------------------------------------SourceTermq
    fq=zeros(2*N,1);
    Sourcetermq=zeros(2*N,1);
    fq(1)=  -4 * pi * pi * cos(2*pi*(X(1)+time));
    k=2;
    for index1=2:2:(2*N - 1)
        fq(index1)=  -4 * pi * pi * cos(2*pi*(X(k)+time));
        fq(index1+1)= - 4 * pi * pi *cos(2*pi*(X(k)+time));
        k=k+1;
    end
    fq(2*N)= - 4 * pi * pi * cos(2*pi*(X(N+1)+time));
    fp= zeros(3 * N ,1);
    fp(1)= - 2 * pi * sin(2*pi*(X(1)+time));
    k1=1;
    for index1=1:N
        fp(3*index1-1)= - 2 * pi * sin(2*pi*(Xdual(k1)+time));
        k1=k1+1;
    end
    k2=2;
    for index1=1:N-1
        fp(3*index1)=  - 2 * pi * sin(2*pi*(X(k2)+time));
        fp(3*index1+1)= - 2 * pi * sin(2*pi*(X(k2)+time));
        k2=k2+1;
    end
    fp(3 * N)=  - 2 * pi * sin(2*pi*(X(N+1)+time));
    for i=1:N
        k=1+3*(i-1);
        k1=1+2*(i-1);
        Sourcetermq(k1:k1+1)=Sourcetermq(k1:k1+1) + Sp' * fp(k:k+2);
    end
    Sourcetermq=Sourcetermq + MMq * fq;
%--------------------------------------------Case 2
else if option==2
    %--------------------------------------------SourceTermp
    fp= zeros(4 * N ,1);
    Sourcetermp=zeros(4 * N ,1);
    fp(1)=- 4 * pi * pi * cos(2*pi*(X(1)+time));
    k1=1;
    for index1=1:N
        fp(4*index1-2)=- 4 * pi * pi * cos(2*pi*(X1(k1)+time));
        fp(4*index1-1)=- 4 * pi * pi * cos(2*pi*(X2(k1)+time));
        k1=k1+1;
    end
    k2=2;
    for index1=1:N-1
        fp(4*index1)=- 4 * pi * pi * cos(2*pi*(X(k2)+time));
        fp(4*index1+1)=- 4 * pi * pi* cos(2*pi*(X(k2)+time));
        k2=k2+1;
    end
    fp(4 * N)=- 4 * pi * pi * cos(2*pi*(X(N+1)+time)) ;
    fq= zeros(3 * N ,1);
    fq(1)=- 2 * pi * sin(2*pi*(X(1)+time));
    k1=1;
    for index1=1:N
        fq(3*index1-1)=- 2 * pi * sin(2*pi*(Xdual(k1)+time));
        k1=k1+1;
    end
    k2=2;
    for index1=1:N-1
        fq(3*index1)= - 2 * pi * sin(2*pi*(X(k2)+time));
        fq(3*index1+1)= - 2 * pi * sin(2*pi*(X(k2)+time));
        k2=k2+1;
    end
    fq(3 * N)= - 2 * pi * sin(2*pi*(X(N+1)+time)) ;
    for i=1:N
        k=1+4*(i-1);
        k1=1+3*(i-1);
        Sourcetermp(k:k+3)=Sourcetermp(k:k+3)+ Sq' * fq(k1:k1+2);
    end
    Sourcetermp = Sourcetermp +  MMp * fp;
    %         --------------------------------------------SourceTermq
    fq= zeros(3 * N ,1);
    Sourcetermq =zeros(3 * N,1);
    fq(1)=- 4 * pi * pi * cos(2*pi*(X(1)+time));
    k1=1;
    for index1=1:N
        fq(3*index1-1)=- 4 * pi * pi * cos(2*pi*(Xdual(k1)+time));
        k1=k1+1;
    end
    k2=2;
    for index1=1:N-1
        fq(3*index1)= - 4 * pi * pi * cos(2*pi*(X(k2)+time));
        fq(3*index1+1)= - 4 * pi * pi* cos(2*pi*(X(k2)+time));
        k2=k2+1;
    end
    fq(3 * N)= - 4 * pi * pi * cos(2*pi*(X(N+1)+time)) ;
    fp= zeros(4 * N ,1);
    fp(1)=- 2 * pi * sin(2*pi*(X(1)+time));
    k1=1;
    for index1=1:N
        fp(4*index1-2)=-2 * pi * sin(2*pi*(X1(k1)+time));
        fp(4*index1-1)=-2 * pi * sin(2*pi*(X2(k1)+time));
        k1=k1+1;
    end
    k2=2;
    for index1=1:N-1
        fp(4*index1)=-2 * pi * sin(2*pi*(X(k2)+time));
        fp(4*index1+1)=-2 * pi * sin(2*pi*(X(k2)+time));
        k2=k2+1;
    end
    fp(4 * N)=-2 * pi * sin(2*pi*(X(N+1)+time)) ;
    for i=1:N
        k=1+4*(i-1);
        k1=1+3*(i-1);
        Sourcetermq(k1:k1+2)=Sourcetermq(k1:k1+2) + Sp' * fp(k:k+3);
    end
    Sourcetermq = Sourcetermq + MMq * fq;
    end
end
end