format long 
alpha=0;
DGPH1Dmesh;
Butchertable;
SSp=zeros(topologyd*N,topologyp*N);%The stiffness matrix for parameter p
SSq=zeros(topologyp*N,topologyd*N);%The stiffness matrix for parameter q
%--------------------------------------------- Assembly of the Stiffness Matrix at Element
%Diagonal Block
KpPH1D;
SSp = kron(eye(N),Sp);
SSq = kron(eye(N),Sq);
% % % % ---------------------------------------------Stiffness Matrix at face
% %  
%%%%Left Left Block
LLpPH1D;
dgelmati=1;
jj=1;
for i=1:N-1% loop through each node point
    for ind1 = 1:topologyd
        dgelmatj=jj;
        for ind2 = 1:topologyp
            SSp(dgelmati,dgelmatj)= SSp(dgelmati,dgelmatj) -  (1-alpha) * LLp(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end 
        dgelmati=dgelmati+1;
    end
    jj=jj+topologyp;
end
dgelmati=1 + (N-1) * topologyd;
jj=1 + (N-1) * topologyp;
for ind1 = 1:topologyd
    dgelmatj=jj;
    for ind2 = 1:topologyp
        SSp(dgelmati,dgelmatj)= SSp(dgelmati,dgelmatj) -  LLp(ind1,ind2);
        dgelmatj=dgelmatj+1;
    end
    dgelmati=dgelmati+1;
end
dgelmati=1;
jj=1;


for i=1:N-1% loop through each node point
    for ind1 = 1:topologyp
        dgelmatj=jj;
        for ind2 = 1:topologyd
            SSq(dgelmati,dgelmatj)= SSq(dgelmati,dgelmatj) -  (alpha) * LLq(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end 
        dgelmati=dgelmati+1;
    end
    jj=jj+topologyd;
end


% % % % % % % % % % % % % Right-Right Block
RRpPH1D;
dgelmati=1+topologyd;
 jj=1+topologyp;
for i=2:N % loop through each node point
    for ind1 = 1:topologyd
        dgelmatj=jj;
        for ind2 = 1:topologyp
            SSp(dgelmati,dgelmatj)= SSp(dgelmati,dgelmatj) + alpha * RRp(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end
        dgelmati=dgelmati+1;
    end
    jj=jj+topologyp;
end
dgelmati=1;
jj=1;
% loop through each node point
for ind1 = 1:topologyp
    dgelmatj=jj;
    for ind2 = 1:topologyd
        SSq(dgelmati,dgelmatj)= SSq(dgelmati,dgelmatj) +  RRq(ind1,ind2);
        dgelmatj=dgelmatj+1;
    end
    dgelmati=dgelmati+1;
end

dgelmati=1+topologyp;
 jj=1+topologyd;
for i=2:N % loop through each node point
    for ind1 = 1:topologyp
        dgelmatj=jj;
        for ind2 = 1:topologyd
            SSq(dgelmati,dgelmatj)= SSq(dgelmati,dgelmatj) +  (1-alpha)* RRq(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end
        dgelmati=dgelmati+1;
    end
    jj=jj+topologyd;
end
% % % % % % % % % % % % % Left Right Block
LRpPH1D;%Left-Right Block matrix for each internal face
dgelmati=1;
jj=1+topologyd;
for i=1:N-1 % for face elements
    for ind1 = 1:topologyp
        dgelmatj=jj;
        for ind2 = 1:topologyd
            SSq(dgelmati,dgelmatj)= SSq(dgelmati,dgelmatj)- (1-alpha) * LRq(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end
        dgelmati=dgelmati+1;
    end 
    jj=jj+topologyd;
end 
dgelmati=1;
jj=1+topologyp;
for i=1:N-1 % for face elements
    for ind1 = 1:topologyd
        dgelmatj=jj;
        for ind2 = 1:topologyp
            SSp(dgelmati,dgelmatj)= SSp(dgelmati,dgelmatj)- alpha * LRp(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end
        dgelmati=dgelmati+1;
    end 
    jj=jj+topologyp;
end
% % % % % % % % % % % % 
% % % % % % % % % % Right-Left Block
RLpPH1D;%Right-Left Block Matrix for each internal face 
dgelmati=1+topologyd;
jj=1;
for i=1:N-1 % for face elements
    for ind1 = 1:topologyd
        dgelmatj=jj;
        for ind2 = 1:topologyp
            SSp(dgelmati,dgelmatj)= SSp(dgelmati,dgelmatj)+ (1-alpha)* RLp(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end
        dgelmati=dgelmati+1;
    end 
    jj=jj+topologyp;
end 
dgelmati=1+topologyp;
jj=1;
for i=1:N-1 % for face elements
    for ind1 = 1:topologyp
        dgelmatj=jj;
        for ind2 = 1:topologyd
            SSq(dgelmati,dgelmatj)= SSq(dgelmati,dgelmatj) + (alpha)* RLq(ind1,ind2);
            dgelmatj=dgelmatj+1;
        end
        dgelmati=dgelmati+1;
    end 
    jj=jj+topologyd;
end
% % % %
starttime=time(1);
k1p=zeros(topologyd * N ,1);
k2p=zeros(topologyd * N ,1);
k3p=zeros(topologyd * N ,1);
k4p=zeros(topologyd * N ,1);
k1q=zeros(topologyp * N ,1);
k2q=zeros(topologyp * N ,1);
k3q=zeros(topologyp * N ,1);
k4q=zeros(topologyp * N ,1);
[eff0p,eff0q]=PH1Dinputs(Y,mu,starttime,X,X1,X2,N,Xdual,option);
neweff0p=eff0p;
neweff0q=eff0q;
MassmatrixPH1D;
[L1,U1] = lu(MMp);
[L2,U2] = lu(MMq);
MMp1=inv(U1) *inv(L1);
MMq1=inv(U2) *inv(L2);
for i=1:length(time)-1
    newtime=time(i);
    newS1p = SSp * neweff0q;
    newS1q = SSq * neweff0p;
    gq1= -2 * pi* sin(2*pi*(X(1)+newtime));
    gp2=-2 * pi* sin(2*pi*(X(N+1)+newtime));
    newS1p(1) = newS1p(1)+  gq1 * 1;
    newS1q(topologyp*N)=newS1q(topologyp*N) -   gp2 * 1;
    [Sourcetermp,Sourcetermq]=SourceTerm(newtime,X,Xdual,X1,X2,MMp,Sp,Sq,MMq,N,option);
    newS1p=newS1p + Sourcetermp;
    newS1q=newS1q + Sourcetermq;
   k1p=MMp1 * newS1p;
    k1q=MMq1 * newS1q;
    neweff0p1=  neweff0p + (dt * a(2,1) *k1p) ;
    neweff0q1=  neweff0q + (dt * a(2,1) *k1q);
    newtime1 = newtime +  cc(2) * dt ;
% % % % % % %     ------------------------------------------ k2 calculation
    newS1p = SSp * neweff0q1;
    newS1q = SSq * neweff0p1;
    gq1= -2 * pi* sin(2*pi*(X(1)+newtime1));
%     gq2= -2 * pi* sin(2*pi*(X(N+1)+newtime1));
%     gp1=-2 * pi* sin(2*pi*(X(1)+newtime1));
    gp2=-2 * pi* sin(2*pi*(X(N+1)+newtime1));
% gq1=X(1);
% 
%     gp2=1;
    newS1p(1) = newS1p(1)+   gq1 * 1;
%     newS1p(topologyd * N) =  newS1p(topologyd * N) -  gq2 * 1;
%     newS1q(1) = newS1q(1)+  gp1 * 1;
    newS1q(topologyp*N)=newS1q(topologyp*N) -    gp2 * 1;
    [Sourcetermp,Sourcetermq]=SourceTerm(newtime1,X,Xdual,X1,X2,MMp,Sp,Sq,MMq,N,option);
    newS1p=newS1p + Sourcetermp;
    newS1q=newS1q + Sourcetermq;
    k2p=MMp1 * newS1p;
    k2q=MMq1 * newS1q;
    neweff0p2 = neweff0p + (dt * (a(3,1)+a(3,2)) * k2p);
    neweff0q2 = neweff0q + (dt * (a(3,1)+a(3,2)) * k2q);
    newtime2 = newtime +  cc(3) * dt ;
% % % % % %     % % % % %      %------------------------------------------ k3 calculation
    newS1p = SSp * neweff0q2;
    newS1q = SSq * neweff0p2;
    gq1= -2 * pi* sin(2*pi*(X(1)+newtime2));
%     gq2= -2 * pi* sin(2*pi*(X(N+1)+newtime2));
%     gp1=-2 * pi* sin(2*pi*(X(1)+newtime2));
    gp2=-2 * pi* sin(2*pi*(X(N+1)+newtime2));

% gq1=X(1);
%  
%     gp2=1;
    newS1p(1) = newS1p(1)+   gq1 * 1;
%     newS1p(topologyd * N) =  newS1p(topologyd * N) -  gq2 * 1;
%     newS1q(1) = newS1q(1)+  gp1 * 1;
    newS1q(topologyp*N)=newS1q(topologyp*N) -   gp2 * 1;
    [Sourcetermp,Sourcetermq]=SourceTerm(newtime2,X,Xdual,X1,X2,MMp,Sp,Sq,MMq,N,option);
    newS1p=newS1p + Sourcetermp;
    newS1q=newS1q + Sourcetermq;
   k3p=MMp1 * newS1p;
    k3q=MMq1 * newS1q;
    neweff0p3 = neweff0p + (dt *(a(4,1)+a(4,2)+a(4,3)) * k3p);
    neweff0q3 = neweff0q + (dt *(a(4,1)+a(4,2)+a(4,3)) * k3q);
    newtime3 = newtime +  cc(4) * dt ;
% % % % %     % % % % %     %------------------------------------------ k4 calculation
    newS1p = SSp * neweff0q3;
    newS1q = SSq * neweff0p3;
    gq1= -2 * pi* sin(2*pi*(X(1)+newtime3));
%     gq2= -2 * pi* sin(2*pi*(X(N+1)+newtime3));
%     gp1=-2 * pi* sin(2*pi*(X(1)+newtime3));
    gp2=-2 * pi* sin(2*pi*(X(N+1)+newtime3));
% gq1=X(1);
%     gp2=1;
    newS1p(1) = newS1p(1)+   gq1 * 1;
%     newS1p(topologyd * N) =  newS1p(topologyd * N) -  gq2 * 1;
%     newS1q(1) = newS1q(1)+  gp1 * 1;
    newS1q(topologyp*N)=newS1q(topologyp*N) -    gp2 * 1;
    [Sourcetermp,Sourcetermq]=SourceTerm(newtime3,X,Xdual,X1,X2,MMp,Sp,Sq,MMq,N,option);
    newS1p=newS1p + Sourcetermp;
    newS1q=newS1q + Sourcetermq;
    k4p=MMp1 * newS1p;
    k4q=MMq1 * newS1q;
    neweff0p=neweff0p + dt * ( b(1) * k1p + b(2) * k2p + b(3) * k3p + b(4) * k4p );
    neweff0q=neweff0q + dt * ( b(1) * k1q + b(2) * k2q + b(3) * k3q + b(4) * k4q );
end
answerp=neweff0p;
answerq=neweff0q;
if option==1
    proanswerp= zeros(2*N+1,1);
    proanswerp(1)=answerp(1);
    k=1;
    for i=2:2:2*N+1
        proanswerp(i) = answerp(3*k-1);
        k=k+1;
    end
    k=1;
    for i=3:2:2*N
        proanswerp(i) = (answerp(3*k) + answerp(3*k+1))/2;
        k=k+1;
    end
    proanswerp(2*N+1)=answerp(3*N);
    actualanswerp=zeros(2*N+1,1);
    finaltime=1;
    for i=1:2*N+1
         actualanswerp(i)=- 2 * pi * sin(2*pi*(Xquad(i)+finaltime));
% actualanswerp=1;
    end
    proanswerq= zeros(N+1,1);
    proanswerq(1) = answerq(1);
    k=2;
    for i=2:N;
        proanswerq(i) = (answerq(k) + answerq(k+1))/2;
        k=k+2;
    end
    proanswerq(N+1) = answerq(2*N);
    actualanswerq= zeros(N+1,1);
    for i=1:N+1
          actualanswerq(i)=  -2 * pi * sin(2*pi*(X(i)+finaltime));
% actualanswerq(i)=X(i);
    end
else if option==2
        proanswerp= zeros(3*N+1,1);
        proanswerp(1)=answerp(1);
        k=1;
        for i=2:3:3*N+1
            proanswerp(i) = answerp(4*k-2);
            proanswerp(i+1)=answerp(4*k-1);
            k=k+1;
        end
        k=1;
        for i=4:3:3*N
            proanswerp(i) = (answerp(4*k) + answerp(4*k+1))/2;
            k=k+1;
        end
        proanswerp(3*N+1)=answerp(4*N);
        actualanswerp=zeros(3*N+1,1);
        finaltime=1;
        for i=1:3*N+1
            actualanswerp(i) = -2 * pi * sin(2*pi*(Xcub(i)+finaltime)) ;
% actualanswerp(i)=finaltime^2;
        end
        proanswerq= zeros(2*N+1,1);
        proanswerq(1)=answerq(1);
        k=1;
        for i=2:2:2*N+1
            proanswerq(i) = answerq(3*k-1);
            k=k+1;
        end
        k=1;
        for i=3:2:2*N
            proanswerq(i) = (answerq(3*k) + answerq(3*k+1))/2;
            k=k+1;
        end
        proanswerq(2*N+1)=answerq(3*N);
        actualanswerq= zeros(2*N+1,1);
        for i=1:2*N+1
            actualanswerq(i)=- 2 * pi * sin(2*pi*(Xquad(i)+finaltime));
        end
    end
end
[L2errorp,L2errorq,Linferrorp,Linferrorq]=errorp(answerp,answerq,N,X,h,option);
L2p(nn)=L2errorp;
L2q(nn)=L2errorq;
Linfp(nn)=Linferrorp;
Linfq(nn)=Linferrorq;