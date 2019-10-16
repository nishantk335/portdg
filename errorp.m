function [L2errorp,L2errorq,Linferrorp,Linferrorq]=errorp(answerp,answerq,N,X,h,option)
if option==1
    finaltime=1;
    jj=1;
    jj1=1;
    sump=0;
    sumq=0;
    maxerrorp=0;
    maxerrorq=0;
    for i=1:N
        hh= (X(i)+X(i+1))/2;
        c=zeros(4,1);
        t=zeros(4,1);
        c(1)=0.3478548451374538;
        c(2)=0.6521451548625461;
        c(3)=0.6521451548625461;
        c(4)=0.3478548451374538;
        t(1)=-0.8611363115940526;
        t(2)=-0.3399810435848563;
        t(3)=0.3399810435848563;
        t(4)=0.8611363115940526;
        phid=zeros(3,4);
        for index1=1:4
            phid(1,index1)= (t(index1)*(t(index1)-1))/2;
            phid(2,index1)= -((t(index1)+1)*(t(index1)-1));
            phid(3,index1)= (t(index1)*(t(index1)+1))/2;
        end
        phip=zeros(2,4);
        for index1=1:4
            phip(1,index1)= (1/2)*(1-t(index1));
            phip(2,index1)= (1/2)*(1+t(index1));
        end
        for j=1:4
            % the solution at each gaussian point
            sum1p(j)=answerp(jj)*phid(1,j)+answerp(jj+1)*phid(2,j)+answerp(jj+2)*phid(3,j);%numerical solution
            sum2p(j)=-2*pi*sin(2*pi*(((h/2)*t(j) +hh)+finaltime));%exact Solution
        end
        errorp=abs(sum2p-sum1p);
        errorpp=max(errorp);
        if(errorpp>maxerrorp)
            maxerrorp=errorpp;
        end
        for j=1:4
            % the solution at each gaussian point
            sum1q(j)=answerq(jj1)*phip(1,j)+answerq(jj1+1)*phip(2,j);%numerical solution
            sum2q(j)=-2*pi*sin(2*pi*(((h/2)*t(j) +hh)+finaltime));%exact Solution

        end
        errorq=abs(sum2q-sum1q);
        errorqq=max(errorq);
        if(errorqq>maxerrorq)
            maxerrorq=errorqq;
        end
        elerrorp(i)=(h/2)*((c(1)*(sum2p(1)-sum1p(1))^2)+(c(2)*(sum2p(2)-sum1p(2))^2)+(c(3)*(sum2p(3)-sum1p(3))^2)+(c(4)*(sum2p(4)-sum1p(4))^2));
        elerrorq(i)=(h/2)*((c(1)*(sum2q(1)-sum1q(1))^2)+(c(2)*(sum2q(2)-sum1q(2))^2)+(c(3)*(sum2q(3)-sum1q(3))^2)+(c(4)*(sum2q(4)-sum1q(4))^2));
        sump=sump+abs(elerrorp(i));
        sumq=sumq+abs(elerrorq(i));
        jj=jj+3;
        jj1=jj1+2;
    end
    Linferrorp=maxerrorp;
    Linferrorq=maxerrorq;
    L2errorp=sqrt(sump);
    L2errorq=sqrt(sumq);
else if option==2
        finaltime=1;
        jj=1;
        jj1=1;
        sump=0;
        sumq=0;
         maxerrorp=0;
    maxerrorq=0;
        for i=1:N
            hh= (X(i)+X(i+1))/2;
            c=zeros(4,1);
            t=zeros(4,1);
            c(1)=0.3478548451374538;
            c(2)=0.6521451548625461;
            c(3)=0.6521451548625461;
            c(4)=0.3478548451374538;
            t(1)=-0.8611363115940526;
            t(2)=-0.3399810435848563;
            t(3)=0.3399810435848563;
            t(4)=0.8611363115940526;
            phid=zeros(4,4);
            for index1=1:4
                phid(1,index1)= - (2/3)*((t(index1)+1/2)*(t(index1)-1/2)*(t(index1)-1));
                phid(2,index1)=   (4/3)*((t(index1)+1)*(t(index1)-1/2)*(t(index1)-1));
                phid(3,index1)= - (4/3)*((t(index1)+1)*(t(index1)+1/2)*(t(index1)-1));
                phid(4,index1)=   (2/3)*((t(index1)+1)*(t(index1)+1/2)*(t(index1)-1/2));
            end
            phip=zeros(3,4);
            for index1=1:4
                phip(1,index1)= (t(index1)*(t(index1)-1))/2;
                phip(2,index1)= -((t(index1)+1)*(t(index1)-1));
                phip(3,index1)= (t(index1)*(t(index1)+1))/2;
            end
            for j=1:4
                % the solution at each gaussian point
                sum1p(j)=answerp(jj)*phid(1,j)+answerp(jj+1)*phid(2,j)+answerp(jj+2)*phid(3,j)+answerp(jj+3)*phid(4,j);%numerical solution
                sum2p(j)=-2*pi*sin(2*pi*((h/2*t(j)+hh)+finaltime));%exact Solution
            end
            errorp=abs(sum2p-sum1p);
            errorpp=max(errorp);
            if(errorpp>maxerrorp)
                maxerrorp=errorpp;
            end
            for j=1:4
                % the solution at each gaussian point
                sum1q(j)=answerq(jj1)*phip(1,j)+answerq(jj1+1)*phip(2,j)+answerq(jj1+2)*phip(3,j);%numerical solution
                sum2q(j)=-2*pi*sin(2*pi*((h/2*t(j)+hh)+finaltime));%exact Solution
            end
            errorq=abs(sum2q-sum1q);
            errorqq=max(errorq);
            if(errorqq>maxerrorq)
                maxerrorq=errorqq;
            end
            elerrorp(i)=(h/2)*((c(1)*(sum2p(1)-sum1p(1))^2)+(c(2)*(sum2p(2)-sum1p(2))^2)+(c(3)*(sum2p(3)-sum1p(3))^2)+(c(4)*(sum2p(4)-sum1p(4))^2));
            elerrorq(i)=(h/2)*((c(1)*(sum2q(1)-sum1q(1))^2)+(c(2)*(sum2q(2)-sum1q(2))^2)+(c(3)*(sum2q(3)-sum1q(3))^2)+(c(4)*(sum2q(4)-sum1q(4))^2));
            sump=sump+abs(elerrorp(i));
            sumq=sumq+abs(elerrorq(i));
            jj=jj+4;
            jj1=jj1+3;
        end
          Linferrorp=maxerrorp;
    Linferrorq=maxerrorq;
        L2errorp=sqrt(sump);
        L2errorq=sqrt(sumq);
    end
end