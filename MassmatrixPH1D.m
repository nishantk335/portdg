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
switch(option)
    case 1
        Mp=zeros(3,3);
        Mq=zeros(2,2);
        MMp=sparse(3*(N),3*(N));
        MMq=sparse(2*(N),2*(N));
        phi=zeros(3,4);
        for index1=1:4
            phi(1,index1)= (t(index1)*(t(index1)-1))/2;
            phi(2,index1)= -((t(index1)+1)*(t(index1)-1));
            phi(3,index1)= (t(index1)*(t(index1)+1))/2;
        end
        
        Mq(1,1)=h/3;
        Mq(1,2)=h/6;
        Mq(2,1)=h/6;
        Mq(2,2)=h/3;
        % loop through each node point
            for index1=1:3
                for index2=1:3
                    for index3=1:4
                        Mp(index1,index2)=  Mp(index1,index2) +  h/2 * ( c(index3)  * phi(index1,index3)* phi(index2,index3));
                    end
                end
            end
        dgelmati=1;
        jj=1;
        for i=1:N
            for ind1 = 1:3
                dgelmatj=jj;
                for ind2 = 1:3
                    MMp(dgelmati,dgelmatj)= MMp(dgelmati,dgelmatj) + Mp(ind1,ind2);
                    dgelmatj=dgelmatj+1;
                end
                dgelmati=dgelmati+1;
            end
            jj=jj+3;
        end
        
        dgelmati=1;
        jj=1;
        for i=1:N % loop through each node point
            for ind1 = 1:2
                dgelmatj=jj;
                for ind2 = 1:2
                    MMq(dgelmati,dgelmatj)= MMq(dgelmati,dgelmatj) + Mq(ind1,ind2);
                    dgelmatj=dgelmatj+1;
                end
                dgelmati=dgelmati+1;
            end
            jj=jj+2;
        end
    case 2
        Mp=zeros(4,4);
        Mq=zeros(3,3);
        MMp=sparse(4*(N),4*(N));
        MMq=sparse(3*(N),3*(N));
        phip=zeros(4,4);
        for index1=1:4
            phip(1,index1)= - (2/3)*((t(index1)+1/2)*(t(index1)-1/2)*(t(index1)-1));
            phip(2,index1)=   (4/3)*((t(index1)+1)*(t(index1)-1/2)*(t(index1)-1));
            phip(3,index1)= - (4/3)*((t(index1)+1)*(t(index1)+1/2)*(t(index1)-1));
            phip(4,index1)=   (2/3)*((t(index1)+1)*(t(index1)+1/2)*(t(index1)-1/2));
        end
        for index1=1:4
            for index2=1:4
                for index3=1:4
                    Mp(index1,index2)=  Mp(index1,index2) +  h/2 * ( c(index3) * phip(index1,index3)* phip(index2,index3));
                end
            end
        end
        phiq=zeros(3,4);
        for index1=1:4
            phiq(1,index1)= (t(index1)*(t(index1)-1))/2;
            phiq(2,index1)= -((t(index1)+1)*(t(index1)-1));
            phiq(3,index1)= (t(index1)*(t(index1)+1))/2;
        end
        for index1=1:3
            for index2=1:3
                for index3=1:4
                    Mq(index1,index2)=  Mq(index1,index2) +  h/2 * ( c(index3) * phiq(index1,index3)* phiq(index2,index3));
                end
            end
        end
        dgelmati=1;
        jj=1;
        for i=1:N % loop through each node point
            for ind1 = 1:4
                dgelmatj=jj;
                for ind2 = 1:4
                    MMp(dgelmati,dgelmatj)= MMp(dgelmati,dgelmatj) + Mp(ind1,ind2);
                    dgelmatj=dgelmatj+1;
                end
                dgelmati=dgelmati+1;
            end
            jj=jj+4;
        end
        dgelmati=1;
        jj=1;
        for i=1:N % loop through each node point
            for ind1 = 1:3
                dgelmatj=jj;
                for ind2 = 1:3
                    MMq(dgelmati,dgelmatj)= MMq(dgelmati,dgelmatj) + Mq(ind1,ind2);
                    dgelmatj=dgelmatj+1;
                end
                dgelmati=dgelmati+1;
            end
            jj=jj+3;
        end
        
end


        
