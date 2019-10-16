c=zeros(3,1);
t=zeros(3,1);
c(1)=5/9;
c(2)=8/9;
c(3)=5/9;
t(1)=-sqrt(3/5);
t(2)=0;
t(3)=sqrt(3/5);
switch (option)
    case 1
        dphid=zeros(3,3);
        for index1=1:3
            dphid(1,index1)= (2/h)*(t(index1)-1/2);
            dphid(2,index1)= -(4/h)*(t(index1));
            dphid(3,index1)= (2/h)*(t(index1)+1/2);
        end
        phip=zeros(2,3);
        for index1=1:3
            phip(1,index1)= (1/2)*(1-t(index1));
            phip(2,index1)= (1/2)*(1+t(index1));
        end
        Sp=zeros(3,2);
        for index1=1:3
            for index2=1:2
                for index3=1:3
                    Sp(index1,index2)=  Sp(index1,index2) +  h/2 * ( c(index3) * dphid(index1,index3)* phip(index2,index3));
                end
            end
        end
        phid=zeros(3,3);
        for index1=1:3
            phid(1,index1)= (t(index1)*(t(index1)-1))/2;
            phid(2,index1)= -((t(index1)+1)*(t(index1)-1));
            phid(3,index1)= (t(index1)*(t(index1)+1))/2;
        end
        Sq=zeros(2,3);
        dphip(1)=-1/h;
        dphip(2)=1/h;
        for index1=1:2
            for index2=1:3
                for index3=1:3
                    Sq(index1,index2)=  Sq(index1,index2) +  h/2 * ( c(index3) * phid(index2,index3)* dphip(index1));
                end
            end
        end
    case 2 
       dphid=zeros(4,3);
        for index1=1:3
            dphid(1,index1)= -(4/(3*h))*((t(index1)+1/2)*(t(index1)-1/2) + (t(index1)-1/2)*(t(index1)-1) +(t(index1)+1/2)*(t(index1)-1));
            dphid(2,index1)= (8/(3*h))*((t(index1)+1)*(t(index1)-1/2) + (t(index1)-1/2)*(t(index1)-1) +(t(index1)+1)*(t(index1)-1));
            dphid(3,index1)= -(8/(3*h))*((t(index1)+1)*(t(index1)+1/2) + (t(index1)+1/2)*(t(index1)-1) +(t(index1)+1)*(t(index1)-1));
            dphid(4,index1)= (4/(3*h))*((t(index1)+1)*(t(index1)+1/2) + (t(index1)-1/2)*(t(index1)+1) +(t(index1)+1/2)*(t(index1)-1/2));
        end
        phip=zeros(3,3);
        for index1=1:3
            phip(1,index1)= (t(index1)*(t(index1)-1))/2;
            phip(2,index1)= -((t(index1)+1)*(t(index1)-1));
            phip(3,index1)= (t(index1)*(t(index1)+1))/2;
        end
        Sp=zeros(4,3);
        for index1=1:4
            for index2=1:3
                for index3=1:3
                    Sp(index1,index2)=  Sp(index1,index2) +  h/2 * ( c(index3) * dphid(index1,index3)* phip(index2,index3));
                end
            end
        end
        phid=zeros(4,3);
         for index1=1:3
            phid(1,index1)= - (2/3)*((t(index1)+1/2)*(t(index1)-1/2)*(t(index1)-1));
            phid(2,index1)=   (4/3)*((t(index1)+1)*(t(index1)-1/2)*(t(index1)-1));
            phid(3,index1)= - (4/3)*((t(index1)+1)*(t(index1)+1/2)*(t(index1)-1));
            phid(4,index1)=   (2/3)*((t(index1)+1)*(t(index1)+1/2)*(t(index1)-1/2));
         end
        dphip=zeros(3,3);
        for index1=1:3
            dphip(1,index1)= (2/h)*(t(index1)-1/2);
            dphip(2,index1)= -(4/h)*(t(index1));
            dphip(3,index1)= (2/h)*(t(index1)+1/2);
        end
        Sq=zeros(3,4);
        for index1=1:3
            for index2=1:4
                for index3=1:3
                    Sq(index1,index2)=  Sq(index1,index2) +  h/2 * ( c(index3) * phid(index2,index3)* dphip(index1,index3));
                end
            end
        end
end