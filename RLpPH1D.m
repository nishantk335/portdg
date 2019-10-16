switch(option)
    case 1
        RLp=zeros(3,2);
        RLp(1,2)=1;
        RLq=zeros(2,3);
        RLq(1,3)=1;
    case 2
        RLp=zeros(4,3);
        RLp(1,3)=1;
        RLq=zeros(3,4);
        RLq(1,4)=1;
end
