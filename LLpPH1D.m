switch(option)
    case 1
        LLp=zeros(3,2);
        LLp(3,2)=1;
        LLq=zeros(2,3);
        LLq(2,3)=1;
    case 2
        LLp=zeros(4,3);
        LLp(4,3)=1;
        LLq=zeros(3,4);
        LLq(3,4)=1;
end