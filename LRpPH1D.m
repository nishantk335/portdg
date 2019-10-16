switch(option)
    case 1
        LRp=zeros(3,2);
        LRp(3,1)=1;
        LRq=zeros(2,3);
        LRq(2,1)=1;
    case 2
        LRp=zeros(4,3);
        LRp(4,1)=1;
        LRq=zeros(3,4);
        LRq(3,1)=1;
end
