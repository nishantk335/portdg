format long
x1=0;
x2=1;
prompt='Enter the number of elements';
N=input(prompt);
prompt='Enter the number of refinement';
nn=input(prompt);
X=linspace(0,1,N+1);
h=X(2)-X(1);
time=linspace(0,1,10001);
dt=time(2)-time(1);
Y=1;
mu=1;
for i=1:N
    Xdual(i)=X(i) + (X(i+1) - X(i))/2;
end
for i=1:N
    X1(i)=X(i) + (X(i+1) - X(i))/4;
    X2(i)=X(i) + (3*(X(i+1) - X(i)))/4;
end
Xquad=zeros(2*N+1,1);
Xcub= zeros(3*N+1,1);
k=1;
for i=1:2:2*N+1
    Xquad(i)=X(k);
    k=k+1;
end
k=1;
for i=2:2:2*N
    Xquad(i)=Xdual(k);
    k=k+1;
end
k=1;
for i=1:3:3*N+1
    Xcub(i)=X(k);
    k=k+1;
end
k=1;
for i=2:3:3*N
    Xcub(i)=X1(k);
    Xcub(i+1)=X2(k);
    k=k+1;
end
k=1;
for i=2:length(X)
    Xtot(k) = (X(i) + X(i-1))/2;
    k=k+1;
end
prompt='Enter 1 for linear\quad and 2 for quad\cub';
option=input(prompt);
if option==1
    topologyd=3;
    topologyp=2;
else if option==2
        topologyd=4;
        topologyp=3;
    end
end