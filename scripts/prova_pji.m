clear all;
R = 10;
NN = 50;
L = 0.1
t = linspace(0,pi,NN);
m = NaN(NN);

for i=1:NN
    for j=1:NN
        p(j,i) = p_ji(t(j),t(i),linspace(0,pi,NN),L,R,m);
    end
    disp(i);
    fflush(stdout);
end
