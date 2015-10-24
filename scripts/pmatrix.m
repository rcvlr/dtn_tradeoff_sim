function p = pmatrix(NN, L, R)

t = linspace(0,pi,NN);
m = NaN(NN);

for i=1:NN
    for j=1:NN
        p(j,i) = p_ji(t(j), t(i), t, L, R, m);
    end
end
