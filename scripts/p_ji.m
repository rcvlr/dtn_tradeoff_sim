# this function calculates the probability to go from direction 
# t_j to t_i in multihop with a static scenario, i.e., nodes do not move.
# t_i and t_j must be positive. If you need to calculate for negative
# values, gives the abs value, due to symmetry the result will be correct.
#
# Input:    t_j, t_i: directions
#           t, usually is t = linspace(0, pi, N). N is the number of
#           discretization interval of the positive axis. So if you have 
#           100 values from -pi to pi, N will be 50.
#           L: density, R: tx range.
#           m: for internal use.
function [p m] = p_ji(t_j, t_i, t, L, R, m)

if (abs(t_i) >= abs(t_j))
    p = 0;
    return
end

N = numel(t);
p_i = @(t_i) (1-exp(-L*pi*R^2/(2*N-2)))*exp(-L*R^2*abs(t_i));

% get i and j
for i=1:N
    if (t(i) == t_i)
        break
    end
end    
for j=1:N
    if (t(j) == t_j)
        break
    end
end

% recursive part
pint = 0;
for k=(i+1):(j-1)
    if (isnan(m(j,k)))
        [a m] = p_ji(t(j),t(k),t, L, R, m);
        m(j,k) = a;
    else
        a = m(j,k);
    end
    
    if (isnan(m(k,i)))
        [b m] = p_ji(t(k),t(i),t, L, R, m);
        m(k,i) = b;
    else
        b = m(k,i);
    end
    
    pint += (a*b);    
end

p = p_i(t_i) + 2*pint;
