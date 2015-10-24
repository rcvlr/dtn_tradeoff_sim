function [f q t]=limprob0(R0, V0, L, R, N)
%~ clear all;

%~ R0 = 1;        % turning rate
%~ V0 = 1;        % speed
%~ L  = 0.1;      % density
%~ R  = 1;        % tx range
%~ N  = 100;      % discretization

beta  = 2*V0*L*R/pi;

rgamma = @(t_i,t_j) beta*(abs(sin((t_i-t_j)/2)));

t = linspace(-pi,pi,N);   % state space
Q = zeros(size(t));       % Q-matrix

for i=1:N         % row
  for j=1:N       % col
    if (j == i) continue end
    if ((abs(t(j)) >= abs(t(i))) || ((i==j) || (j == (numel(t)-i+1))))
      Q(i,j) = R0/N;
    else
      Q(i,j) = R0/N + 2*pi/N*rgamma(t(i),t(j));
    end
  end
end

% diagonal entries
Q = Q + diag(-sum(Q,2));

[U,S,V] = svd(Q);
f = U(:,N)/sum(U(:,N));

% check that the sum of the pdf is 1
if (abs(sum(f)-1) > 0.001)
	disp('error: integral of pdf is not 1')
end

%~ figure,
%~ plot(t,f,'-o');
%~ stem(t,f)
%~ xlim([-pi pi])
%~ ylim([0 1.2*max(f)])
%~ title(S(N,N))
%~ legend('f(\theta)')

q = -diag(Q);
t = t';
