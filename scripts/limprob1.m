function [f q t]=limprob1(R0, V0, L, R, N, Pmat)
%~ clear all;

%~ R0 = 1;        % turning rate
%~ V0 = 1;        % speed
%~ L  = 0.1;      % density
%~ R  = 1;        % tx range
%~ N  = 100;      % discretization

beta  = 2*V0*L*R/pi;

rgamma = @(t_i,t_j) (beta*(abs(sin((t_i-t_j)/2))));
pnot   = @(t_i) (exp(-L*R^2*abs(t_i)));

t = linspace(-pi,pi,N);   % state space
Q = zeros(size(t));       % Q-matrix

tt = t(ceil(N/2):end);
for j=1:N         % row
    % get jj
    for jj=1:numel(tt)
        if ((abs(t(j))-tt(jj)) <= 2*eps)
            break
        end
    end 
    
    %~ printf('t_j = %f, jj = %d\n',t(j),jj);
    for i=1:N       % col
        if (i == j)
            %~ printf('t_i = %f\n',t(i));
            %~ input('press any key to proceede');
            %~ printf('-------------------------\n');
            continue 
        end
        
        % get ii
        for ii=1:numel(tt)
            if ((abs(t(i))-tt(ii)) <= 2*eps)
                break
            end
        end
        
        %~ if ((abs(t(i)) >= abs(t(j))))
        if (ii > jj)
            %~ printf('t_i = %f\n',t(i));          
            a = 0;
            for k=(ii+1):numel(tt) a+=Pmat(k,ii); end
            %~ Q(j,i) = R0/(N-1)*(pnot(tt(ii)-tt(jj)) + exp(L*R^2*tt(jj))*a);
            Q(j,i) = R0/(N-1)*(pnot(tt(ii)-tt(jj)) + 2.8*a);
            
        else
            a = R0/(N-1);
            
            b1 = 0;
            %~ printf('t_i = %f, ii = %d\n',t(i),ii);
            for k=(ii+1):(jj-1)
                %~ printf('P(%f --> %f) = P(%d --> %d)\n',tt(k),t(i), k,ii);
                b1 += Pmat(k,ii);
            end
            b = (L*pi*R^2*R0/(N-1)) * (pnot(t(i)) + b1);
             
            c1 = 0;
            for k=(ii+1):(jj-1)
                c1 += rgamma(t(k),t(j))*Pmat(k,ii);
            end
            c = (rgamma(t(i),t(j))*pnot(t(i)) + c1)*2*pi/(N-1);
            
            Q(j,i) = a+b+c;
            
        end
        %~ input('press any key to proceede');
        %~ printf('-------------------------\n');
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
