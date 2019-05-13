%function metrics=metrics(node_pars,FR_type,potential_type,cost_type,...
%                                           density_type,num_pars,verbose)
% FILE 8
% 
% This function calculates the performance metrics following closely the
% derivations in the paper. The calculations will make more sense if you 
% read the paper. 
% 
% Input Parameters
% ================
% node_pars=[lambda vo r0] is a vector of parameters describing the nodes,
%    as follows:
%      * lambda: the density of nodes.
%      * v0: the speed of nodes.
%      * r0: the rate with which nodes change their direction.
% FR_type: the type of FR used. see the boundary.m routine for details.
% potential_type: the type of potential used. See the potential routine for
%     details.
% cost_type: the type of cost used. See the cost routine for details. 
% density_type: the type of direction density type used. See the
%     direction_density routine for details. 
% num_pars=[N B L S] is a vector of parameters specifying the accuracy of 
%     numerical calculations, as follows:
%       * N is the number of angle samples in the range [-pi,pi). The interval 
%               [-pi,pi) is broken in N subintervals and each angle is
%                at the center of one of these. 
%       * B is such that the FR fully belongs in the box [-B B -B B]
%       * L is the number of samples of both x and y in the range [-B B]
%       * S is the number of parts in which we break the boundary (first used in Part 11)
% verbose: a parameters that specifies how much output the routine produces.
%    * verbose=0: Nothing is plotted/printed
%    * verbose=1: Results are printed but not plotted
%    * verbose=2: Results are printed and plotted and verifications are made
%
% Output Parameters
% =================
% metrics: A structure containing rates.
%     metrics.thetas(N,1): the angles used in the angle discretization
%     metrics.locations(M,4): the locations used in the location discretization, in the
%         following format:
%         * (locations(j,1),locations(j,2)) are the cartesian coordinates
%         * (locations(j,3),locations(j,4)) are the polar coordinates
%     metrics.EN: E(N;theta,r) (see paper)
%     metrics.PE: P_E(theta,r) (see paper)
%     metrics.rA: rA(theta,theta')
%     metrics.rAtotal: rA(theta)
%     metrics.rB: rB(theta,theta',r')
%     metrics.rBtotal: rB(theta)
%     metrics.rC: rC(theta,theta',r')
%     metrics.rCtotal: rC(theta)
%     metrics.rD: rD(theta,theta',s)
%     metrics.rDtotal: rD(theta)
%     metrics.rDhat: rDhat(theta,theta',r')
%     metrics.rDhattotal: rDhat(theta)
%     metrics.rtotal: r(theta)
%     metrics.EXW: Expected progress per stage due to wireless transmissions
%     metrics.EC: Expected cost per stage
%     metrics.EDelta: Expected duration of each stage
%     metrics.EXB: Expected progress per stage due to psysical transport
%     metrics.Vp: the packet speed
%     metrics.Cp: the packet cost
%     metrics.limitingEX: this is a late addition made during the first
%       revision. It gives the expected progress in a transmission in the
%       limiting desnity regime. see the details in the paper. 
%     metrics.limitingEC: this is a late addition made during the first
%       revision. It gives the expected cost in a transmission in the
%       limiting desnity regime. see the details in the paper.
%      
% Comment on use of indices: 
% ==========================
%      * 1<=i<=N describes thetas.
%      * 1<=j<=M describes locations.
%      * 1<=k<=N*M descibes pairs of thetas and locations.  
%      * 1<=m<=S describes locations on boundary curves.
%

function metrics=metrics(node_pars,FR_type,potential_type,cost_type,...
                                           density_type,num_pars,verbose)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART A
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=node_pars(1);
v0=node_pars(2);
r0=node_pars(3);
N=num_pars(1); 
B=num_pars(2); 
L=num_pars(3); 
S=num_pars(4); 
dA=(2*B/L)^2;  % This is the incremental volume of area integrals. 
               % Each point is associated with a square such that the 
               % point is at the center of the square. All squares fall 
               % completely within the box.
dtheta=2*pi/N; % This is the incremental volume of angle integrals.
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART B
% We create the thetas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We create thetas
thetas=-pi+(pi/N)*(2*(1:1:N)-1); % ATTN, we do not want to have either -pi or pi
thetas=(thetas-fliplr(thetas))/2; % we use this trick to ensure that 
            % thetas are *perfectly* symmetric, otherwise results are not
            % symmetric
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART C
% We create the locations (x_j,y_j) and M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We create the grid of points
xx=-B+(B/L)*(2*(1:1:L)-1);
xx=(xx-fliplr(xx))/2; % we use this trick to enforce symmetry
yy=xx;
if verbose>=1, disp(['Grid points created: ' num2str(L^2)]); end

% We keep some points in the grid
xin=zeros(L^2,1);
yin=zeros(L^2,1);
M=0;
for j1=1:L
    for j2=1:L
        if inFR(FR_type,xx(j1),yy(j2))
             M=M+1;
             xin(M)=xx(j1);
             yin(M)=yy(j2);
        end
    end
end
if verbose>=1, disp(['Grid points added: ' num2str(M)]); end

xin=xin(1:M,1);
yin=yin(1:M,1);
r=sqrt(xin.^2+yin.^2);
phi=atan2(yin,xin);
locations=[xin yin r phi];

% We plot those points
if verbose>=2
   newfig(1,[2 2 9 9]); hold on; 
   scatter(r.*cos(phi),r.*sin(phi),30,'black','filled');
   xbound=zeros(1000,1);
   ybound=zeros(1000,1);
   for i=1:1000
       angle=-pi+(i-1)*2*pi/999;
       radius=feval('boundary',FR_type,angle);
       xbound(i)=radius*cos(angle);
       ybound(i)=radius*sin(angle);
   end
   plot(xbound,ybound,'Linewidth',2,'Color','red');
   axis([-B-0.1 B+0.1 -B-0.1 B+0.1]);
   axis equal;
   xlabel('$x$','Fontsize',10);
   ylabel('$y$','Fontsize',10);
   title('The boundary curve and grid points used','Fontsize',10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART D
% Late addition during first revision. We calculate limiting EX and EC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metrics.limitingEX=sum(xin)/M;

metrics.limitingEC=(sum(xin.^2)+sum(yin.^2))/M;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1
% We calculate E(N;theta,r) and P_E(theta,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose>=1, disp('PART 1: Starting calculations for E(N;theta,r) and P_E(theta,r)'); end

% We calculate an integral, I1, we will be using frequently. 
% This holds with perfect accuracy, when density is uniform:
% I1=2*abs(thetas);
% However, we use the following to ensure that results with slow and fast 
% methods match. 
I1=zeros(N,1);
for i=1:N
    int=0;
    for iprime=1:N
        if abs(thetas(iprime))<abs(thetas(i))
            int=int+direction_density(density_type,thetas(iprime));
        end
    end
    I1(i)=int*dtheta;
end
 
% We calculate another integral, I2, we will also be using frequently.
I2=zeros(M,1);
for j=1:M
    x=locations(j,1);
    y=locations(j,2);
    int=0;
    for jprime=1:M 
       xprime=locations(jprime,1);
       yprime=locations(jprime,2);
       if ~inFR(FR_type,xprime+x,yprime+y)
            int=int+1;
       end
    end
    I2(j)=int*dA;
end
%end;

% Now we calculate EN and PE
if potential_type==1
   EN=I1*I2'*lambda;
else
   EN=zeros(N,M);

   pot_check=zeros(N,M);
   for iprime=1:N
       for jprime=1:N
          pot_check(iprime,jprime)=potential(potential_type,thetas(iprime),...
                                    locations(jprime,1),locations(jprime,2));
       end
   end
   
   for i=1:N
       if verbose>=1, disp([num2str(i) '/' num2str(N)]); end
       pot_threshold=potential(potential_type,thetas(i),0,0);
       for j=1:M
          for iprime=1:N
              for jprime=1:M
                if pot_check(iprime,jprime)>pot_threshold
                    if ~inFR(FR_type,locations(jprime,1)+locations(j,1),locations(jprime,2)+locations(j,2))
                        EN(i,j)=EN(i,j)+direction_density(density_type,thetas(iprime));
                    end
                end
              end
          end
       end
   end
   EN=lambda*dA*dtheta*EN;
end

PE=exp(-EN);

if verbose>=1, disp('Ended calculations for E(N;theta,r) and P_E(theta,r)'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2
% We calculate g(theta',r';theta,r),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=zeros(N,M,N,M);

if verbose>=1, disp('PART 2: Starting calculations for g'); end

if potential_type==1
    outerbound=floor(N+1/2);
else
    outerbound=N;
end
for iprime=1:outerbound % ATTN, If the direction density function is not an even function, errors are introduced in this loop
    if verbose>=1, disp([num2str(iprime) '/' num2str(outerbound)]); end
    if potential_type==1
        rng=(1:iprime-1); % not needed: (N-iprime+2:N)];
    else
        rng=(1:1:N);
    end
    for jprime=1:M
    %disp(['  INNER LOOP: ' num2str(jprime) '/' num2str(M)]);
        xprime=locations(jprime,1);
        yprime=locations(jprime,2);    
        pot2=potential(potential_type,thetas(iprime),xprime,yprime);
        for i=rng
            pot1=potential(potential_type,thetas(i),0,0);
            if (pot1<pot2) % we only consider transmitting to nodes with 
                           % strictly better potential. Equality can
                           % only happen due to discretization
               for j=1:M
                  % WE HAVE NOW FIXED theta',r',theta,r
                  x=locations(j,1);
                  y=locations(j,2);                
              
                  if ~inFR(FR_type,xprime+x,yprime+y)
                    
                     % first, we calculate the integral
                     if potential_type==1
                        int=I1(iprime)*I2(j); % we avoid the costly integral calculation 
                     else
   
                        % we will have to calculate the integral
                        int=0;
                        for i2prime=1:N
                           for j2prime=1:M
                              x2prime=locations(j2prime,1);
                              y2prime=locations(j2prime,2);
                              pot3=potential(potential_type,thetas(i2prime),...
                                             x2prime,y2prime);
                              if (pot3>pot2) % here, you make a mistake due 
                                 % to discretization no matter how you pick
                                 % the inequality, because there may be more
                                 % than 1 candidates!
                                 if ~inFR(FR_type,x2prime+x,y2prime+y)
                                     int=int+direction_density(density_type,thetas(i2prime));
                                 end
                              end      
                           end
                        end
                        int=int*dA*dtheta;            
                     end % we calculated the integral
                     
                     g(iprime,jprime,i,j)=lambda*direction_density(density_type,thetas(iprime))*exp(-lambda*int);                          
                     
                     if potential_type==1 % we take advantage of symmetries
                         g(iprime,jprime,N-i+1,j)=g(iprime,jprime,i,j);                        
                         g(N-iprime+1,jprime,i,j)=g(iprime,jprime,i,j);                      
                         g(N-iprime+1,jprime,N-i+1,j)=g(iprime,jprime,i,j);
                     end
                     
                  end % finished FR condition
               end %finished j
            end % finished potential condition
        end % finished i
    end % finished jprime 
end % finished iprime
if verbose>=1, disp('Ended calculations for g'); end
   
% We verify that for fixed theta and r, PE(theta,r) plus the integral of g
% is unity. 

if verbose>=2
    check=zeros(N,M);
    for i=1:N
        for j=1:M
            for iprime=1:N
                for jprime=1:M
                    check(i,j)=check(i,j)+g(iprime,jprime,i,j);
                end % finished jprime
            end % finished iprime
            check(i,j)=check(i,j)*dtheta*dA+PE(i,j);
        end % finished j
    end % finished i
   newfig(2,[4 3 9 9]); hold on;
   surf(ones(N,M)-check);
   title('$1-P_E(\theta,\mathbf{r})-\int_{-\pi}^\pi \int_{\mathcal{F}} g(\theta'',\mathbf{r}'';\theta,\mathbf{r})\, dA''d\theta'' $');
   xlabel('$\mathbf{r}$');
   ylabel('$\theta$');
   view(3);
   disp(['Maximum error of part 2:' num2str(max(max(abs(ones(N,M)-check))))]);
   clear check;
end % finished check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3
% We calculate rA(theta,theta')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose>=1, disp('PART 3: We start the calculation of rA(theta,thetaprime) and rA(theta)'); end
rA=zeros(N,N);
rAtotal=zeros(N,1);

for i=1:N
    if verbose>=1, disp([num2str(i) '/' num2str(N)]); end
    pot1=potential(potential_type,thetas(i),0,0);     
    for iprime=1:N
        pot3=potential(potential_type,thetas(iprime),0,0);

        % We calculate rA for fixed theta and theta'
        if pot3>=pot1
            rA(i,iprime)=r0*direction_density(density_type,thetas(iprime));
        else
            % we will have to calculate an integral
            int=0;                   
            for i2prime=1:N
                for j2prime=1:M
                    pot2=potential(potential_type,thetas(i2prime),locations(j2prime,1),locations(j2prime,2));
                    if (pot1>=pot2) && (pot2>pot3) % If candidate is as good as old direction, but better than new direction, 
                                                    % we allow it. makes better sense in dicretization
                        int=int+direction_density(density_type,thetas(i2prime));
                    end
                end
            end
            rA(i,iprime)=r0*direction_density(density_type,thetas(iprime))... 
               *exp(-dA*dtheta*int*lambda); 
        end
        rAtotal(i)=rAtotal(i)+rA(i,iprime);
    end
end
rAtotal=rAtotal*dtheta;

if verbose==2
    
   % We prepare the figure
      
   newfig(3,[6 4 9 9]); hold on;

   surf1=surf(thetas,thetas,rA);
   set(surf1,'FaceColor', 'texturemap'); % makes sure coloring is symmetric
   ylabel('$\theta$','Fontsize',10);
   xlabel('$\theta''$','Fontsize',10);
   title('Incremental rate $r_A(\theta,\theta'')$','Fontsize',10);
   thetaaxes(1,1);
end

if verbose>=1, disp('The calculation of rA(theta,thetaprime) and rA(theta) is finished'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4
% We calculate rB(theta,theta',r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose>=1, disp('PART 4: We start the calculation of rB(theta,thetaprime,rprime) and rB(theta)'); end

% We calculate two integrals
AreaFR=M*dA;

I3=zeros(N,1);
for i=1:N
    int=0;
    for iprime=1:N
        if abs(thetas(i))<abs(thetas(iprime))
            int=int+direction_density(density_type,thetas(iprime));
        end
    end
    I3(i)=int*dtheta;
end

I4=zeros(N,N);
for i=1:N
    for iprime=1:N
        int=0;
        for i2prime=1:N
            a1=abs(thetas(i));
            a2=abs(thetas(i2prime));
            a3=abs(thetas(iprime));
            if a1<=a2 && a2<a3    %if a1<=a2 && a2<=a3, 
                int=int+direction_density(density_type,thetas(i2prime));
            end
        end
        I4(i,iprime)=int*dtheta;
    end
end

rB=zeros(N,N,M);
rBtotal=zeros(N,1);

for i=1:N
    if verbose>=1, disp([num2str(i) '/' num2str(N)]); end
    pot1=potential(potential_type,thetas(i),0,0);
    for iprime=1:N
        for jprime=1:M
            pot2=potential(potential_type,thetas(iprime),locations(jprime,1),locations(jprime,2));
            if pot1>=pot2
                
                if potential_type==1
                    coeff1=I3(iprime);
                    coeff2=exp(-lambda*AreaFR*I4(i,iprime));
                else
                    
                    % We calculate 1-D integral
                    int=0;
                    for i2prime=1:N 
                        pot3=potential(potential_type,thetas(i2prime),0,0);
                        if pot2>pot3
                            int=int+direction_density(density_type,thetas(i2prime));
                        end
                    end
                    coeff1=int*dtheta;
                    
                    % We calculate 3-D integral
                    int=0;
                    for i3prime=1:N
                        for j3prime=1:N
                           pot4=potential(potential_type,thetas(i3prime),locations(i3prime,1),locations(i3prime,2));
                           if pot1>=pot4 && pot4>pot2
                               int=int+direction_density(density_type,thetas(i3prime));
                           end
                        end
                    end
                    coeff2=exp(-lambda*int*dtheta*dA);
            
                end
                                                
                % We wrap up the calculation for given triplet
                rB(i,iprime,jprime)=r0*lambda*direction_density(density_type,thetas(iprime))*coeff1*coeff2;
                rBtotal(i)=rBtotal(i)+rB(i,iprime,jprime);
            end
        end
    end
end
rBtotal=rBtotal*dtheta*dA;

if verbose>=1, disp('The calculation of rB(theta,thetaprime,rprime) and rB(theta) is finished'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 5
% We calculate rC(theta,theta',r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose>=1, disp('PART 5: We start the calculation of rC(theta,thetaprime,rprime) and rC(theta)'); end
rC=zeros(N,N,M);
rCtotal=zeros(N,1);

for i=1:N
    if verbose>=1, disp([num2str(i) '/' num2str(N)]); end
    pot1=potential(potential_type,thetas(i),0,0);
    for iprime=1:N
        for jprime=1:M 
           xprime=locations(jprime,1);
           yprime=locations(jprime,2);           
           pot2=potential(potential_type,thetas(iprime),xprime,yprime);
           if pot2>pot1 % must be strictly better to do a transmission
               int=0;
               for i2prime=1:N
                  pot3=potential(potential_type,thetas(i2prime),xprime,yprime);
                  if pot3<pot1
                      int=int+direction_density(density_type,thetas(i2prime));
                  end
               end
               coef1=lambda*r0*direction_density(density_type,thetas(iprime));
               rC(i,iprime,jprime)=coef1*int*dtheta;
           end
           rCtotal(i)=rCtotal(i)+rC(i,iprime,jprime);
        end
    end
end
rCtotal=rCtotal*dA*dtheta;

if verbose>=2
   newfig(4,[8 5 9 9]); hold on;
   plot_2d_fun(rC(round(N/4),round(N/4)+1,:),locations,[-0.5 2 -1.5 1.5]);
   title(['$r_C(' num2str(thetas(round(N/4))) ',' num2str(thetas(round(N/4)+1)) ',\mathbf{r})$']);
end

if verbose>=1, disp('The calculation of rC(theta,thetaprime,rprime) and rC(theta) is finished'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 6
% We calculate rD(theta,theta',s) and rD(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose>=1, disp('PART 6: We start the calculation of rD(theta,theta'',s) and rD(theta)'); end

% We create the parametrization
s=(1/(2*S))+(0:1:S-1)/S; % we avoid both s=0 and s=1, [0,1] interval breaks in S equals pieces. 

rD=zeros(N,N,S);
rDtotal=zeros(N,1);

% POTENTIAL 1 CASE
if potential_type==1
    [~,~,tx,ty,der]=threshold_curve(potential_type,FR_type,-pi/2,0,s);

    for i=1:N
       if verbose>=1, disp([num2str(i) '/' num2str(N)]); end
       for iprime=better_angles(N,i)
          for m=1:S
             j=sqrt(-1);
             veccom=exp(j*thetas(i))-exp(j*thetas(iprime));
             vx=real(veccom);
             vy=imag(veccom);
             inprod=tx(m)*vx+ty(m)*vy;     
             if inprod>0
                coef1=lambda*v0*direction_density(density_type,thetas(iprime));
                rD(i,iprime,m)=coef1*inprod*der(m);
                rDtotal(i)=rDtotal(i)+rD(i,iprime,m)*dtheta/S; 
             end
          end
       end
    end

% OTHER CASES
else 

   for i=1:N
      if verbose>=1, disp([num2str(i) '/' num2str(N)]); end
      for iprime=1:N
         [~,~,tx,ty,der]=threshold_curve(potential_type,FR_type,thetas(i),thetas(iprime),s);
         for m=1:S
            j=sqrt(-1);
            veccom=exp(j*thetas(i))-exp(j*thetas(iprime));
            vx=real(veccom);
            vy=imag(veccom);
            inprod=tx(m)*vx+ty(m)*vy;     
            if inprod>0
               coef1=lambda*v0*direction_density(density_type,thetas(iprime));
               rD(i,iprime,m)=coef1*inprod*der(m);
               rDtotal(i)=rDtotal(i)+rD(i,iprime,m)*dtheta/S;
            end
         end
      end
   end

end

if verbose>=1, disp('The calculation of rD(theta,thetaprime,s) and rD(theta) is finished'); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 7
% We calculate rDhat(theta,theta',r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose>=1, disp('PART 7: We start the calculation of rDhat(theta,theta'',r) and rD(theta)'); end
rDhat=zeros(N,N,M);
rDhattotal=zeros(N,1);

if potential_type==1
   [x,y,~,~,~]=threshold_curve(potential_type,FR_type,-pi/2,0,s);
   j=knnsearch(locations(:,1:2),[x' y']);
end

for i=1:N
   for iprime=better_angles(N,i)
        
       if potential_type~=1
           [x,y,~,~,~]=threshold_curve(potential_type,FR_type,thetas(i),thetas(iprime),s);
           j=knnsearch(locations(:,1:2),[x' y']);
       end
       
       for m=1:S
          rDhat(i,iprime,j(m))=rDhat(i,iprime,j(m))+rD(i,iprime,m)/(S*dA);
          rDhattotal(i)=rDhattotal(i)+rD(i,iprime,m)*dtheta/S;
       end
   end
end

if verbose>=2
   newfig(5,[10 6 9 9]); hold on;
   plot_2d_fun(rDhat(round(N/4),round(3*N/4)-2,:),locations,[-0.5 2 -1.5 1.5]);
   title(['$\hat{r}_D(' num2str(thetas(round(N/4))) ',' num2str(thetas(round(3*N/4)-2)) ',\mathbf{r})$']);
end

if verbose>=1 
   disp('The calculation of rDhat(theta,thetaprime,r) is finished'); 
   disp(['The integral of rD is ' num2str(dtheta*dtheta*(1/S)*sum(sum(sum(rD))))]);
   disp(['The integral of rDhat is ' num2str(dtheta*dtheta*dA*sum(sum(sum(rDhat))))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 8
% We calculate rtotal(theta) and draw plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rtotal=rAtotal+rBtotal+rCtotal+rDhattotal;

if verbose>=2
   newfig(6,[12 7 9 9]); hold on; 
   plot(thetas,rAtotal,'Linestyle','-');
   plot(thetas,rBtotal,'Linestyle','-.');
   plot(thetas,rCtotal,'Linestyle',':');
   plot(thetas,rDtotal,'Linestyle','--');
   plot(thetas,rAtotal+rBtotal,'-.r','Linewidth',2);
   plot(thetas,rtotal,':b','Linewidth',2);
   plot(thetas,rDhattotal,':r','Linewidth',2);
   title('Aggregate rates','Fontsize',10);
   hhh=legend('$r_\mathcal{A}(\theta)$','$r_\mathcal{B}(\theta)$','$r_\mathcal{C}(\theta)$',...
      '$r_\mathcal{D}(\theta)$','$r_\mathcal{A}(\theta)+r_\mathcal{B}(\theta)$',...
      '$r_\mathcal{A}(\theta)+r_\mathcal{B}(\theta)+r_\mathcal{C}(\theta)+r_\mathcal{D}(\theta)$',...
      '$\hat{r}_\mathcal{D}(\theta)$');
   set(hhh,'Interpreter','latex');
   thetaaxes(1,0);
   xlabel('$\theta$','Fontsize',10);
end

if verbose>=1, disp('The calculation of rDhat(theta,thetaprime,r) is finished'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 8
% We calculate the matrix of the eigenvalue problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose>=1, disp('PART 8: We start the construction of the eigenvalue matrix'); end

Kernel=zeros((M+1)*N,(M+1)*N);

for i1=1:N %psi_B
    for i2=1:N
        Kernel(i1,i2)=dtheta*rA(i2,i1)/rtotal(i2);
    end 
    
    for jprime=1:M
        i2=(jprime-1)*N+i1;
        Kernel(i1,N+i2)=PE(i1,jprime);
    end
end
 
for i1=1:M*N  % psi_W
    
    j=floor((i1-1)/N)+1;
    i=i1-(j-1)*N;
    
    for i2=1:N 
        Kernel(i1+N,i2)=dA*dtheta*(rB(i2,i,j)+rC(i2,i,j)+rDhat(i2,i,j))/rtotal(i2);
    end

    for i2=1:M*N 
        jprime=floor((i2-1)/N)+1;    
        iprime=i2-(jprime-1)*N;
        Kernel(i1+N,i2+N)=dtheta*dA*g(i,j,iprime,jprime);
    end
    
end

clear g; 

if verbose>=2
   newfig(7,[14 8 9 9]); hold on; 
   plot(sum(Kernel));   
   title('The column sums of the stochastic matrix');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 9
% We calculate PsiB and PsiW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose>=1, disp('PART 9: We start the calculation of the steady state distribution'); end

% We find the distribution
[V,~]=eigs(Kernel,1);
distribution=V(:,1)/(sum(V(:,1)));
 
psiB=distribution(1:N,1)/dtheta;

psiW=zeros(N,M);
for index=1:N*M
    j=floor((index-1)/N)+1;
    i=index-(j-1)*N;
    psiW(i,j)=distribution(N+index);
end
psiW=psiW/(dtheta*dA);

if verbose>=2
   newfig(8,[16 9 9 9]); hold on; 
   plot(thetas,psiB);   
   title('$\psi_B$');
   thetaaxes(1,0);
   
   newfig(9,[18 10 9 9]); hold on; 
   plot_2d_fun(psiW(floor(N/4),:),locations,[-0.5 2 -1.5 1.5]);
   title(['$\psi_W(' num2str(thetas(round(N/4))) ',\mathbf{r})$']);
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 10
% We calculate all expectations and metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose>=1, disp('PART 10: We start the calculation of expectations and metrics'); end

EXW=0;
EC=0;
EDelta=0;
EXB=0;

for i=1:N
    for j=1:M
        EXW=EXW+psiW(i,j)*locations(j,1)*dtheta*dA;
        EC=EC+psiW(i,j)*cost(cost_type,locations(j,1),locations(j,2))*dtheta*dA;
    end
end

for i=1:N
    EDelta=EDelta+psiB(i)*dtheta/rtotal(i);
    EXB=EXB+psiB(i)*v0*cos(thetas(i))*dtheta/rtotal(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART D
% We create the output structures. Used for keeping the number of output
% variables small
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metrics.thetas=thetas;
metrics.locations=locations;
metrics.EN=EN; 
metrics.PE=PE;
metrics.rA=rA;  
metrics.rAtotal=rAtotal; 
metrics.rB=rB;  
metrics.rBtotal=rBtotal; 
metrics.rC=rC;  
metrics.rCtotal=rCtotal;
metrics.rD=rD;  
metrics.rDtotal=rDtotal; 
metrics.rDhat=rDhat;  
metrics.rDhattotal=rDhattotal; 
metrics.rtotal=rtotal; 
metrics.EXW=EXW;
metrics.EC=EC;
metrics.EDelta=EDelta;
metrics.EXB=EXB;
metrics.Vp=(EXB+EXW)/EDelta;
metrics.Cp=EC/(EXB+EXW);



