%~ clear all

% ----------------------------------------------------------------------
% run the analysis
% ----------------------------------------------------------------------
% function [f]=limprob0(R0, V0, L, R, N)
%~ R0 = 0.1;
%~ V0 = 1;
%~ L = 0.1;
%~ R = 1;
%~ N = 101;
[f q t] = limprob1(R0, V0, L, R, N, Pmat);
%~ [f q t] = limprob0(R0, V0, L, R, N);

% ----------------------------------------------------------------------
% read the stages file
% ----------------------------------------------------------------------
fp = fopen('log_stages.txt');
% use textscat to read the columns of log_stages
a = textscan(fp,'%f %f %f %f %f %f','Delimiter','\t');
fclose(fp);

% skip fist (column name)
dur       = a{2}(2:end);
dir       = a{3}(2:end);
txDist    = a{4}(2:end);
txCost    = a{5}(2:end);
numHops   = a{6}(2:end);

% remove NaN
nanidx = find(isnan(txDist));
fprintf('%u NaN found in txDist\n', numel(nanidx));
txDist(nanidx)=[];
numHops(nanidx)=[];
numHops(numHops==0)=[];
dir(nanidx)=[];
dur(nanidx)=[];
txCost(nanidx)=[];

% ----------------------------------------------------------------------
% Number of hops
% ----------------------------------------------------------------------
printf("average number of hops = %f\n", mean(numHops));
printf("max number of hops = %f\n", max(numHops));
%~ figure, hist(numHops, [1:max(numHops)]);
%~ waitforbuttonpress();


% ----------------------------------------------------------------------
% stage duration
% ----------------------------------------------------------------------
%~ [dur_hist, dur_x] = hist(dur,100);
%~ figure('Visible','off');
%~ % we expect an exponential distribution
%~ lambda_dur = 1/mean(dur);
%~ dur_pdf = lambda_dur * exp(-lambda_dur.*dur_x); % exp pdf
%~ stem(dur_x, dur_pdf, 'd');
%~ hold all
%~ stem(dur_x, dur_hist/trapz(dur_x,dur_hist));

%~ legend('model','sim')
%~ xlabel('$\Delta$')
%~ ylabel('$f_\Delta(\delta)$')
%~ print -dpdflatex "-S400,300" "delta.tex"

%~ % calculate the analitic mean speed duration
%~ dur_mod = 1/sum(f.*q);
%~ printf("average sim. stage duration = %f\n", 1/lambda_dur);
%~ printf("average mod. stage duration = %f\n", dur_mod);


% ----------------------------------------------------------------------
% average speed due to wireless tx
% ----------------------------------------------------------------------

% calculate Vw from simulation
vw_sim = sum(txDist)/sum(dur);

% calculate Vw from the model
%~ vw_mod = (V0*L*R^2/2) * sum(f .* (sign(t).*(t.*cos(t)-sin(t))) );

printf("average sim V_w = %f\n", vw_sim); 
%~ printf("average model V_w = %f\n", vw_mod); 

% ----------------------------------------------------------------------
% average cost
% ----------------------------------------------------------------------
c_sim = sum(txCost)/sum(dur);
%~ c_mod = 4*V0*L*R^3/pi * sum(f.*(1-cos(t)));
printf("average sim C = %f\n", c_sim);
%~ printf("average model C = %f\n", c_mod);

% ----------------------------------------------------------------------
% average speed due to transport
% ----------------------------------------------------------------------
% read simulation results
fp = fopen('log.txt');
a = textscan(fp,'%f %f %f');
fclose(fp);

theta = a{2}-(a{2}-a{1})/2;
fsim = a{3}/sum(a{3});

% plot
figure('Visible','off');
stem(t,f)
hold all
stem(t,fsim, 'd')

xlim([-pi pi]);
legend('model','sim')
xlabel('$\theta_i$')
ylabel('$\pi_i$')
print -dpdflatex "-S400,300" "limDist.tex"

vp_sim = V0 * sum(cos(dir).*dur)/sum(dur);
vp_mod = V0 * sum(f.*cos(theta));

printf("average sim V_p = %f\n", vp_sim);
printf("average model V_p = %f\n", vp_mod);


% ----------------------------------------------------------------------
% average speed
% ----------------------------------------------------------------------
printf("average sim V = %f\n", vp_sim+vw_sim);
%~ printf("average model V = %f\n", vp_mod+vw_mod);
% ----------------------------------------------------------------------
% save results on text file
% ----------------------------------------------------------------------
%~ fid = fopen("results.txt","w+");
%~ fprintf(fid, "sim E[Delta]\tmodel E[Delta]\tsim V_w\tmod V_w\tsim V_t\tmodel V_t\tsim C_p\tmodel C_p\n%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
             %~ 1/lambda_dur, dur_mod, vw_sim, vw_mod, vp_sim, vp_mod, c_sim, c_mod);
