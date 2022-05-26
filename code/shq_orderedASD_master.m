clear
close all
clc
set(0,'defaulttextinterpreter','latex')
params.tol = 0.2;
params.contraction = 0.5;
params.stepsize = 1.0;
params.maxsearches = 50;
params.break_tol = 1e-10;
params.niters = 2000;


rng(0);
d = 10;
k = 3;
% sig = 0.1;
level = 1;
rank = k;
% rank = d;
% step = sig;

sigs = [1e-5,1e-4,1e-3,1e-2,1e-1,0.25,0.5,0.75,1];
% sigs = [0,1e-5];
% sigs = [10];
% sigs = [1];
num_trials = 100;

results_log = zeros(length(sigs),num_trials,6);
stats_log = zeros(length(sigs),num_trials,2);
for i=1:length(sigs)
    sig = sigs(i);
    step = 1;
    
    for j=1:num_trials
        
        if(mod(j,10)==0 || j == 1)
            fprintf("Sigma: %f,  Trial: %i \n",sig,j)
        end
        
        %%%% Generate the Mi
        M = generateOrderedASD(d,rank,k,level,sig,step);
        M = normalizeM(M);

        %%%% Measure the Mi stats
        rho = zeros(d,1);
        for m=1:k
            Mi = M{m};
            rho(m) = max((2*diag(abs(Mi)) - sum(abs(Mi),2) )./ sum(abs(Mi),2));
        end
        stats_log(i,j,1) = mean(rho);
        
        
        c = combnk(1:k,2);
        commute_norm = zeros(length(c),1);
        for m=1:length(c)
            commute_norm(m) = norm(M{c(m,1)}*M{c(m,2)}  - M{c(m,2)}*M{c(m,1)},'fro');
        end
        
        stats_log(i,j,2) = max(commute_norm);
        
        %%%% Run experiment
        [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
        
        U0 = orth(randn(d,k));
        fxn = @(U) F(U,M);
%         [Ustga,optval_stga] = runStGA(M,U0,params,fxn);
        [Umm,optval_mm] = runStMM(M,U0,params,fxn);

%         gap = fxn(Uhat) - optval_mm(end);
        gap = abs(cvx_optval - optval_mm(end));

        [result,nu_hat] = feasibilityCheck_CVX(Umm,M); 
        
        Udist = norm(abs(Umm'*Uhat) - eye(k),'fro') / norm(eye(k),'fro');
        
        sonc = checkSONC(Umm,M);

        results_log(i,j,:) = [proj_err, Xi_err, Udist, gap, result, sonc];
    end

end

% figure,
% semilogx(sigs,mean(stats_log(:,:,1),2));
% figure,
% semilogx(sigs,mean(stats_log(:,:,2),2));

%%

marker_size = 75;
fontsize = 20;

fig1 = figure('Visible', 'off');
ax1 = axes('Parent', fig1);
hold(ax1, 'on');
fig2 = figure('Visible', 'off');
ax2 = axes('Parent', fig2);
hold(ax2, 'on');
fig3 = figure('Visible', 'off');
ax3 = axes('Parent', fig3);
hold(ax3, 'on');

tight_thresh = 1e-5;

for i=1:length(sigs)
   
    for j=1:num_trials
        
        sig = stats_log(i,j,2);
        result = results_log(i,j,5);
%         proj_err = results_log(i,j,1); %%change this
        proj_err = results_log(i,j,2); %%change this
        tight = proj_err < tight_thresh;
        
%         color = pointColor(result,tight);
        if(~tight)
            marker_type = 'x';
            color = 'red';
        elseif(result)
            marker_type = '^';
            color = 'blue';
        else
            marker_type = 'o';
            color = 'black';
        end

        scatter(ax1,sig,proj_err,marker_size,marker_type,'MarkerEdgeColor',color);
        scatter(ax2,sig,results_log(i,j,3),marker_size,marker_type,'MarkerEdgeColor',color);
        scatter(ax3,sig,results_log(i,j,4),marker_size,marker_type,'MarkerEdgeColor',color);

    end
end

hold(ax1, 'off');
hold(ax2, 'off');
hold(ax3, 'off');

grid(ax1, 'on');
grid(ax2, 'on');
grid(ax3, 'on');

set([fig1, fig2, fig3], 'Visible', 'on');

set(ax1, 'YScale', 'log','XScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')
set(ax2, 'YScale', 'log','XScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')
set(ax3, 'YScale', 'log','XScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')

figure(1),xlabel('Max commuting distance','Interpreter','latex'),title('SDP Projection error')
figure(2),xlabel('Max commuting distance','Interpreter','latex'),title('$\bf{U}$ distance')
figure(3),xlabel('Max commuting distance','Interpreter','latex'),title('Objective value gap')
%%

% figure
% t=linspace(0,10,100);
% plot(t,sin(t));
% hold on;

marker_size = 10;
figure(1)
set(gca,'FontSize',25)
hold on
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'xr','MarkerSize',marker_size);
h(2) = plot(NaN,NaN,'^b','MarkerSize',marker_size);
h(3) = plot(NaN,NaN,'ok','MarkerSize',marker_size);
leg1 = legend(h, 'Stationary: not tight','Stationary: tight','Global Opt.','Location','northwest');
set(leg1,'Interpreter','latex')
% figure()
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'xr');
% h(2) = plot(NaN,NaN,'^b');
% h(3) = plot(NaN,NaN,'ok');
% legend(h, 'red','blue','black');

%%
tight_counts = zeros(length(sigs),2);
% num_stat_pts = zeros(length(sigs),1);
for i=1:length(sigs)
%     tight_counts(i,1) = 0;
    tight_counts(i,2) = sum(results_log(i,:,1) > tight_thresh);
    
    for j=1:num_trials 
        result = results_log(i,j,5);
        proj_err = results_log(i,j,1);
        tight = proj_err < tight_thresh;
        
        if(tight && result)
            tight_counts(i,1) =  tight_counts(i,1) + 1;
        end

    end
end


figure(5),
bar(tight_counts ./ num_trials * 100);
% ylim([0 1])
ax = gca;
% ax.XTick = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '0.25' '0.5' '0.75' '1'};
set(gca,'xticklabel',{'$10^{-5}$' '$10^{-4}$' '$10^{-3}$' '$10^{-2}$' '$10^{-1}$' '0.25' '0.5' '0.75' '1'},'TickLabelInterpreter','latex')
set(ax,'FontSize',fontsize)
xlabel('$\sigma$','Interpreter','latex')
ylabel('Fraction of trials','Interpreter','latex')
title('Global certification')
leg1 = legend('Stationary: tight', 'Stationary: not tight','Location','northwest');
set(leg1,'Interpreter','latex')
% set(gcf,'Position',[100 100 450 350])


% figure(6),
% semilogx(sigs,mean(stats_log(:,:,1),2));
% figure(7),
% semilogx(sigs,mean(stats_log(:,:,2),2));




function M = generateOrderedASD(d,rank,k,level,sig,step)
M = {};
n = 10*d;
A = sqrt(sig)*randn(d,n);
N = 1/n*A*A';
% D = level*diag(sort(rand(d,1),'descend'));
D = level*diag([sort(rand(rank,1),'descend'); zeros(d-rank,1)]);
M{k} = D + N;
for i=k-1:-1:1
    A = sqrt(sig)*randn(d,n);
    N = 1/n*A*A';
    D = level*diag([sort(rand(rank,1),'descend'); zeros(d-rank,1)]);
    M{i} = M{i+1} + D + step*N;
%     ispsd(M{i});
end

checkOrdered(M)

end


function M = normalizeM(M)
    k = length(M);
    norms = zeros(k,1);
    for i=1:k
        norms(i) = norm(M{i},2);
    end
    max_norm = max(norms);
    for i=1:k
        M{i} = M{i} ./ max_norm;
    end
end


function M = generateASD(d,k,level,sig)
M = {};
for i=1:k
    A = randn(d);
    N = sig*A*A';
    D = level*diag(sort(rand(d,1),'descend'));
    M{i} = sym(D + N);
end

end


function checkOrdered(M)
k = length(M);
    for i=1:k-1
        if(~ispsd(M{i} - M{i+1},1e-7))
            error("M's are not ordered PSD")
        end
    end
end

function result = sym(A)
result = 0.5*(A + A');
end

function result = checkSONC(U,M)
%%%% return 1 if SONC is satisfied; 0 otherwise

     [~,k] = size(U);
     c = combnk(1:k,2);
     temp = zeros(length(c),1);
     for k=1:length(c)
         i = c(k,1); j = c(k,2);
         temp(k) = (U(:,i)'*M{i}*U(:,i) + U(:,j)'*M{j}*U(:,j) > U(:,j)'*M{i}*U(:,j) + U(:,i)'*M{j}*U(:,i));
     end
     result = all(temp > 0);
end
% function colorString = pointColor(result)
%     if(result)
%         colorString = 'red';
%     else
%         colorString = 'black';
%     end
% end
