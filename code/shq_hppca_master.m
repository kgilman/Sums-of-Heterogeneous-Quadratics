clear
close all
clc
set(0,'defaulttextinterpreter','latex')
params.tol = 0.2;
params.contraction = 0.5;
params.stepsize = 1.0;
params.maxsearches = 50;
params.break_tol = 1e-10;
params.niters = 10000;

rng(0);
% d = 10;
d = 50;

% k = 3;
k = 5;
% lambda = [4;2;1];
lambda = linspace(1,4,k);

% v = [0.5;1;2;4];

v = [1;4];

L = length(v);

% ns = [2,5,10,20,50,100,200,500,1000];
% ns = [2,3,4,5,10];
% ns = [1:10,25,50,125,250,500];
ns = [3, 5, 10, 25, 50, 100, 200, 500];

num_trials = 100;
results_log = zeros(length(ns),num_trials,5);
stats_log = zeros(length(ns),num_trials,2);
for i=1:length(ns)
    
%     n = ns(i)*ones(1,L);
    n = [ns(i), 4*ns(i)];
    
    for j=1:num_trials
        
        if(mod(j,10)==0 || j == 1)
            fprintf("n: %i,  Trial: %i \n",ns(i),j)
        end
        
        U = orth(randn(d,k));
        [M,~] = hppca_problem(U,lambda,n,v);
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

        gap = cvx_optval - optval_mm(end);

        [result,nu_hat] = feasibilityCheck_CVX(Umm,M);

        Udist = norm(abs(Umm'*Uhat) - eye(k),'fro') / k;

        results_log(i,j,:) = [proj_err, Xi_err, Udist, gap, result];
    end
end

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

for i=1:length(ns)
    for j=1:num_trials
        
        sig = stats_log(i,j,2);
        result = results_log(i,j,5);
        proj_err = results_log(i,j,1);
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
        scatter(ax3,sig,abs(results_log(i,j,4)),marker_size,marker_type,'MarkerEdgeColor',color);

    end
end

hold(ax1, 'off');
hold(ax2, 'off');
hold(ax3, 'off');

grid(ax1, 'on');
grid(ax2, 'on');
grid(ax3, 'on');

set([fig1, fig2, fig3], 'Visible', 'on');
% set([fig2, fig3], 'Visible', 'on');


set(ax1, 'YScale', 'log','XScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')
set(ax2, 'YScale', 'log','XScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')
set(ax3, 'YScale', 'log','XScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')

% set(ax1, 'YScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')
% set(ax2, 'YScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')
% set(ax3, 'YScale', 'log','FontSize',fontsize,'TickLabelInterpreter','latex')

figure(1),xlabel('Max commuting distance'),title('SDP Projection error')
figure(2),xlabel('Max commuting distance'),title('U distance')
figure(3),xlabel('Max commuting distance'),title('Objective value gap')

figure(2)
set(gca,'FontSize',25)

marker_size = 10;
figure(3)
set(gca,'FontSize',25)
hold on
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'xr','MarkerSize',marker_size);
h(2) = plot(NaN,NaN,'^b','MarkerSize',marker_size);
h(3) = plot(NaN,NaN,'ok','MarkerSize',marker_size);
leg1 = legend(h, 'Stationary: not tight','Stationary: tight','Global Opt.','Location','southwest');
set(leg1,'Interpreter','latex')

%%

tight_counts = zeros(length(ns),2);
% num_stat_pts = zeros(length(sigs),1);
for i=1:length(ns)
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
% ax.XTick = {'1e-5' '1e-4' '1e-3' '1e-2' '1e-1' '0.25' '0.5' '0.75' '1'};
% set(gca,'xticklabel',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '25' '50' '125' '250' '500' },'TickLabelInterpreter','latex')
set(gca,'xticklabel',{ '3' '5' '10' '25' '50' '100' '200' '500'},'TickLabelInterpreter','latex')
set(ax,'FontSize',25)
xlabel('$n = [n_1, 4n_1]$')
ylabel('Fraction of trials')
% title('$d = 10, \lambda = [4,2,1]; v = [0.5,1,2,4]$')
title("Global certification")
leg1 = legend('Stationary points: tight', 'Stationary points: not tight','Location','northeast');
set(leg1,'Interpreter','latex')

% figure(6),
% semilogx(ns,mean(stats_log(:,:,1),2));
% title('Diagonal dominance')
% figure(7),
% semilogx(ns,mean(stats_log(:,:,2),2));
% title('Max commuting distance')


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



