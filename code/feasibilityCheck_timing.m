clear
close all
clc

params.tol = 0.2;
params.contraction = 0.5;
params.stepsize = 1.0;
params.maxsearches = 50;
params.break_tol = 1e-9;
params.niters = 2000;

rng(0);
% d = 10;
% k = 3;
% lambda = [4;2;1];
v = [1;4];

% n = [200,800];
n = [100,400];

% ds = [10,50,100,200,300,400,500,600,700,800,900,1000];
% ds = [10,50,100];

ds = 10:10:200;
% ds = 10:10:30;
ranks = [3,5,7,10];

num_trials = 10;
time_stats_log = zeros(length(ds),2,2);

for i=1:length(ds)
    for j=1:length(ranks)
        sprintf('Dimension %i, rank %i',ds(i),ranks(j))
        time_log = zeros(num_trials,2);
        results = zeros(num_trials,1);

        d = ds(i);
        k = ranks(j);
        U = orth(randn(d,k));
        lambda = linspace(1,4,k);
        [M,~] = hppca_problem(U,lambda,n,v);

        for t=1:num_trials

        %%%% Run experiment
        tStart = tic;
            [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
        time_log(t,1) = toc(tStart);

        U0 = orth(randn(d,k));
        fxn = @(U) F(U,M);

    %     gap = cvx_optval - optval_stga;

        tStart = tic;
            [Umm,optval_mm] = runStMM(M,U0,params,fxn);
            [result,nu_hat] = feasibilityCheck_CVX(Umm,M);
        time_log(t,2) = toc(tStart);
        results(t) = result;
        end

        time_stats_log(i,j,1,1) = median(time_log(:,1));
        time_stats_log(i,j,1,2) = std(time_log(:,1));

        time_stats_log(i,j,2,1) = median(time_log(:,2));
        time_stats_log(i,j,2,2) = std(time_log(:,2));
%     results_log{i} = results;
    end
end


%%
% Load the data and plot the results

load timing_results_sdp_v_certificate_extended.mat 

ds = 10:10:200;
% ds = 10:10:30;
ranks = [3,5,7,10];

set(0,'defaulttextinterpreter','latex')
figure,
hold on
grid on


% sdp_styles = ['-','--',':','-.'];
% certificate_styles = ['-d','--d',':d','-.d'];

colors = {
    [0 0.4470 0.7410],
    [0.8500 0.3250 0.0980],
    [0.4940 0.1840 0.5560],
    [0.4660 0.6740 0.1880]
}

for j=1:length(ranks)
    sdp_name = sprintf('SDP: rank %i',ranks(j));
    errorbar(ds(1:end-2),time_stats_log(1:end-2,j,1,1),time_stats_log(1:end-2,j,1,2),'-s','MarkerSize',10,'Color',colors{j},'DisplayName',sdp_name);
end

for j=1:length(ranks)
    cert_name = sprintf('StMM + Cert: rank %i',ranks(j));
    errorbar(ds(1:end-2),time_stats_log(1:end-2,j,2,1),time_stats_log(1:end-2,j,2,2),'-^','MarkerSize',10,'Color',colors{j},'DisplayName',cert_name);
end

% legend("SDP","StMM + Global Certificate",'Location','northwest')
leg1 = legend('show','location','eastoutside');
set(leg1,'interpreter','latex')
set(gca,'FontSize',20)
set(gca, 'YScale', 'log')
xlabel('Data dimension')
set(gca,'TickLabelInterpreter','latex')
% xlim([ds(1) ds(17)])
title('Computation time (s)')
set(gcf,'Position',[100 100 650 400])
xlim([10,180])



