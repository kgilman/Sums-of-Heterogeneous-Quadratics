close all
clear all
clc



% dims = 300;
% dims = 10;
dims = 20;
rank = [5,10];
num_trials = 20;

%% HPPCA

% ns = [5, 10, 20, 50, 100, 200];
% ns = [20, 50];
% ns = [50];
ns = [10];
vs = [3];


results_hppca = zeros(length(dims),length(rank),length(ns), length(vs), num_trials);
time_hppca = zeros(length(dims),length(rank),length(ns), length(vs), num_trials);
toctime = 0;

fprintf('\n HPPCA')
for i=1:length(dims)
    d = dims(i);
    for j=1:length(rank)
        k = rank(j);
        lambda = linspace(1,4,k);
        
        for m=1:length(ns)
        
            n = [ns(m); 4*ns(m)];

            for p = 1:length(vs)
                v = [1; vs(p)];
                for t = 1:num_trials
                    
                    tic
                    
                    if(t==1 || mod(t,10)==0)
                        fprintf('\n Elapsed time %s',toctime);
                        fprintf('\n dim: %i, rank: %i, n = %i, v = %i, trial: %i \n',d,k,sum(n),vs(p),t);
                    end
                    U = orth(randn(d,k));
                    [M,~] = hppca_problem(U,lambda,n,v);
                    [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
                    
                    toctime = toc;
                    results_hppca(i,j,m,p,t) = Xi_err;
                    time_hppca(i,j,m,p,t) = toctime;
                    
                    save(sprintf("results_hppca_extended_large_d%02d_k%02d_n%02d.mat",d,k,sum(n)), 'results_hppca', 'dims', 'rank', 'ns', 'vs', 'num_trials', 'time_hppca')
                end

            end
        
        end
    end
end

% save('results_hppca_extended_large_d.mat', 'results_hppca', 'dims', 'rank', 'ns', 'vs', 'num_trials', 'time_hppca');
