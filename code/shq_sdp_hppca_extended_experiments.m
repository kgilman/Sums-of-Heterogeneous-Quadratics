close all
clear
clc


dims = 10:10:50;
% dims = 50;
rank = [3,5,7,10];
num_trials = 100;
% num_trials = 1;

%% HPPCA
results_hppca = zeros(length(dims),length(rank),num_trials);

ns = [5, 10, 20, 50, 100, 200];
% n = [50,200];
vs = [1;2;3;4];

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
                    if(t==1 || mod(t,10)==0)
                        fprintf('\n dim: %i, rank: %i, n = %i, v = %i, trial: %i',d,k,ns(m),vs(p),t);
                    end
                    U = orth(randn(d,k));
                    [M,~] = hppca_problem(U,lambda,n,v);
                    [proj_err,Xi_err,cvx_optval,Uhat,X,nu,Z,Y] = solve_sdp_CVX(M);
                    results_hppca(i,j,m,p,t) = Xi_err;
                end

            end
        
        end
    end
end

save('results_hppca_extended.mat', 'results_hppca');
