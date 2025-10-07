close all
clear
clc

thresh = 1e-5;


%% HPPCA
% load results_hppca.mat results_hppca

% load results_hppca_extended_large_d300_k10_n250.mat
% load results_hppca_extended_large_d300_k05_n50.mat
load results_hppca_extended_large_d300_k10_n50.mat
% test = max(results_hppca,[],3)

num_trials = 20;

% array = 1 - sum(results_hppca > thresh,3) ./ num_trials;

% results_hppca_squeeze = squeeze(results_hppca);
results_hppca_squeeze = results_hppca(1,:,1,1,:);

array = 1 - sum(results_hppca_squeeze > thresh,5) ./ num_trials;
% dims = {'d=10';'d=20';'d=30';'d=40';'d=50'};
dims = {'d=300'};
r5 = array(:,1);
r10 = array(:,2);


T = table(r5,r10,'RowNames',dims)
writetable(T,'results_hppca_larger_d.dat','WriteRowNames',true) 