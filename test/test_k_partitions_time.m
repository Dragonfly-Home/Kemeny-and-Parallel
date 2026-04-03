% test_k_partitions_time.m
% Specifically designed to test the impact of different number of partitions (k_partitions) on the runtime of the parallel Kemeny algorithm
clear; clc;

% Keep custom paths if any
addpath('../src');

% ==================================================
% 1. Environment warm-up
% ==================================================
fprintf('Initializing Parallel Environment (Warm-up)...\n');
if isempty(gcp('nocreate'))
    parpool(); 
end
fprintf('Parallel Environment Ready.\n\n');

mat_name = 'cage10.mat'; % You can replace this with the matrix you want to test here

fprintf('Loading and preprocessing matrix: %s\n', mat_name);

load(['../matrix/', mat_name]);

if exist('Problem', 'var')
    A = spones(Problem.A);
    if ~issymmetric(A)
        A = max(A, A'); 
    end
    G = graph(A);
    [bin, binsize] = conncomp(G);
    idx = binsize(bin) == max(binsize); 
    SG = subgraph(G, idx);
    A = SG.adjacency();
end

n = size(A, 1);

% Use dissect preprocessing to obtain optimal bandwidth reordering
p = dissect(A);
A = A(p, p);

% Convert to transition probability matrix
P = spdiags(A*ones(n,1), 0, n, n) \ A;

fprintf('Computing stationary distribution (pi_vec)...\n');
[pi_vec, ~] = eigs(P', 1, 'largestabs', 'MaxIterations', 10000);
pi_vec = pi_vec / sum(pi_vec);

% ==================================================
% 3. Test the impact of different k_partitions on runtime
% ==================================================
% Define the list of k partitions to test (suggested to include multiples of CPU cores)
k_list = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 24, 32, 48, 64, 80]; 

options.solver = 'lu';

fprintf('\nPerformance test of: %s\n', mat_name);
fprintf('%-15s | %-15s\n', 'k_partitions', 'Time (seconds)');
fprintf(repmat('-', 1, 35));
fprintf('\n');

% Loop through different number of partitions
for i = 1:length(k_list)
    k = k_list(i);
    
    % Call the parallel method
    [~, comp_time] = kemeny_k_parallel(P, k, pi_vec, options);
    
    % Output results to the console
    fprintf('%-15d | %-15.4f\n', k, comp_time);
end

fprintf(repmat('-', 1, 35));
fprintf('\nTest completed.\n');