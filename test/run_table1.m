clear; clc;

addpath('../src');

% ==================================================
% Environment warm-up
% ==================================================
fprintf('Initializing Parallel Environment (Warm-up)...\n');
if isempty(gcp('nocreate'))
    parpool(); 
end
fprintf('Parallel Environment Ready.\n\n');
% ==================================================

% Define matrices to test and their corresponding optimal partition numbers (k_partitions)
% Ensure that 'matrices' and 'k_list' have strictly the same length
matrices = {'big.mat', 'ca-CondMat.mat', 'ca-HepPh.mat', 'gre_1107.mat', 'L-9.mat', 'minnesota.mat', 'nopoly.mat', 'pesa.mat', 'Pisa.mat', 'tuma1.mat', 'USpowerGrid.mat', 'usps_norm_5NN.mat', 'vsp_c-30_data_data.mat', 'wing_nodal.mat', };
k_list   = [9, 6, 6, 6, 12, 6, 10, 10, 5, 12, 6, 4, 6, 3];          % Corresponding partition number k for the matrices above

% Expand table header to accommodate the new k column, 4 time comparisons, and 3 error comparisons
fprintf('%-20s | %-6s | %-4s | %-12s | %-35s | %-35s\n', ...
    'Matrix', 'n', 'k', 'Kemeny(P)', 'Time (s): Dir / Rec / DirRec / Par', 'Rel. Error: Rec / DirRec / Par');
fprintf(repmat('-', 1, 126));
fprintf('\n');

for m = 1:length(matrices)
    mat_name = matrices{m};
    k_partitions = k_list(m); % Get the specific partition number k for the current matrix
    
    clear Problem A matrixname G SG bin binsize idx p P pi_vec; 
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
    
    P = spdiags(A*ones(n,1), 0, n, n) \ A;
    [pi_vec, ~] = eigs(P', 1, 'largestabs', 'MaxIterations', 10000);
    pi_vec = pi_vec / sum(pi_vec);

    % 1. Direct Method (Exact inversion / linear system solve)
    tic;
    kemeny_dir = kemenydirect(P);
    t_dir = toc;
    
    % 2. Bini's Recursive Method (Pure GMRES iterative solver)
    tic;
    try
        kemeny_rec = recursivekemeny(P, pi_vec); 
    catch
        kemeny_rec = -1; 
    end
    t_rec = toc;

    % 3. Bini's Dir-Rec Method (Recursion based on exact LU)
    tic;
    try
        kemeny_dirrec = recursivekemenydirect(P, pi_vec); 
    catch
        kemeny_dirrec = -1; 
    end
    t_dirrec = toc;

    % 4. Proposed K-Partitions Method (Multi-process parallel method based on Theorem 3.1)
    options.solver = 'lu'; 
    
    % Pass the specific k_partitions for the current matrix
    [kemeny_par, t_par] = kemeny_k_parallel(P, k_partitions, pi_vec, options);
    
    % Error calculation
    if kemeny_rec ~= -1
        err_rec = abs(kemeny_dir - kemeny_rec) / kemeny_dir;
    else
        err_rec = NaN;
    end
    
    if kemeny_dirrec ~= -1
        err_dirrec = abs(kemeny_dir - kemeny_dirrec) / kemeny_dir;
    else
        err_dirrec = NaN;
    end

    err_par = abs(kemeny_dir - kemeny_par) / kemeny_dir;
    
    % Formatted output (added %-4d for k)
    fprintf('%-20s | %-6d | %-4d | %-12.2f | %-6.2f / %-6.2f / %-6.2f / %-6.2f | %-10.2e / %-10.2e / %-10.2e\n', ...
        strrep(mat_name, '.mat', ''), n, k_partitions, kemeny_dir, t_dir, t_rec, t_dirrec, t_par, err_rec, err_dirrec, err_par);
end