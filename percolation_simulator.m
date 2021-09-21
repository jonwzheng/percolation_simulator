%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  authored by Jonathan Zheng                              %
%  jonzheng@mit.edu    |    December 9th, 2020             %
%                                                          %
%  This code simulates percolation on a 2-D lattice.       %
%  A given lattice of dimension length N is instantiated,  %
%  then xdead (% fraction of defects) * N of those sites   %
%  are labeled as defects. All valid bond sites            %
%  (those connected to non-defects and non-edges) have     %
%  probability p related by "prob_factor" to form.         %
%  Bond info is stored in an N^2 x N^2 adjacency matrix,   %
%  which is converted to a graph; from this graph, the     %
%  script automatically analyzes whether percolation has   %
%  occurred. This whole process repeats Nsamples times.    %
%  Also included: functions to measure the mean cluster    %
%  size & chance of belonging to an infinite cluster.      %
%
%  N.B. May be able to improve computational efficiency    %
%  by recasting bonds and atoms in a sparse matrix.        %

% ******************************************************** %
%% ************** USER FLAGS (change these) ************* %%
N = 20; % # sites per row and column (total N^2 sites)
Nsamples = 5; % # of times to randomly create and test matrix
xdead = 0.03; % fraction of dead proteins.
flags.vis = 0; % visualize each trial lattice                 (1 on, 0 off)
flags.cluster_analysis = 0; %  do cluster & bond analysis     (1 on, 0 off)
flags.verb = 1; % command window messages indicating progress (1 on, 0 off).
flags.ignore_ones = 1; % ignore nodes of size 1 in cluster size  
                       %    determination: highly recommended (1 on, 0 off)
flags.probtype = 1;    % choose whether to calculate probability of 
                       %    bond formation using Boltzmann probability
                       %    or explicitly-given probability   (1 explicit,
                       %                                       0 Boltzmann)

%% Types of User Simulations
%% Example 1: Sweep over & sample a probability range
%sample_range = [0 0.2 0.4 (1-logspace(-0.25,-0.3,20)) logspace(-0.3,-0.25, 20) .6 .8 1];
%sample_range = [0 0.4 0.7 0.72 :0.02:1] ;
sample_range = linspace(0,1,11);
%sample_range = logspace( log10(0.5), log10(0.51), 16);
percolate_matrix = zeros(length(sample_range),Nsamples);
pinf_matrix = zeros(length(sample_range),Nsamples);
clusters = cell(1, length(sample_range));

for E = 1:length(sample_range)
    [percolate_matrix(E,:), clusters{E}, pinf_matrix(E,:)] = ...
        percolation_mc(sample_range(E), N, Nsamples, xdead, flags);
end

if flags.probtype == 0
    p_vector = exp(-sample_range);
elseif flags.probtype == 1
    p_vector = sample_range;
end

% percolation analysis
mean_percolation = mean(percolate_matrix, 2);
figure(1)
plot(p_vector, mean_percolation,'Color','b','LineWidth',3)
xlabel("Probability of bond formation"); ylabel("% trials percolated"); grid on;
title(['Percolation threshold (lattice length = ' num2str(N) ', x_{dead} = ' num2str(xdead)...
    ', n_{samples} = ' num2str(Nsamples) ')'])
xlim([0 1]); ylim([0 1])
hold on;
plot([0.5 0.5], [0 1], '--','Color','k','LineWidth',2);
set(gca,'FontSize',12)

% avg. cluster size analysis
avg_sizes = avg_cluster(clusters, Nsamples, N);
figure(10)
scatter(p_vector, avg_sizes, 'filled')
xlabel("P(bond formation)"); ylabel("Avg. cluster size"); grid on;
title(['Avg. cluster sizes (dim length = ' num2str(N) ', xdead = ' num2str(xdead)...
    ', nsamples = ' num2str(Nsamples) ')'])

% order parameter analysis: 
alpha = power_fit(p_vector, avg_sizes, 0.5);
% infinite bond analysis

pinf_avgs = zeros(length(sample_range),1);
for a = 1:length(sample_range) 
    b = pinf_matrix(a,:);  
    c = b(b ~= -1); % eliminate non-percolated samples.
    if ~isempty(c)
        pinf_avgs(a) = mean(c);
    end
end

beta = power_fit(p_vector, pinf_avgs, 0.5);


plot(p_vector, pinf_avgs,'-o','Color','b','LineWidth',3,'MarkerSize',4)
xlabel("P(bond formation)"); ylabel("Chance of belonging to infinite cluster"); grid on;
title(['p_{infinite} (dim length = ' num2str(N) ', x_{dead} = ' num2str(xdead)...
        ', n_{samples} = ' num2str(Nsamples) ')'])
xlim([0,1])
hold on;
plot([0.5 0.5], [0 1], '--','Color','k','LineWidth',2);
set(gca,'FontSize',12)


%% Example 2: Visualize a large system with bound proteins
N = 30;
flags.vis = 1; %
Nsamples = 1;
x_bound = 0.01;
percolation_mc(.7, N, 1, x_bound, flags);

disp("Simulation complete.")

%% Simulation Functions
function [percolated, clusters, fracinf] = percolation_mc(prob_factor, N, nsamples, xdead, flags)
% Input
%   prob_factor: either explicitly given probability of bond forming or
%               ratio of E to T used for Boltzmann bond formation,
%               choice dictated by flags.probtype [1 x 1 double]
%   N:          length of one direction of lattice
%   nsamples:   number of MC samples to take
%   xdead:     proportion of non-denatured sites (Ratio) [1 x 1 double]
%   flags:      user input for certain options: 
%                   vis - data visualization (1 on, 0 off)
%                   cluster_analysis - whether to analyze (1 on, 0 off)
%                   verbosity - whether to show a lot of output in
%                                       the terminal (1 on, 0 off)
%                   probtype - type of probability supplied to function
%                                       (1 explicit, 0 Boltzmann)
%
% Output
%   percolated: vector of whether each trial is percolated  [nsamples x 1]
%   clusters:   cell array containing size of clusters in simulation.
%   fracinf:    vector containing % chance of being in infinite cluster

    % Unpack parameters
    disp("Instantiating simulation...")
    vis = flags.vis;
    cluster_analysis = flags.cluster_analysis;
    verb = flags.verb;

    % Acceptance criterion: based on Boltzmann energy
    if flags.probtype == 0
        Paccept = exp(-prob_factor);
    elseif flags.probtype == 1
        Paccept = prob_factor;
    else
        error("Invalid flags.probtype. Choose either 0 (Boltzmann) or 1 (explicit)")
    end
    if Paccept < 0 || Paccept > 1
        error("Probability of bond formation must be between 0 or 1.")
    end
    
    disp(['P(accept): ' num2str(Paccept)])
    disp(['Lattice size: ' num2str(N) ' by ' num2str(N)])
    disp(['Total size: ' num2str(N.^2)])
    disp("*****************************************")

    fracinf = ones(1,nsamples) .* -1;
    percolated = ones(1, nsamples);
    clusters = cell(nsamples,1);
    
    % Loop through samples
    for t = 1:nsamples
        if verb 
            disp(['*** Trial ' num2str(t) ' ***'])
            disp(['Creating lattice of size ' num2str(N) ' by ' num2str(N) '...' ])
        end
        % Create lattice of sites as N x N matrix, introducing unbound sites
        sites = ones(1,N.^2);
        nbound = floor(N.^2 * xdead); 
        l = 0;
        while l < nbound % ensure unbound sites are all introduced
            randrow = randi(N);
            randcol = randi(N);
            n = atom2bond(randrow, randcol, N);
            if sites(n) == 1
                sites(n) = 0;
                l = l + 1;
            end
        end 
        if verb
            disp([num2str(l) ' bound sites were created.'])
        end
        % Create lattice of bonds as N^2 x N^2 matrix
        bonds = eye(N^2); % row element's interaction on column element 
                          % may be improvable via sparsity
        % Algorithm: Iterate through all atoms
        %            For each atom, iterate through bottom + right bonds
        %               (avoid double counting)
        %            Post-process the matrix at the end to ensure that
        %               the matrix is upper diagonal.
        for j = 1 : N % iterating through y positions
            for i = 1 : N % iterating through x positions
                n = atom2bond(i,j,N);
                bond_ind = allowable_bonds(i, j, N, sites);                
                for k = 1:length(bond_ind)
                    if sites(n) == 1 && sites(bond_ind(k)) == 1 
                        bonds(n, bond_ind(k)) = double( rand(1) < Paccept);
                    end
                end
                
            % end of column index loop.
            end
                        
        % end of row index loop.
        end
        
        if verb; disp( ['Finished sampling trial number ' num2str(t)] ); end
        
        % Post-processing
        bonds_processed = diag_A(bonds);
        
        % see whether system is percolated
        percolated(t) = check_percolated(bonds_processed, N, flags); 

        % get cluster size
        if cluster_analysis 
            clusters{t} = get_clusters(bonds_processed, flags);
        end
        
        % get p of being in an infinite cluster if system is percolated
        if percolated(t) && cluster_analysis
            fracinf(t) = strength_analysis(bonds, flags);
        end
        
        if verb; disp("   "); end
    % end of sampling loop.
    end
    disp("Done with all trials.")
    disp("*****************************************")    
    
    if vis; datavis(sites, bonds_processed, N); end 
    
    disp("Done!")
    disp(" ")
end

%% Secondary functions

function neighboring_posns  = allowable_bonds(i, j, N, sites)
% algorithm for getting allowable bonds based on lattice position
% Input
%   i: x coordinate
%   j: y coordinate
%   N: size of lattice in 1dim
%   atoms: array of atoms indicating whether they're bound or unbound
% Output
%   neighboring_posns: 2 x 1 list of bonds for a given site

    neighboring_posns = zeros(1,2);
    allowed_mat = [1, 1];

    % Check if neighboring atom is even allowed or not (bound/unbound)
    % - only counts right & bottom neighbors( avoids double counting)
    
    if i < N
        if sites(atom2bond(i+1,j,N) ) == 0 % check to the right
            allowed_mat(2) = 0; %
        end
        imat = [i, i+1]; % just check the right
    else
        imat = [i, i];
        allowed_mat(2) = 0;
    end
            
    if j < N
        if sites(atom2bond(i,j+1,N) ) == 0
            allowed_mat(1) = 0; % check below; with 4 dims, index = 2
        end 
        jmat = [j+1, j]; % just check the one below.
    else
        jmat = [j, j];
        allowed_mat(1) = 0;
    end 
    for z = 1:2 % this is a dumb way of doing it, but it works .. probably
        if allowed_mat(z) == 1
            neighboring_posns(z) = atom2bond(imat(z), jmat(z), N);
        end 
    end
    neighboring_posns = nonzeros(neighboring_posns);
end

function node_sizes = get_clusters(bonds, flags)
    % Warning: May not give accurate results for bound proteins.
    A = graph(bonds,'upper');

    % idea: do depth first search. get number of edges
    %   this way. keep track of nodes traversed in search.
    %   repeat, until no more nodes are remaining

    n_nodes = length(bonds);
    all_nodes = 1 : 1 : n_nodes; % enumerate all nodes out
    traversed_nodes = [];
    node_sizes = [];
    i = 1;
    while ~isempty(all_nodes)
        searched_nodes = dfsearch(A,i);
        if length(searched_nodes) == 1 && flags.ignore_ones
        else
            node_sizes = [node_sizes length(searched_nodes)];
        end
        traversed_nodes = [traversed_nodes; searched_nodes];
        all_nodes = setdiff(all_nodes, traversed_nodes);
        i = min(all_nodes);
    end
    
    if flags.verb; disp("Finished cluster size analysis."); end
end

function is_percolated = check_percolated(bonds, N, flags)
    A = graph(bonds,'upper');
    if flags.vis
        figure(3)
        plot(A)
    end
    
    % basic algorithm is to see if there's actually a finite path 
    %   from ANY point on 1 edge to the other edge.
    % To be fully percolated, this must span from x = 1 to x = N,
    %   note this does not consider, say, left edge to the top edge,
    % Check vertical edges first.
    i = 1;
    paths = [];
    while i <= N && isempty(paths)
        n = atom2bond(i, 1, N);
        i2 = 1;
        while i2 <= N && isempty(paths)
            n2 = atom2bond(i2, N, N);
            i2 = i2 + 1;
            paths = shortestpath(A, n, n2); % returns [] if no path exists
        end
        i = i + 1;
    end
    
    j = 1;
    % Check horizontal edges.
    while j <= N && isempty(paths)
        n = atom2bond(1, j, N);
        j2 = 1;
        while j2 <= N && isempty(paths)
            n2 = atom2bond(N, j2, N);
            j2 = i2 + 1;
            paths = shortestpath(A, n, n2);
        end
        j = i + 1;
    end

    is_percolated = ~isempty(paths);
    if flags.verb
        if(is_percolated); disp("* Percolated system! *")
        else; disp("No percolation."); end
    end
end        


function size_infinite_clusters = get_all_infinite_clusters(bonds, flags)
    A = graph(bonds,'upper');
   
    N = sqrt(length(bonds));
    % Check horizontal edges first.
    exhausted_start = 0;
    size_infinite_clusters = [];
    infinite_clusters = [];
    
    untraversed_start = 1 : N : N*(N-1)+1;
    untraversed_end = N : N : N^2;
    while ~ exhausted_start    
        exhausted_edge = 0;
        n = min(untraversed_start); % pick the lowest start index not yet traversed.
        locally_untraversed_end = untraversed_end; 
        if ~isempty(n)        
            while ~ exhausted_edge
                n2 = min(locally_untraversed_end);
                if ~isempty(n2)
                    nodesearch = shortestpath(A, n, n2);
                    if ~isempty(nodesearch) % found an infinite path!
                        atoms_to_rem = dfsearch(A, n);
                        untraversed_start = setdiff(untraversed_start, atoms_to_rem);
                        untraversed_end = setdiff(untraversed_end, atoms_to_rem);
                        locally_untraversed_end = setdiff(locally_untraversed_end, atoms_to_rem);
                        infinite_clusters = [infinite_clusters; atoms_to_rem];
                        size_infinite_clusters = [size_infinite_clusters length(atoms_to_rem)];
                    else
                        locally_untraversed_end = setdiff(locally_untraversed_end, n2);
                    end
                else
                    exhausted_edge = 1;
                end
            end 
            untraversed_start = setdiff(untraversed_start, n);
        else
            exhausted_start = 1;
        end
    end 
    
    % Repeat for vertical edges, with more info than before!
    exhausted_start = 0;
    
    untraversed_start = setdiff(1 : 1 : N, infinite_clusters);
    untraversed_end = setdiff(N*(N-1)+1 : 1 : N^2 , infinite_clusters);
    
    while ~ exhausted_start    
        exhausted_edge = 0;
        n = min(untraversed_start); % pick the lowest start index not yet traversed.
        locally_untraversed_end = untraversed_end; 
        if ~isempty(n)        
            while ~ exhausted_edge
                n2 = min(locally_untraversed_end);
                if ~isempty(n2)
                    nodesearch = shortestpath(A, n, n2);
                    if ~isempty(nodesearch) % found a path!
                        atoms_to_rem = dfsearch(A, n);
                        untraversed_start = setdiff(untraversed_start, atoms_to_rem);
                        untraversed_end = setdiff(untraversed_end, atoms_to_rem);
                        locally_untraversed_end = setdiff(locally_untraversed_end, atoms_to_rem);
                        infinite_clusters = [infinite_clusters; atoms_to_rem];
                        size_infinite_clusters = [size_infinite_clusters length(atoms_to_rem)];
                    else
                        locally_untraversed_end = setdiff(locally_untraversed_end, n2);
                    end
                else
                    exhausted_edge = 1;
                end
            end 
            untraversed_start = setdiff(untraversed_start, n);
        else
            exhausted_start = 1;
        end
    end 
    if flags.verb
        disp(['Number of infinite clusters: ' num2str(length(size_infinite_clusters))]);
    end
end        

function alpha = power_fit(pvector, sizes, pcrit)
    % Based on provided pcritical values, attempts to 
    % find scaling behavior near critical point.
    p = abs(pcrit - pvector);
    coeffs = polyfit(log(p), log(sizes), 1);
    alpha = coeffs;
end

function frac_inf = strength_analysis(bonds, flags)
    % Strength def: fraction of any given site belonging to an infinite
    %       cluster. 
    Nbonds = length(bonds);
    frac_inf = sum(get_all_infinite_clusters(bonds, flags)) ./ Nbonds;    
end 

function datavis(sites, bonds, N)
    % Visualizes a network using its (1) lattice modeland (2) graph model
    disp("Visualizing sites...")
    figure(2)
    lattice = zeros(N,N);
    for row = 1:N
        lattice(row,:) = sites(1 + (row-1)*N : row*N);
    end 
    [X, Y] = meshgrid(1:N, 1:N);
    
    for i = 1:N
        color = 0.2+0.1*lattice(N+1-i,:);
        scatter(X(i,:),Y(i,:), [], color, 'filled');
        hold on;
    end
    colormap copper % default is also OK
    xlim([0 N+1]); ylim([0 N+1]);
    
    disp("")
    disp("Visualizing bonds...")
    for p = 1 : N^2 % iterate through all the rows of "bonds", essentially
        [xstart, ystart] = bond2atom(p,N);
        % process each node this point is connected to
        end_indices = find(bonds(p,:));
        for ind = 1:length(end_indices)
            [x,y] = bond2atom(end_indices(ind),N);
            if ((y - ystart) >= (N-1)) % check y overlap
                y = ystart - 1;
            elseif (ystart - y >= (N-1))
                y = ystart + 1;
            end
            
            if (x - xstart >= (N-1))
                x = xstart - 1;
            elseif (xstart - x >= (N-1))
                x = xstart + 1;
            end
            % bandaid fix...
            plot([xstart, x], N-[ystart, y]+1,'b','LineWidth',1.5)
        end
        % end of finding all nn
    end
    disp("Done visualizing.")    
    disp("*****************************************")    

end

%% Helper functions

function avgsizes = avg_cluster(clusters, Nsamples, N)
% reads in a cell of arrays of cluster sizes and converts to a vector of
%       cluster sizes, [1 x trials]
    M = length(clusters);
    avgsizes = ones(M, 1);
    for i = 1 : M
        cell_sample = ones(1,Nsamples);
        obj_to_extract = clusters{i};
        if isempty(obj_to_extract) % if only nodes of size 1
            cell_sample(i) = 0;
        else
            for sample = 1:Nsamples
                cell_sample(sample) = mean(cell2mat(obj_to_extract(sample)));
            end 
        end 
        avgsizes(i) = mean(cell_sample);
    end
    avgsizes = avgsizes ./ N.^2 ;
end


function [i, j] = bond2atom(n,N)
% convert n index to i,j indices.
    i = mod(n,N);
    if i == 0
        i = N;
    end
    j = ceil(n / N);
end 

function n = atom2bond(i,j,N)
% Convert i,j indices to n index. No PBC.
    n = i + N*(j-1);
end

function diagA = diag_A(A)
% converts a nonsymm matrix into a diagonal matrix minus its 1 indices
    size_A = length(A);
    low_A = tril(A);
    upA = triu(A);
    diagA = upA + (low_A - 2*eye(size_A))' ;
end
    

%% Obsolete functions
function n = atom2bond_PBC(i,j,N)
% Convert i,j indices to n index
% Works with PBC
    if j == 0 % above top edge
        %check corner cases
        if i == 0 % top left corner
            n = N.^2;
        elseif i == N+1 % top right corner
            n = N*(N-1) + 1;
        else
            n = i + N*(N-1);
        end 
    elseif j == N+1 % below bottom edge
        if i == 0 % bottom left corner
            n = N;
        elseif i == N+1 % bottom right corner
            n = 1;
        else
            n = i;
        end 
    elseif i == 0 % left of left edge
        n = N + N*(j-1);
    elseif i == (N+1) % right of right edge
        n = 1 + N*(j-1);
    else
        n = i + N*(j-1);
    end 
end
