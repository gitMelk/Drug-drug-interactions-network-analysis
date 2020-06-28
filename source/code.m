%#########################################################################
%#########################################################################
%############################   PART I   #################################
%#########################################################################
%#########################################################################

close all
clear all
clc
%% PRE-PROCESSING ########################################################
% Import the file "drugs_edges.txt" created with the Python script 
% "data_gen_DDI.py". 
% The file  contains 2 colums, the source node and the 
% target node, so it represents the edges of the network.
G = importdata('drugs_edges.txt', '\t', 1);

% adjacency matrix
N = max(max(G.data));
A = sparse(G.data(:,2),G.data(:,1),ones(size(G.data,1),1),N,N);

% Simple check to see if there are some self-loops. 
for i=1:N
    if(A(i,i) == 1) 
        disp("self-loop");
    end
end

% build undirected network
A = A+A'; 
clear G; % won't need this G again

% remove nodes which are NOT connected
pos = find(sum(A)~=0);
A = A(pos,pos);
N = size(A,1);

% Find the largest connected component
 e1 = [1;zeros(N-1,1)];
 exit = false;
 while(~exit)
     e1_old = e1;
     e1 = 1*(A*e1+e1>0);
     exit = (sum(e1-e1_old)==0);
 end
 pos = find(e1);
 % The largest component is now the graph on which the analysis will be
 % made. 
 A = A(pos,pos);
 N = size(A,1);

%% LINKS #################################################################
% These are the edges of the network after the pre-processing. 
counter = 1;
for x=1:N
    for y=1:x
        if(A(x,y)== 1)
           links(counter,1) = x;
           links(counter,2) = y;
           counter = counter + 1;
        end
    end
end
links_original = links;
%% DEGREE DISTRIBUTION ###################################################

d = full(sum(A,2)); %in-degree = out-degree
mean_deg = mean(d); 
second_momentum = mean(d.^2);
third_momentum = mean(d.^3);
disp(['Average degree <d>: ' num2str(mean_deg)])
disp(['Second momentum degree <d^2>: ' num2str(second_momentum)])
disp(['Third momentum degree <d^3>: ' num2str(third_momentum)])
d = d(d>0); % avoid zero degrees
k = unique(d); % degree samples
pk = histc(d,k)'; % counts occurrences
pk = pk/sum(pk); % normalize to 1

% Cumulative distribution
Pk = cumsum(pk,'reverse');

% CCDF plot
figure(1)
loglog(k,Pk,'.')
grid on
xlabel('k')
ylabel('CCDF')
title('Logarithmic CCDF plot')

%% GRAPH PROPIETIES ######################################################
% From A we get the new graph to extract its propieties 
G = graph(A); 

L = numedges(G); % # edges
distance = distances(G);
mean_dist = sum(sum(distance))/N^(2);
diameter = max(max(distance));

% Print propieties
disp(['Number of nodes: ' num2str(N)]);
disp(['Number of edges: ' num2str(L)]);
disp(['Diameter: ' num2str(diameter)]);
% Print distances
disp(['Mean distance - real: ' num2str(mean_dist)]);
disp(['Mean distance - random net: ' num2str(log(N)/log(mean_deg))]);
disp(['Mean distance - power law net: ' num2str(log(log(N)))]);


%% ML FITTING ############################################################
kmin = 210; % this is empirically derived from the graph
d_ML = d(d>=kmin); % restrict range
ga = 1+1/mean(log(d_ML/kmin)); % estimate the exponent
disp(['Gamma ML = ' num2str(ga)])

c1 = (ga-1)*kmin^(ga-1); % constant c
C1 = sum(pk(k>=kmin)); % constant C
                       % correction factor taking into account for  
                       % the fractional impact of the degrees k>=kmin

% fitting in the CCDF signal
s3 = C1*c1/(ga-1)*k.^(1-ga);

figure(1)
% Data
hold on
loglog(k,s3);
xlabel('k')
ylabel('CCDF')
title('ML fittings')
legend('data','ML')
ylim([0.0007 1])


%% CLUSTERING ############################################################
E_bet_nodes = diag(A*triu(A)*A); % # of edges between the nodes
C = E_bet_nodes*2./(d.*(d-1)); % Clustering fomula 
% To avoid NaN, it happens when nod deg = 1
C(isnan(C)) = 0;
mean_c = mean(C);
disp(['Average Clustering Coefficient = ' num2str(mean_c)])

figure
plot(d,C,'.')
hold on
plot(d,ones(length(d),1)*mean_c)
grid on
xlabel('k')
ylabel('C(k)')
legend('Clustering Coefficient','Average Clustering Coefficient')

%% ASSORTATIVITY #########################################################
% R-S implementation: Rewire edges, the original degree distribution is 
% manteigned.
% randomperm do not reperat numbers 
links(:,2) = links(randperm(length(links)),2); 

% New adjency matrix, might have sel-floops (not a problem for this task)
A_as = sparse(links(:,2),links(:,1),ones(size(links,1),1),N,N);
A_as = A_as'+A_as;

% Averages deg. of neighbours
k_tmp = (A*d)./d;
k_tmp_as = (A_as*d)./d;

% Extract averages for each value of k
u = unique(d);
for k = 1:length(u)
    knn(k) = mean(k_tmp(d==u(k))); 
    knn_as(k) = mean(k_tmp_as(d==u(k)));
end

% Do the linear fitting
p = polyfit(log(u'),log(knn),1);
pr = polyfit(log(u'),log(knn_as),1);
disp(['Assortativity factor ' num2str(p(1))])
disp(['Assortativity factor random rewiring ' num2str(pr(1))])


figure
loglog(d,k_tmp,'.','Color',[125/255 125/255 125/255]);
hold on
loglog(u,exp(p(2)+log(u)*p(1)),'b-');
%loglog(u,exp(pr(2)+log(u)*pr(1)),'b-'); % Random data(log-bin)
loglog(u,knn,'.','Color',[255/255 0/255 0/255]);
loglog(u,knn_as,'.','Color',[3/255 253/255 86/255]);
hold off
grid
xlabel('k')
ylabel('k_{nn}')
legend('Data','Data(log-bin)','Linear Fitting','R-S')
title('Assortativity')
ylim([3 327])
xlim([0 444])
%% ROBUSTNESS ############################################################
inhomogeneity = second_momentum/mean_deg; % > 2
disp(['Inhomogeneity ratio: ' num2str(inhomogeneity)])

% Random failures
clone_A = graph(A);
P_inf_zero = 1; % all nodes connected 
P_inf_random = zeros(N+1,1);
bincounts = 0;

for i = 1:N
    if (i==1)
         P_inf_random(i) = 1;
    else
    tmp = randi(numnodes(clone_A)); % randomly select a node [1 max]
    clone_A = rmnode(clone_A,tmp); % remove the node
    all_bins = conncomp(clone_A); % find connected subgraphs
    unique_bins = unique(all_bins); % find unique bins
    % then, create the histogram in the range of unique bins
    bincounts = histc(all_bins,unique_bins); 
    P_inf_random(i) = max(bincounts)/N; % # nodes Giant Component over N
    end
end

% Attack failures
clone_A = graph(A); % rebuild the graph
P_inf_atk = zeros(N+1,1);
d_atk = d;
bincounts = 0;

for i = 1:N
    if (i==1)
        P_inf_atk(i) = 1;
    else
    [~,tmp] = max(d_atk); % find highest degree node
    d_atk(tmp) = []; % remove it from the degree list
    clone_A = rmnode(clone_A,tmp); % remove the node
    all_bins = conncomp(clone_A); % find connected subgraphs
    numBins = unique(all_bins); % find unique bins
    % then, create the histogram in the range of unique bins
    bincounts = histc(all_bins,numBins);
    P_inf_atk(i) = max(bincounts)/N; % # nodes Giant Component over N
    end
end
f_del = 0:1/N:1; % Deleted nodes (fraction)
P_random = P_inf_random/P_inf_zero;
P_atk = P_inf_atk/P_inf_zero;

figure
% 100 samples for the plot
plot(f_del(1:N/100:N),P_random(1:N/100:N),'.-', ...
    'Color',[40/255 171/255 31/255]); 
hold on
plot(f_del(1:N/100:N),P_atk(1:N/100:N),'.-', ...
    'Color',[126/255 3/255 189/255]);
hold on
plot(linspace(0,1,2), linspace(0,1,2)*0+0.1, 'r-')
title('Robusteness')
xlabel('f')
ylabel('P_{\infty}(f) / P_{\infty}(0)')
legend('Random Failures','Attack Failures','0.1 thresh')
xlim([0 1])
set(gca,'XTick',(0:0.25:1))
ylim([0 1])
set(gca,'YTick',(0:0.25:1))
%%
%#########################################################################
%#########################################################################
%############################   PART II   ################################
%#########################################################################
%#########################################################################
%%
links = links_original;
%% NAME FILTERING ########################################################

% Import name-nodes from the file created with python.  
% Need to get rid of the names of the deleted nodes.
fid = fopen('drugs_nodes.txt');
% Read all lines & collect in cell array
txt = textscan(fid,'%s','delimiter','	'); 
node_names_raw = txt{1,1}(4:2:end);

counter = 1;
for i=1:size(node_names_raw,1)
    node_names_tmp{i} = node_names_raw(i);
    if (find(links == i)) 
        continue        
    else
        name_to_delete{counter} = node_names_raw(i);
        counter = counter + 1;
    end
end
name_to_delete = name_to_delete';
node_names_tmp = node_names_tmp';
for i=1:size(name_to_delete,1)
    tmp = find(strcmp([node_names_tmp{:}], [name_to_delete{i}]));
    node_names_tmp(tmp) = [];
end

for i=1:size(node_names_tmp,1)
    node_names(i) = node_names_tmp{i};
end
node_names = node_names';
%% DEGREE AND BETWEENNESS RANKING ########################################
G = graph(A);
deg_temp = full(sum(A,2)); %in-degree = out-degree
centrality_temp = centrality(G, 'betweenness');
for i=1:length(centrality_temp)
   centrality_temp(i,2)  = i;
   deg_temp(i,2) = i;
end
sorted_centrality = sortrows(centrality_temp,1,'descend');
sorted_temp = sorted_centrality(:,2);

disp('Top 5 nodes, betweenness:')
disp(node_names(sorted_temp(1:5)))

sorted_deg = sortrows(deg_temp,1,'descend');
sorted_temp = sorted_deg(:,2);

disp('Top 5 nodes, degree:')
disp(node_names(sorted_temp(1:5)))

%clear G  
%% RANKING ###############################################################
%% PAGERANK %%
% A is an undirected graph
M = A*sparse(diag(1./sum(A)));
c = 0.85;
q = ones(N,1)/N;
disp('Computing PageRank - Linear system solution')
tic
rank_pag = sparse((eye(N)-c*M)/(1-c))\q;
rank_pag = rank_pag/sum(rank_pag);
toc

iterations_pag = 35;
disp('Computing PageRank - Power iteration')
tic
speed_pag =[];
% Initial guess: uniform probability distribution
p0 = ones(N,1)/N;
for k = 1:iterations_pag
    % vector's value at k-th iteration
    p0 = c*M*p0+(1-c)*q;
    % Normalization (stochastic)
    p0 = p0/sum(p0);
    % Speed for the convergence
    speed_pag(k) = norm(p0-rank_pag)/sqrt(N);
end
toc

disp('Extracting PageRank Eigenvalues')
tic
% Eigenvector of M
eigs_pag = eigs(M,size(M,1));
toc

% Top 5 nodes with PageRank
[~, rank_pag_top] = sort(rank_pag,'descend');
disp('Top 5 nodes with PageRank:')
disp(node_names(rank_pag_top(1:5)))

% Check if the solutions are the same
[~, rank_pag_top_est] = sort(p0,'descend');
guard = 1;
for i=1:5
   if(rank_pag_top(i) ~= rank_pag_top_est(i)) 
       guard = 0;
   end
end
if(guard == 1)
    disp('Power iteration = Linear solution')
else
    disp('Power iteration != Linear solution')
end

% Results plot
%figure
set(0,'defaultTextInterpreter','latex') % to use LaTeX format
% The biggest eigenvalue = 1 
% The convergence is proportional to the 2nd largest eigenvalue
subplot(1,2,1)
eigenvalue2 = (c*abs(eigs_pag(2))).^(1:iterations_pag);
semilogy([speed_pag; eigenvalue2/eigenvalue2(end)*speed_pag(end)]')
grid
xlim([0 iterations_pag])
legend('Power iteration','Second Eigenvalue')
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('PageRank convergence')

% Plot  M eigenvalues
subplot(1,2,2)
plot(eigs_pag,'x')
hold on
plot(exp(2i*pi*(0:0.001:1)))
grid
xlim([-1 1])
legend('Eigenvalues','Unit circle')
title('PageRank eigenvalues')
%% HITS %%
M = A*A';
disp('Computing HITS - Eigenvalue Extraction')
tic
[pp,ee] = eigs(M,2); % Compute eigenvalues
toc
p = -pp(:,1)/norm(pp(:,1));

iterations_hits = 25;
disp('Computing HITS - Power iteration')
tic
N = size(M,1);
p0 = ones(N,1)/sqrt(N);
speed_hits = [];
for k = 1:iterations_hits
    p00 = p0;
    p0 = A*(A'*p0);
    p0 = p0/norm(p0);
    speed_hits(k) = norm(p0-p00)/sqrt(N);
end
toc 

% Results plot
figure
eig_ratio = (ee(2,2)/ee(1,1)).^(1:iterations_hits);
semilogy([speed_hits;eig_ratio*speed_hits(end)/eig_ratio(end)]')
grid
legend('Power iteration','Second eigenvalue')
xlabel('k [iteration \#]')
ylabel('$\|r_k - r_\infty\|$')
title('HITS Convergence')

% Top 5 nodes with HITS
[~, rank_hits_top] = sort(p0,'descend');
disp('Top 5 nodes with HITS:')
disp(node_names(rank_hits_top(1:5)))

%% Comparison %%
figure
subplot(1,2,1)
stem([p/sum(p),-rank_pag/sum(rank_pag)],'.')
grid
legend('HITS','PageRank')
title('PageRank vs HITS')

subplot(1,2,2)
plot(p/sum(p),rank_pag/sum(rank_pag),'x')
grid
xlabel('HITS score')
ylabel('Pagerank score')
title('PageRank vs HITS ')

%% COMMUNITY DETECTION ###################################################
% SPECTRAL APPROACH
% Build Laplacian matrix
d = full(sum(A)); % Degree vector
Di = spdiags(1./sqrt(d'),0,N,N);
L = speye(N) - Di*A*Di; % Laplacian

% Compute and plot the eigenvalues
eig_lap = eigs(L,N); % eigenvalues
figure
subplot(1,2,1)
plot(eig_lap,'x')
grid
title('Eigenvalues (of the normalized Laplacian)')
% Zoom
subplot(1,2,2)
plot(eig_lap,'x')
axis([1300,1550,0,0.3])
grid
title('Eigenvalues (Zoom of the plot)')


% extract eigenvectors
[V,DD] = eigs(L,3,'SA');
Vv = Di*V; % Normalize eigenvectors
v1 = Vv(:,2)/norm(Vv(:,2)); % Fiedler's vector
v2 = Vv(:,3)/norm(Vv(:,3)); 

% Sweep wrt the ordering identified by v1
% Reorder the adjacency matrix
[v1_sorted, pos] = sort(v1);
A_reordered= A(pos,pos);

% Evaluate the conductance measure
a = sum(triu(A_reordered));
b = sum(tril(A_reordered));
d = a+b;
D = sum(d);
assoc = cumsum(d);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:end-1);

% Plot conductance Measure
figure
plot(conduct,'x-')
grid
ylabel('Conductance $\Phi(A_i)$')
xlabel('Node ordered according to Fiedler''s eigenvector')
title('Conductance')

% There are seven community
% Identify 1 minimum and 5 local minima, so 6 thresholds
min_1 = 135;
min_2 = 171;
min_3 = 279;
min_4 = 561;
min_5 = 930;
min_6 = 1424;


threshold_1 = mean(v1_sorted(min_1:min_1+1));
threshold_2 = mean(v1_sorted(min_2:min_2+1));
threshold_3 = mean(v1_sorted(min_3:min_3+1));
threshold_4 = mean(v1_sorted(min_4:min_4+1));
threshold_5 = mean(v1_sorted(min_5:min_5+1));
threshold_6 = mean(v1_sorted(min_6:min_6+1));

% C1 = sort(pos(1:min_1));
% C2 = sort(pos(min_1+1:min_2));
% C3 = sort(pos(min_2+1:min_3));
% C4 = sort(pos(min_3+1:min_4));
% C5 = sort(pos(min_4+1:min_5));
% C6 = sort(pos(min_5+1:min_6));
% C7 = sort(pos(min_6+1:end));

% figure
% spy(A([C1;C2;C3;C4;C5;C6;C7],[C1;C2;C3;C4;C5;C6;C7]))


% Display data
disp(['Minimum conductance: ' num2str(conduct(min_4))])
disp(['   Cheeger''s lower bound: ' num2str(.5*DD(2,2))])
disp(['   Cheeger''s upper bound: ' num2str(sqrt(2*DD(2,2)))])

disp(['   Community size #1: ' num2str(min_1) '  ' ...
    '(' num2str((min_1)/N*100) '%)'])
disp(['   Community size #2: ' num2str(min_2-min_1) '  ' ...
    '(' num2str((min_2-min_1)/N*100) '%)'])
disp(['   Community size #3: ' num2str(min_3-min_2) '  ' ...
    '(' num2str((min_3-min_2)/N*100) '%)'])
disp(['   Community size #4: ' num2str(min_4-min_3) '  ' ...
    '(' num2str((min_4-min_3)/N*100) '%)'])
disp(['   Community size #5: ' num2str(min_5-min_4) '  ' ...
    '(' num2str((min_5-min_4)/N*100) '%)'])
disp(['   Community size #6: ' num2str(min_6-min_5) '  ' ...
    '(' num2str((min_6-min_5)/N*100) '%)'])
disp(['   Community size #7: ' num2str(N-min_6) '  ' ...
    '(' num2str((N-min_6)/N*100) '%)'])

% Show network with partition
figure
plot(v1,v2,'.')
grid
hold on
plot(threshold_1*[1,1],ylim,'r-')
plot(threshold_2*[1,1],ylim,'r-')
plot(threshold_3*[1,1],ylim,'r-')
plot(threshold_4*[1,1],ylim,'r-')
plot(threshold_5*[1,1],ylim,'r-')
plot(threshold_6*[1,1],ylim,'r-')
xlabel('$v_{N-1}$')
ylabel('$v_{N-2}$')
ylim([-0.2 0.2])
xlim([-15e-4 10e-4])
title('2D representation by Laplacian eigenvectors')

%% LINK PREDICTION #######################################################
%close all

AUC_all = zeros(5,1); % initialize the AUC result vector 
precision_all = zeros(5,1); % initialize the precision result vector
L = N*0.1; % Top 10% links, for the precision

% Create the probe and test sets
G = graph(A);
all_links = numedges(G);
test_size = 0.8;
probe_size = round(all_links*(1-test_size));

for i = 1:probe_size
    edge_remove = randi(all_links); % Suppress a random link
    G = rmedge(G, edge_remove);
    all_links = numedges(G);
end 

test_set = full(adjacency(G));
probe_set = A - test_set;
results_counter = 1;

%% SIMILARITY methods

% Neighbour based 
% Common Neighbours - CM
s_matix_CN = test_set*test_set;

AUC_all(results_counter) = compute_AUC(probe_set, A, s_matix_CN);
precision_all(results_counter) = compute_precision(probe_set, ...
    test_set, s_matix_CN, L);
results_counter = results_counter + 1;

% Adamic Adar - AA
deg = sum(test_set,1); % Degree in row
weight = 1./log(deg); % Weight of each node
W = repmat(weight,N,1);
A_weighted = test_set.*W; % Apply weight to the matrix
A_weighted(isnan(A_weighted)) = 0;
s_matrix_AA = A_weighted*test_set; % Similarity matrix

AUC_all(results_counter) = compute_AUC(probe_set, A, s_matrix_AA);
precision_all(results_counter) = compute_precision(probe_set, ...
    test_set, s_matrix_AA, L);
results_counter = results_counter + 1;

% Resource Allocation - RA
weight = 1./deg; % Weight of each node
W = repmat(weight,N,1);
A_weighted = test_set.*W; % Apply weight to the matrix
A_weighted(isnan(A_weighted)) = 0;
s_matrix_RA = A_weighted*test_set; % Similarity matrix

AUC_all(results_counter) = compute_AUC(probe_set, A, s_matrix_RA);
precision_all(results_counter) = compute_precision(probe_set, ...
    test_set, s_matrix_RA, L);
results_counter = results_counter + 1;


% Path Based

% Local Path - LP
beta = 0.5;
s_matrix_LP = test_set^2 + beta*test_set^3; % Similarity matrix 

AUC_all(results_counter) = compute_AUC(probe_set, A, s_matrix_LP);
precision_all(results_counter) = compute_precision(probe_set, ...
    test_set, s_matrix_LP, L);
results_counter = results_counter + 1;

% Kats
alpha = 0.2;
I = eye(N);    
s_matrix_Kats = (I-alpha*A)^-1 - I; % Similarity matrix

AUC_all(results_counter) = compute_AUC(probe_set, A, s_matrix_Kats);
precision_all(results_counter) = compute_precision(probe_set, ...
     test_set, s_matrix_Kats, L);
results_counter = results_counter + 1;


% Random walk based 

% Random Walk with Restart - RWR 
weight = 1./deg; % Weight of each node
W = repmat(weight,N,1);
M = test_set.*W; % Matrix with weights for the PageRank
M(isnan(M)) = 0;
c = 0.85;

I = eye(N);                                                        
iteration_num = 35; 
s_matrix_RWR = zeros(N);  

% Power iteration
for k=1:iteration_num
	s_matrix_RWR =  c*M*s_matrix_RWR + (1-c)*I;
    s_matrix_RWR(isnan(s_matrix_RWR)) = 0; % remove NaN
    % Normalize the Matrix (it is composed by stochastic vectors)
    s_matrix_RWR = s_matrix_RWR/sum(s_matrix_RWR,1);
end

s_matrix_RWR = s_matrix_RWR + s_matrix_RWR'; % Similarity matrix

AUC_all(results_counter) = compute_AUC(probe_set, A, s_matrix_RWR);
precision_all(results_counter) = compute_precision(probe_set, ...
    test_set, s_matrix_RWR, L);


% Display results
disp_names = {(' CN '), (' AA '), (' RA '), (' LP '), ...
    (' Kats '), (' RWR ')};
disp(disp_names);
disp('AUC:')
disp(AUC_all');
disp('Precision:')
disp(precision_all');
%% Auxiliary functions ###################################################
function AUC = compute_AUC(probe_set, A, S)
    % S is symmetric -> use the upper triangle of the matrix
    A = A + eye(size(A,1));
    A_v = A(triu(true(size(A))));
    S_v = S(triu(true(size(S))));
    probe_vec = probe_set(triu(true(size(probe_set))));
    
    S_active = S_v(probe_vec==1); % Similarity active edges (probe set)
    S_inactive = S_v(A_v == 0); % Similarity inactive edges (probe set)
    
    % Compute auc
    AUC = 0;
    for A=1:length(S_active)
        for i=1:length(S_inactive)
            if (S_active(A) > S_inactive(i))
                AUC = AUC + 1;
            end
        end
    end
    AUC = AUC/(size(S_active,1)*size(S_inactive,1));
end

function precision = compute_precision(probe_set, test_set, S, L)
    % S is symmetric -> use the upper triangle of the matrix
    S_v = S(triu(true(size(S)), -1));
    probe_v= probe_set(triu(true(size(probe_set)), -1));
    test_v = test_set(triu(true(size(test_set)), -1));
    
    S_v = S_v(test_v == 0); % Similarity inactive links (test set)
    probe_v = probe_v(test_v == 0);
    
    [~, idx] = sort(S_v,'descend'); % Sort S_vec
    probe_v = probe_v(idx); % Coherenty sort the probe
    % Compute precision
    L = min([L,length(probe_v)]);
    precision = sum(probe_v(1:L))/L;
end
