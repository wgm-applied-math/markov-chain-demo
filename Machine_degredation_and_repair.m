%% Markov chain - Machine degredation and repair

% Have a machine that can be in one of four states:
% 1 - like new
% 2 - broken in
% 3 - damaged
% 4 - unusable

StateNames = ["like new", "broken in", "damaged", "unusable"];

% General wear and tear causes transitions 1 -> 2 -> 3 -> 4.
% A damaged machine (3) is sometimes selected for repair, 
% which returns it to (2) broken in.
% An unusable machine (4) will frequently be refurbished,
% which returns it to (1) like new.

P = [
    0.0 1.0 0.0 0.0
    0.0 0.9 0.1 0.0
    0.0 0.1 0.8 0.1
    0.9 0.0 0.0 0.1
];
N = 4;

%% Stationary distribution

% Since this MC is regular, it has a unique stationary
% distribution, which is a left-eigenvector for eigenvalue 1.
% The eig function computes right-eigenvectors (columns of V),
% left-eigenvectors (columns of W, *shrug*).
% The diagonal matrix D contains the eigenvalues.
[V, D, W] = eig(P);

% Get the eigenvalues in order of increasing magnitude:
[d, ind] = sort(diag(abs(D)));

% Now ind(end) is the index of the largest eigenvalue,
% which by the Perron-Frobenius theorem is 1 for a
% stochastic matrix.

% Pull out the corresponding left-eigenvector.
% This is a column of W, even though a left-eigenvector
% must be a row.  Hence the transpose here.
w_unscaled = W(:,ind(end))';

% Scale it to a distribution vector
w = w_unscaled / sum(w_unscaled);

%% Finding all expected first passage times

% Vector of first return times
tau_d = 1.0 ./ w;

% Matrix of copies of rows like w
A = ones([N,1]) * w;

% Fundamental matrix
Z = inv(eye(N) - P + A);

% Zero out the non-diagonal elements of Z.
% diag(square matrix) = its diagonal as a vector
% diag(vector) = square matrix of 0s but with the vector entries
% going down the diagonal
% Which leads to this bizarre formula:
Zd = diag(diag(Z));

% Matrix of first passage times
tau = (eye(N) - Z + ones([N,N]) * Zd) * diag(tau_d);

% This should be equal to tau
tau_check = ones([N,N]) + P * (tau - diag(tau_d));

% And tau should have tau_d on the diagonal.

%% Sample trajectory

StartState = 1;
NumSteps = 50;

Trajectory = MCTrajectory(P, StartState, NumSteps);

% Line plot

fig = figure();
t = tiledlayout(1,1);
ax = nexttile(t);

p = plot(ax, Trajectory, "x-");

% Make sure there's room for all states vertically
ylim(ax, [0.5, N+0.5]);

% Histogram of states

fig = figure();
t = tiledlayout(1,1);
ax = nexttile(t);

h = histogram(ax, Trajectory, BinEdges=0.5:(N+0.5));

% Pie chart

fig = figure();

% The histogram function produces bin counts
counts = h.BinCounts;

% Has to be a fig, not axes. Shrug.
piechart(fig, counts, StateNames);


%% Picture of the network of states

% Construct a directed graph. This function call interprets P as a weighted
% adjacency matrix, which is perfect for a Markov chain. We can also name
% the states with a vector of strings.
G = digraph(P, StateNames);

fig = figure();
t = tiledlayout(1,1);
ax = nexttile(t);

plot(ax, G, EdgeLabel=G.Edges.Weight, MarkerSize=10);