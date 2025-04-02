%% Markov chain - Machine degredation and repair

% Have a machine that can be in one of four states:
% 1 - like new
% 2 - broken in
% 3 - damaged
% 4 - unusable

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