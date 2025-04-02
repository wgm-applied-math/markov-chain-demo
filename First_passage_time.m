%% First passage time

% Following the inventory example in 28.6
% We have either 0, 1, 2, or 3 cameras in stock
% at the end of the week.
% Demand for cameras has a Poisson distribution
% with rate parameter 1. 
% If demand exhausts the number in stock,
% we just lose those customers and transition to 
% state 0.
% If there are 0 cameras in stock at the end of the week,
% we order 3 cameras to be delivered the next
% week.

% Each of these is probability of demand being for 
% 0, 1, 2, cameras, or 3 or more:
pd0 = exp(-1);
pd1 = exp(-1);
pd2 = exp(-1)/2;
pd3p = 1 - (pd0 + pd1 + pd2);

P = [           pd3p  pd2  pd1  pd0
      (pd1+pd2+pd3p)  pd0    0    0
          (pd2+pd3p)  pd1  pd0    0
                pd3p  pd2  pd1  pd0
    ];
N = size(P, 1);

[V, D, W] = eig(P);
[d, ind] = sort(diag(abs(D)));
w_unscaled = W(:,ind(end))';
w = w_unscaled / sum(w_unscaled);

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