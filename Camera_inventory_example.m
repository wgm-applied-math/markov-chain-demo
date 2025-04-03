%% Camera inventory example

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
M = N - 1;

%% Finding the stationary distribution as an eigenvector

[V, D, W] = eig(P);
[d, ind] = sort(diag(abs(D)));
w_unscaled = W(:,ind(end))';
w = w_unscaled / sum(w_unscaled);

%% Expected cost

% Demand each week has a Poisson distribution
% with rate parameter 1,
% so Prob(D = k) = exp(-1) / k!

% The cost to order z cameras is 10 + 25 z.
% The cost of each disappointed customer is 50.

% The actual cost of next week given the state at the end of this week:
function c = Cost(M, x, dNext)
    if x == 0
        % At the end of the week we have 0 in stock.
        % We order M cameras.
        % The cost for the following week includes that
        % plus any shortage cost.
        c = 10 + 25*M + 50*max([dNext-M, 0]);
    else
        % At the end of the week we have x in stock.
        % We don't place an order.
        % The cost for the following week is just shortage.
        c = 50*max([dNext-x, 0])
    end
end

% The average cost of next week given the state at the end of this week:
function c = ExCost(M, x)
    if x == 0
        % Ordering cost
        c = 10 + 25*M;
        % Approximate mean shortage cost
        for d = (1+M):20
            c = c + 50*(d - M)*exp(-1)/factorial(d);
        end
    else
        % No order.
        c = 0;
        % Approximate mean shortage cost
        for d = (x+1):20
            c = c + 50*(d - x)*exp(-1)/factorial(d);
        end
    end
end

% Then the expected cost per week in the long term is
cLongTerm = 0;
for x = 0:M
    cLongTerm = cLongTerm + w(1+x)*ExCost(M, x);
end

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