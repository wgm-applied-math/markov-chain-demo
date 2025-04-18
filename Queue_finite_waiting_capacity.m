%% Queue with finite waiting capacity

% Suppose a queue has room for m customers to wait.
% At each time step,
% the probability that a customer arrives is pA
% and 
% the probability that a customer departs is pD.
% There is at most one event,
% an arrival or departure
% during each time step.
% The state is the number of customers waiting.
% If the waiting area is full and someone else
% arrives, the state is unchanged.

pA = 0.08;
pD = 0.09;

% m = 3, so four states 0, 1, 2, 3
PRegular = [
    1-pA       pA        0     0
      pD  1-pA-pD       pA     0
       0       pD  1-pA-pD    pA
       0        0       pD  1-pD
    ];

%% Expected time to full

% To calculate the expected time until
% the waiting area is full (state 3, index 4)
% convert state 3 to an absorbing state
% and use the fundamental matrix.

PAbsOriginal = PRegular;
PAbsOriginal(4,:) = [0 0 0 1];

% Re-order states
IndexReorder = [4, 1, 2, 3];

P = PAbsOriginal(IndexReorder, IndexReorder);

NAbsorbing = 1;
NTransient = 3;

R = P(NAbsorbing+1:end, 1:NAbsorbing);
Q = P(NAbsorbing+1:end, NAbsorbing+1:end);
N = inv(eye(NTransient) - Q);
B = N * R;
S = N * ones([NTransient,1]);