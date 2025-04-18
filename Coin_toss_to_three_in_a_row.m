%% Markov chains - Tossing coin until 3 heads or 3 tails

% States are
% start H  HH HHH   T  TT TTT
POriginal = [
    0 1/2   0   0 1/2   0   0
    0   0 1/2   0 1/2   0   0
    0   0   0 1/2 1/2   0   0
    0   0   0   1   0   0   0
    0 1/2  0   0   0  1/2   0
    0 1/2  0   0   0    0 1/2
    0   0  0   0   0    0   1
    ];

IndexReorder = [4, 7, 1, 2, 3, 5, 6];

P = POriginal(IndexReorder, IndexReorder);

NAbsorbing = 2;
NTransient = 5;

R = P(NAbsorbing+1:end, 1:NAbsorbing);
Q = P(NAbsorbing+1:end, NAbsorbing+1:end);
N = inv(eye(NTransient) - Q);
B = N * R;
S = N * ones([NTransient,1]);