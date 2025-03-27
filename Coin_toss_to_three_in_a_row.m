%% Markov chains - Tossing coin until 3 heads or 3 tails

% States are
% start H  HH HHH   T  TT TTT
P = [
    0 1/2   0   0 1/2   0   0
    0   0 1/2   0 1/2   0   0
    0   0   0 1/2 1/2   0   0
    0   0   0   1   0   0   0
    0 1/2  0   0   0  1/2   0
    0 1/2  0   0   0    0 1/2
    0   0  0   0   0    0   1
    ];