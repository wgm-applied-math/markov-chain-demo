function Trajectory = MCTrajectory(P, StartState, NumSteps)
%MCTrajectory Make a sample trajectory of a Markov Chain
%   Given a matrix of transition probabilities, a start state, and the
%   number of steps, 

% Number of states.
N = size(P, 1);

% Check that the transition matrix is proper.
assert(N == size(P, 2), "The matrix P must be square");
assert(min(P, [], "all") >= 0, "All entries of P must be >= 0");
assert(all(isapprox(sum(P, 2), 1, "veryloose")), "All rows must sum to 1");
assert(1 <= StartState && StartState <= N, ...
    "The start state must be between 1 and the size of P");

CurrentState = StartState;
Trajectory = zeros([1, NumSteps]);
for j = 1:NumSteps
    % One random number between 1 and N,
    % weighted by the row of P corresponding to CurrentState.
    NextState = randsample(N, 1, true, P(CurrentState,:));
    CurrentState = NextState;
    Trajectory(j) = NextState;
end

end