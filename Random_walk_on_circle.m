%% Markov chain - Random walk on a circle

NStates = 4;

P = zeros(NStates,NStates);

for j = 0:NStates-1
    P(1+j,1+mod(j+1,NStates)) = 1/2;
    P(1+j,1+mod(j-1,NStates)) = 1/2;
end
