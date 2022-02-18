function [bestSoFarFit ,bestSoFarSolution ...
    ]=simpleEA( ...  % name of your simple EA function
    fitFunc, ... % name of objective/fitness function
    T, ... % total number of evaluations
    input) % replace it by your input arguments

% Check the inputs
if isempty(fitFunc)
  warning(['Objective function not specified, ''' objFunc ''' used']);
  fitFunc = 'objFunc';
end
if ~ischar(fitFunc)
  error('Argument FITFUNC must be a string');
end
if isempty(T)
  warning(['Budget not specified. 1000000 used']);
  T = '1000000';
end
eval(sprintf('objective=@%s;',fitFunc));
% Initialise variables
nbGen = 0; % generation counter
nbEval = 0; % evaluation counter
bestSoFarFit = 0; % best-so-far fitness value
bestSoFarSolution = NaN; % best-so-far solution
%recorders
fitness_gen=[]; % record the best fitness so far
solution_gen=[];% record the best phenotype of each generation
fitness_pop=[];% record the best fitness in current population 
%% Below starting your code
pop_size=4;
LowerBound=0;
UpperBound=31;
GeneLen=5;
Populations=randi([LowerBound,UpperBound],pop_size,1);
genotypes=dec2bin(Populations);

% Initialise a population
%% TODO
fitness=objective(Populations);
[Sort_array,index]=sort(fitness,'descend');
fitness_pop=[fitness_pop,Sort_array(1)];
for i=1:pop_size
    if fitness(i)>bestSoFarFit
        bestSoFarFit=fitness(i);
        bestSoFarSolution=Populations(i,:);
    end
end
nbGen=nbGen+1;
nbEval=nbEval+pop_size;
fitness_gen=[fitness_gen,bestSoFarFit];
solution_gen=[solution_gen,bestSoFarSolution];
% Evaluate the initial population
%% TODO
% Start the loop
while (nbEval<T) 
% Reproduction (selection, crossver)
cross_prob=fitness./sum(fitness);
offspringGenes=[];
for i=1:pop_size/2
    parentIndexes=[];
    for j=1:2
        rnd=rand();
        for Index=1:pop_size
            if rnd>sum(cross_prob(1:Index-1))&&rnd<=sum(cross_prob(1:Index))
                break;
            end
        end
        parentIndexes=[parentIndexes,Index];
    end
    crossPoint=randi([1,GeneLen-1]);
    offspringGenes=[offspringGenes;[genotypes(parentIndexes(1),1:crossPoint),genotypes(parentIndexes(2),crossPoint+1:end)]];
    offspringGenes=[offspringGenes;[genotypes(parentIndexes(2),1:crossPoint),genotypes(parentIndexes(1),crossPoint+1:end)]];
end
    
%% TODO
mutation_Prob=1/GeneLen;
for i=1:pop_size
    isMutation=rand(1,GeneLen)<mutation_Prob;
    offspringGenes(i,isMutation)=dec2bin('1'-offspringGenes(i,isMutation))';
end
genotypes=offspringGenes;
Populations=bin2dec(genotypes);
fitness=objective(Populations);
[Sort_array,index]=sort(fitness,'descend');
fitness_pop=[fitness_pop,Sort_array(1)];
for i=1:pop_size
    if fitness(i)>bestSoFarFit
        bestSoFarFit=fitness(i);
        bestSoFarSolution=Populations(i,:);
    end
end
nbGen=nbGen+1;
nbEval=nbEval+pop_size;
fitness_gen=[fitness_gen,bestSoFarFit];
solution_gen=[solution_gen,bestSoFarSolution];
end
% Mutation
%% TODO
figure,plot(1:nbGen,fitness_gen,'b');
title('fitness\_gen');

figure,plot(1:nbGen,solution_gen,'b');
title('solution\_gen');

figure,plot(1:nbGen,fitness_pop,'b');
title('fitness\_pop');





