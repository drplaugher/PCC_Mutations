%% This code was written by Daniel Plaugher 
% -- finds strongly connected components (modules) of networks

clear
clc
close all

%% Add pathways to "MatlabToolBoxes" folders
% addpath('/Users ... /SDDS');
% addpath('/Users ... /BNPBN/');
% addpath('/Users/... /BNPBN/BNPBN/BNPBN');


%% PCC model set up for SDDS

varF=load('comb_pcc_varF.txt'); % variable communications per node
nv=load('comb_pcc_nv.txt'); % number of variable communications per node
F = load('comb_pcc_tt.txt'); % truth table function 


% -- induce mutations
% f(43)=1; % KRAS
% f(44)=0; % SMAD
% f(56)=0; % CycD
% f(64)=0; % TP53


%% Setting up tails and heads   

[row,col]=find(varF~=-1);

t=NaN(1,length(row)); % tail
h=NaN(1,length(row)); % head
for i=1:length(row)
    t(1,i)= varF(row(i),col(i)); %builds vector of tails
    h(1,i) = col(i); % builds vector of heads
end

G = digraph(t,h,[],n);
plot(G,'Layout','layered')
%% Finding modules
[bin,binsize] = conncomp(G); % finds connected components of graph 

mods = find(binsize~=1); % nontrivial modules found

mod_nodes=-1*ones(max(binsize),length(mods)); % nodes in each module in columns
for i=1:length(mods)
    mod_nodes(1:binsize(mods(i)),i)=find(bin==mods(i))';   
end

modTable= [mods; NaN(1,length(mods)); mod_nodes]; % table with modules and their nodes

%% Graphical Analysis

% to find module rankings
p = plot(G);
p.MarkerSize = 7;
p.NodeCData = bin;
colormap(hsv(length(binsize)))

C = condensation(G);
figure('Name','Condensation Graph')
p2 = plot(C);
p2.MarkerSize = 7;
p2.NodeCData = 1:length(binsize);
colormap(hsv(length(binsize)))

% topological sorting
TopSort = toposort(C); % sorts according to topological importance/rank
ModRanks = [mods; find(ismember(TopSort,mods)) ] %finds location of nontrivial mods in TopSort
RowNames = ["modules"; "Ranks"];
TopSortModRanks = [RowNames ModRanks];



