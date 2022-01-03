close all
clear
clc
tic

%% This code was written by Daniel Plaugher and David Murrugarra 
% Boolean model of pacreatic cancer from papers:

% Modeling the Pancreatic Cancer Microenvironment in Search of Control Targets
% -Daniel Plaugher and David Murrugarra 

% Formal Modeling and Analysis of Pancreatic Cancer Microenvironment
% - Qinsi Wang et al

%% Add pathways to "MatlabToolBoxes" folders
% addpath('/Users ... /SDDS');
% addpath('/Users ... /BNPBN/');
% addpath('/Users/... /BNPBN/BNPBN/BNPBN');


% Number of inputs
n = 69; % whole system
p = 2; % number of states
c = 0.9*ones(2,n); % normal propensities for SDDS

%%% entire system
varF=load('comb_pcc_varF.txt'); % variable communications per node
nv=load('comb_pcc_nv.txt'); % number of variable communications per node
F = load('comb_pcc_tt.txt'); % truth table function 


%% Inductions and Control
% --- Mutations ~ TP53-64(off)	CDKN2A-56(off)	SMAD4-44(off)	KRAS-43(on)
% F1 = TruthTable_del_n_temp(F,nv,varF,p, 43,1); % regulates KRASp
% F2 = TruthTable_del_n_temp(F,nv,varF,p, 64,0); % regulates TP53p
% F3 = TruthTable_del_n_temp(F,nv,varF,p, 56,0); % regulates CyclinDp
% F4 = TruthTable_del_n_temp(F,nv,varF,p, 44,0); % regulates SMADp

% --- Node control ~ PI3K(42), PIP3(46), RAF(47), BAX(61),
% F5 = TruthTable_del_n_temp(F4,nv,varF,p, 46,0); % regulates 
% F6 = TruthTable_del_n_temp(F5,nv,varF,p, 47,0); % regulates 
% F7 = TruthTable_del_n_temp(F6,nv,varF,p, 29,0); % regulates ERKS
% F8 = TruthTable_del_n_temp(F7,nv,varF,p, 18,0); % regulates PI3K

% --- Edge control
% F5 = TruthTable_del_a_temp(F4,nv,varF,p,46,51,0); %edge control PIP3p-AKTp
% F6 = TruthTable_del_a_temp(F5,nv,varF,p,47,49,0); %edge control RAFp-MEKp
% F6 = TruthTable_del_a_temp(F5,nv,varF,p,61,66,1); %edge control BAXp-CASPp
% F7 = TruthTable_del_a_temp(F6,nv,varF,p,29,30,0); %edge control ERKs-AP1s
% F8 = TruthTable_del_a_temp(F7,nv,varF,p,18,25,0); %edge control PI3Ks-PIP3s


%% Simulation
nins = 1000; % number of initializations
nsteps=300; % number of steps for SDDS
g=0.01; % noise

[Y,My]=SDDS_simNoise(g,F,varF,nv,p,c,n, nsteps,nins); 
% [Y,My]=SDDS_sim(F,varF,nv,p,c,n, nsteps,nins); % simulation

Ylast=Y(:,end); %  end trajectories for each node -- can interpret as probability of node expression in long term
Phen = [Ylast(33); Ylast(34); Ylast(35); Ylast(36); Ylast(67); Ylast(68); Ylast(69)] %end phenotype readings
%      ('Apops',   'Prols',   'Migs',    'Acts',    'Autp',    'Apop',    'Prop')



%% Graphing 
X = 0:1:nsteps; % time steps
% 
%whole
figure('Name', 'PC Simulation')
plot(X,Y(33,:),'c',X,Y(34,:),'m',X,Y(35,:),'k',X,Y(36,:),'y',X,Y(67,:),'g',X,Y(68,:),'b', X,Y(69,:),'r','LineWidth',1.5,'MarkerSize',10)
legend('Apops','Prols', 'Migs','Acts', 'Autp','Apop','Prop')
xlabel('Time Steps')
ylabel('Average Frequencies')


%% k means  - used to confirm statistical analysis
clear
clc
close



% X= [0    4255.925;
% 1    2223.360;
% 2    1513.170;
% 3    1505.100]; % MDM2 medians
% 
% X= [0    6490.965;
% 1    3652.850;
% 2    3105.320;
% 3    3086.555]; %Bax medians

X=[0     27.25;
1    309.07;
2     37.22;
3     24.20]; % PIK3CD medians



colmax=max(X);
Xnew= X./colmax;% scaled X


[idx,C]= kmeans(X,2);
C_up = C.*colmax;


figure;
hold on
gscatter(X(:,1),X(:,2),idx,'bgm')
plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
title 'Cluster Assignments and Centroids'
hold off


toc


