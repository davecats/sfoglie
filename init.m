%  
% This script prepares the MATLAB environment 
%

% Clean environment
clc
close all
clear all

%  Add paths
addpath('.//')
addpath('./BoundaryElementMethod/')
addpath('./BoundaryElementMethod/geometry/')
addpath('./BoundaryElementMethod/wake/')


addpath('./BoundaryLayerSolution/')
addpath('./BoundaryLayerSolution/Amplification/')
addpath('./BoundaryLayerSolution/Transition')

addpath('./utilities/')
addpath('./utilities/Plotting/')
addpath('./utilities/SolutionProcess/')
addpath('./models/')

addpath('./optimization/')

%addpath('./pgfplots/')


set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');