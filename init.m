%  
% This script prepares the MATLAB environment 
%


% Clean environment
close all
clear all

%  Add paths
addpath('./panel/')
addpath('./Amplification/')
addpath('./geometry/')
addpath('./utilities/')
addpath('./models/')
addpath('./wake/')
addpath('./BoundaryLayer/')
set(groot, 'defaultAxesTickLabelInterpreter','LaTex'); set(groot, 'defaultLegendInterpreter','LaTex');