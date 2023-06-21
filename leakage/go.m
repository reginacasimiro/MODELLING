clear all
clc

load C-TOWN

PDA                     = 1;            % Pressure-Driven and Leakages
n_hours                 = 24;           % number of hours 
DT                      = n_hours*60;   % time interval of the snapshot
ref                     = n_hours*20;   % number of refinements for each hour if tanks

WDN_simulation(WDNname,pipes,nodes,coords,assets,PDA,DT,ref);