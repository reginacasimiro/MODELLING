function WDN_simulation(WDNname,pipes,nodes,coords,assets,PDA,DT,ref)

%   This function sets up the steady-state simulation of a water distribution network (WDN) using "Simulation_Model_WDN.m".
%   Simulation_Model_WDN.m performs Demand-Driven Analysis (DDA) or Pressure-Driven Analysis (PDA) 
%   allowing the prediction of the background leakages by means of Germanopoulos’ model (1985). 
%   Simulation_Model_WDN.m is extracted from WDNetXL (www.hydroinformatics.it) and it is intended for research purpose 
%   in order to explain the modelling strategy reported in Giustolisi et al. (2008a; b). 
%   In fact, it is limited to background leakage and/or hydraulic pressure deficient conditions 
%   modelling and does not implement devices as for example flow-pressure control valves;  check valves, etc.. 
%
%   Suggested References
%   Germanopoulos, G. 1985 A technical note on the inclusion of pressure dependent demand and leakage terms in water supply network models. 
%   Civil Eng. Syst., 2(3), 171–179.
%   Giustolisi, O. & Walski, T.M. 2012 A Demand Components in Water Distribution Network Analysis. 
%   J. Water Resour. Plan. Manage., 138(4), 356 -367.
%   Giustolisi, O., Savic, D.A. & Kapelan, Z. 2008a Pressure-driven demand and leakage simulation for water distribution networks. 
%   J. Hydr. Eng., 134(5), 626–635.
%   Giustolisi, O., Kapelan, Z. & Savic, D.A. 2008b An algorithm for automatic detection of topological changes in water distribution networks. 
%   J. Hydr. Eng., 134(4), 435–446.
% 
%   Other References of further implementations in WDNetXL 
%   Giustolisi, O., Berardi, L., Laucelli, D. 2014 Modeling local water storages delivering customer-demands in WDN models, 
%   J. Hydr. Eng.. J. Hydr. Eng., 140(1), 1–16.
%   Giustolisi, O., Berardi, L., & Laucelli, D. 2012 “Accounting for directional devices in WDN modeling. 
%   J. Hydr. Eng., 138(10), 858-869.
%   Giustolisi, O., Berardi, L., Laucelli, D. & Savic, D.A. 2012 A Computationally efficient modeling method for large size water network analysis.
%   J. Hydr. Eng., 138(4), 313-326.
%   Giustolisi, O., Berardi, L. & Laucelli, D. 2012 Generalizing WDN simulation models to variable tank levels. 
%   J. Hydroinf., 12(3) 562–573.
%   Giustolisi, O. & Laucelli, D. 2011 Water distribution networks pressure-driven analysis using EGGA. 
%   J. Water Resour. Plan. Manage., 137(6), 117-127.
%   Giustolisi, O. 2010 Considering actual pipe connections in WDN analysis. 
%   J. Hydr. Eng. 136(11), 889–900.
%
%   Input Variables
%   WDNname     = name of the network
%   pipes       = matrix with pipes data
%   nodes       = matrix with nodes data
%   coords      = XYZ coordinates
%   assets      = matrix of asset data
%   PDA         = flag for pressure-driven analysis
%   DT          = time interval [min] of the snap shot of the steady-state (for variable tank levels)
%   ref         = mumber of step for rephinement when Tanks exist
%
%   Output variables
%   NO OUTPUT
%
%   Programmed by   Orazio Giustolisi, dICAR, Technical University of Bari
%   LastEditDate    January, 2014
%   e-mail          o.giustolisi@poliba.it
%   web page        www.hydroinformatics.it   

close all
% directory of files
folder_name      = pwd;

% TOPOLOGICAL MATRIX and OTHER
num_pipes       = size(pipes,1);                                            % number of pipes
num_nodes       = size(nodes,1);                                            % number of nodes
not_reservoirs  = nodes(:,2)==0;                                            % not reservoirs 
reservoirs      = not(not_reservoirs);                                      % reservoirs
tank_nodes      = find(reservoirs);                                         % position of reservoirs
int_nodes       = find(not_reservoirs);                                     % position of internals
A_gen           = spalloc(num_pipes,num_nodes,2*num_pipes);                 % matrix of nodes connections
A_gen((1:num_pipes)'+num_pipes.*(pipes(:,2)-1)) = -1;
A_gen((1:num_pipes)'+num_pipes.*(pipes(:,3)-1)) = +1;
A12             = A_gen(:,int_nodes);                                       % matrix of internal nodes
A10             = A_gen(:,tank_nodes);                                      % matrix of reservoir nodes
H0s             = sum(nodes(tank_nodes,1:2),2);                             % tank levels of reservoirs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTTING
set(figure(1000),'Position',[1 1 1280 700],'color',[1 1 1])
axis off
hold on
Xpoint(1,:) = coords(pipes(:,2),1);
Xpoint(2,:) = coords(pipes(:,3),1);
Ypoint(1,:) = coords(pipes(:,2),2);
Ypoint(2,:) = coords(pipes(:,3),2);
Zpoint(1,:) = coords(pipes(:,2),3);
Zpoint(2,:) = coords(pipes(:,3),3);
meanXpoint  = mean(Xpoint);
meanYpoint  = mean(Ypoint);
meanZpoint  = mean(Zpoint);
plot3(Xpoint,Ypoint,Zpoint,'k','erasemode','normal')
Xpoint = coords(reservoirs,1);
Ypoint = coords(reservoirs,2);
Zpoint = coords(reservoirs,3);
str_a = 'Ho_ ';
str_colour = 'b';
text(Xpoint,Ypoint,Zpoint,str_a,'Color',str_colour,...
    'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',10,'FontWeight','bold','FontName','Cambria')
plot3(Xpoint,Ypoint,Zpoint,'sk','MarkerFaceColor','b','MarkerSize',9)
font_size   = 8;
I_pumps = not(pipes(:,9)==0) & not(pipes(:,10)==0);
text(meanXpoint(I_pumps),meanYpoint(I_pumps),meanZpoint(I_pumps),'PM','Color','b',...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',font_size,'FontWeight','bold','FontName','Cambria',...
    'EdgeColor','k')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MODEL INPUT
Lk              = pipes(:,5);                                               % lengths
Dindex          = pipes(:,6);                                               % indices of diameters
Rk              = assets(Dindex,2).*Lk;                                     % hydraulic resistances
nk              = assets(Dindex,3);                                         % exponent of the head loss function
Kmlk            = pipes(:,7);                                               % coeffcient of minor losses
wk              = pipes(:,11);                                              % speed factor of the pumps
ck              = pipes(:,10);                                              % exponents of pump curve
Hok             = (wk.^2).*pipes(:,8);                                      % static head of pumps
rk              = (wk.^(2-ck)).*pipes(:,9);                                 % coefficents of pump curve
alfak           = pipes(:,14);                                              % exponents of the leakage model
betak           = pipes(:,15);                                              % coeffcients of the leakage model
Zs              = nodes(:,1);                                               % elevations
Pmins           = nodes(not_reservoirs,10);                                 % minimum pressure for any service
Pservs          = nodes(not_reservoirs,11);                                 % minimum pressure for correct service
ds_hum          = nodes(not_reservoirs,9);                                  % customer demand
fast            = 1;                                                        % condition for model accuracy - PI<fast*1e-7
maxiter         = 100;                                                      % maximum number of iterations
lstep           = [2.5 10 1];                                               % step sizes for overrelaxation if required and initial lambda

% Tanks
DT              = DT*60;                                                    % DT is in minutes
Hmax            = [];
Hmin            = [];
A_DT            = [];
d_ext           = [];
tank_var_nodes  = nodes(:,5)>0;
if any(tank_var_nodes) && DT>0
    Hmin        = nodes(tank_var_nodes,3)+nodes(tank_var_nodes,1);  % minimum head level 
    Hmax        = nodes(tank_var_nodes,4)+nodes(tank_var_nodes,1);  % maximum head level 
    A_DT        = nodes(tank_var_nodes,5)/DT;                               % tank cross-sectional area divided by time interval
    d_ext       = nodes(tank_var_nodes,12);                                 % fixed demand generally from an external pipe
else
    ref         = 1;
end
tank_var_nodes  = nodes(reservoirs,5)>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1 = msgbox('please wait while computing');
% SIMULATION OF A SNAP-SHOT
% WDN simulation model without VSP or PCV
[dnodes_hum,dnodes_leak,dpipes_leak,Q_pipes,H_nodes,err,t_solver,Ltanks] = Simulation_model_WDN(A12,A10,PDA,ds_hum,H0s,Zs, ...
    tank_nodes,int_nodes,Rk,Kmlk,Hok,rk,ck,nk,Pservs,Pmins,alfak,betak.*Lk,fast,ref*maxiter,lstep,ref,Hmax,Hmin,A_DT,d_ext,tank_var_nodes);

disp(['Error model = ' num2str(err(1:end-1))])
disp(['Global Mass Balance = ' num2str(err(end))])

% This is to store WDN state variables and write mat
P_nodes                             = H_nodes(int_nodes) - Zs(int_nodes);
actual_nodal_demand_hum             = dnodes_hum;
pipe_leakages                       = dpipes_leak;
nodal_pressure(reservoirs,1)        = H0s;
nodal_pressure(not_reservoirs,1)    = P_nodes;
pipe_flows                          = Q_pipes;
% This is to write results in a matlab file
filename = [folder_name '\Results_' WDNname '.mat'];
save(filename,'actual_nodal_demand_hum','nodal_pressure','pipe_flows','pipe_leakages')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTTING AND SAVE DATA
% This section creates a visualzation of the network with status at the end of DT
Pmin            = -10;
hold on
set(gcf,'Name',WDNname,'NumberTitle','off')
% this is to avoid change in plot limits when reservoirs are perifirical
Xpoint          = coords(reservoirs,1);
Ypoint          = coords(reservoirs,2);
Zpoint          = coords(reservoirs,3);
plot3(Xpoint,Ypoint,Zpoint,'b.')
x               = 10*round(log10(size(nodes,1)));
size_point      = max(80-x,10);
scatter3(coords(not_reservoirs,1),coords(not_reservoirs,2),coords(not_reservoirs,3),size_point,P_nodes,'filled');

axis off
load('MyColormaps','mycmap')
set(gcf,'Colormap',mycmap)
max_val = ceil(max(P_nodes));
min_val = fix(min(P_nodes));
if not(isempty(Pmin)) && max_val>Pmin, min_val = Pmin; end
caxis([min_val max_val])
colorbar('West','Xtick',0,'XTickLabel','Pressure [m]')
title(WDNname,'color','b','FontSize',18,'FontWeight','bold','FontName','Cambria')
hold off
delete(h1)

%   This is a function for plotting results of the simulations. It is useful to compact the code.
line_size = 1;
set(figure(50),'Position',[1 1 1280 700],'color',[1 1 1])
set(gcf,'Name','RESULTS','NumberTitle','off')
subplot(2,1,1)
x = (1:sum(not_reservoirs))';
plot(x,P_nodes,'bo--',x,H_nodes(not_reservoirs),'rs--','LineWidth',line_size)
if length(x) < 50, grid MINOR, else grid ON, end
title('NODAL PRESSURE-HEAD','FontSize',12,'FontWeight','bold')
xlabel('Node # (internal)','FontSize',14,'FontWeight','bold','FontName','Times New Roman')
ylabel('Nodal Pressure-Head [m]','FontSize',14,'FontWeight','bold','FontName','Times New Roman')
legend('Pressures','Heads')

subplot(2,1,2)
plot(x,dnodes_hum*1000,'bo-',x,dnodes_leak*1000,'rs-','LineWidth',line_size)

if length(x) < 50, grid MINOR, else grid ON, end
title('ACTUAL NODAL DEMAND COMPONENTS','FontSize',12,'FontWeight','bold')
xlabel('Node # (internal)','FontSize',14,'FontWeight','bold','FontName','Times New Roman')

ylabel('Actual Nodal Demand Components [l/s]','FontSize',14,'FontWeight','bold', 'FontName','Times New Roman')
legend('Customer-based demands','Leakages')

if ref>1
    interval        = linspace(DT/(60*ref),DT/60,ref);
    str_x_axis      = 'Step [min]';
    ref             = length(interval);
    line_size       = 2;
    title_size      = 14;
    axis_size       = 12;
    name_font       = 'Times New Roman';
    label_L         = []; %cell(size(L_tanks,1),1);
    for k = 1:size(Ltanks,2)
        label_L{k} = ['T' num2str(k)];
    end
    set(figure(51),'Position',[1 1 1280 700],'color',[1 1 1])
    hold on
    x = plot([0 interval],Ltanks,'LineWidth',line_size);
    set(gcf,'Name',[WDNname ' WDN - Tank Levels'],'NumberTitle','off')
    
    if ref < 50, grid MINOR, else grid ON, end
    title('TANK LEVELS','FontSize',title_size,'FontWeight','bold','FontName',name_font)
    xlabel(str_x_axis,'FontSize',axis_size,'FontWeight','bold','FontName',name_font)
    ylabel('Level [m]','FontSize',axis_size,'FontWeight','bold', 'FontName',name_font)
    L = legend(x, label_L);
    set(L,'FontSize',axis_size-1,'FontName',name_font);
    hold off
end



h = msgbox('     End of the Analysis','END','help');
uiwait(h)
end

