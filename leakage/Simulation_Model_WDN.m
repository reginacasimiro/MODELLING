function [qnode_act_hum,qnode_leak,qpipe_leak,Q,H,err,t_solver,Ltanks] = ...
    Simulation_Model_WDN(A12,A10,sim_PD,demand,H0,Hg,tank_nodes,int_nodes,RL,ML,Hp,Rp,mp,n,Pserv,Pmin,alfa,betaL,fast,maxiter,lstep,...
    ref,Hmax,Hmin,A_DT,d_ext,tank_var_nodes)

%   This function executes the steady-state simulation of a water distribution network (WDN).
%   The function performs Demand-Driven Analysis (DDA) or Pressure-Driven Analysis (PDA) 
%   allowing the prediction of the background leakages by means of Germanopoulos’ model (1985). 
%   The function is extracted from WDNetXL (www.hydroinformatics.it) and it is intended for research purpose 
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
%
%   Input for GGA (PDA or DDA)
%   A12         = incidence sub-matrix (unknown heads)
%   A10         = incidence sub-matrix (known heads)
%   sym_PD      = Pressure-Driven simulation (1) reordering (2) solver type (3)
%   demand      = matrix of demands in normal working conditions
%   H0          = known heads of reservoir
%   Hg          = nodal elevation (reservoirs comprised)
%   tank_nodes  = indices of source nodes (known heads)
%   int_nodes   = indices of internal nodes (unknown heads)
%   RL          = pipe hydraulic resistance in fully turbulent flow
%   ML          = minor head losses DH = ML*Q^2
%   Hp,Rp,mp    = pump curve Hop-Rp*Q^mp
%   n           = exponent of head loss formulae
%   Pserv,Pmin  = minimum pressure for a correct service or anty service for Wagner's model
%   alfa,betaL  = parameters of the Germanopoulos model
%
%   Input for accuracy and search
%   fast        = condition PI<fast*1e-7
%   maxiter     = maximum number of iterations
%   lstep       = step sizes for overrelaxation if required and initial lambda
%
%   Input for Tanks
%   ref         = number of refinements
%   Hmax,Hmin   = tank_nodes => minimum and maximum tank levels above elevation
%   A_DT        = tanks cross-section area divided by time interval
%   d_ext       = demand from external supply never pressure-dependent
%   tank_var_nodes = logical vector identifying tanks among reservoirs
%
%   Output variables
%   qnode_act_hum,qnode_leak,qpipe_leak = vectors of leaks (pipe and nodal levels) and actual demand
%   [Q H]       = state variables
%   err         = statistics
%
%   Programmed by   Orazio Giustolisi, dICAR, Technical University of Bari
%   LastEditDate    January, 2014
%   e-mail          o.giustolisi@poliba.it
%   web page        www.hydroinformatics.it

warning off MATLAB:nearlySingularMatrix
warning off MATLAB:nearlySingularMatrixUMFPACK
warning off MATLAB:singularMatrix

max_lambda              = 1;                                            % maximum value of over-relaxation parameter
lambda_used             = 0;                                            % how many times is used over-relaxation
lambda                  = lstep(3);                                    	% initial step for H and Q
regolarize_pij          = 1e-3;                                         % this is for regularization
A12                     = sparse(A12);                                  % to be sure that is sparse
A10                     = sparse(A10);                                  % to be sure that is sparse
t_solver                = zeros(maxiter,1);                             % initialize stats for solver

if all((Pserv-Pmin)) && (all(alfa==0) || all(betaL==0))
    sim_PD(1)           = 0;
end

if sim_PD(1)
    init_Q(1)           = .1;
else
    init_Q(1)           = 1;
end

% Assigning memory space to variables to speed up
Njunct                  = length(demand);
Nlink                   = length(RL);
Ntank                   = length(H0);
Ntot                    = Njunct+Ntank;
% Useful for PD analysis
A12G(:,[int_nodes; tank_nodes]) = [A12 A10];

% modification to have n for each pipe
if length(n) == 1,
    n                   = n*ones(Nlink,1);
end

% Useful later on in PDA
A21G                    = A12G';
A12Gabs                 = abs(A12G);
OG                      = abs(A21G);
A21                     = A12';

% Ordersing to reduce fill-in in the following factorization
%Approximate minimum degree permutation
p                       = amd(A21*A12);

% useful to initialize anyway
qpipe_leak              = zeros(Nlink,1);
qpipe_leak_new          = zeros(Nlink,1);
D22                     = zeros(Njunct,1);
qnode_leak              = zeros(Njunct,1);
qnode_leak_new          = zeros(Njunct,1);
qnode_act_hum           = zeros(Njunct,1);
qnode_act_hum_new       = zeros(Njunct,1);
H_new                   = zeros(Njunct,1);
d22                     = zeros(Njunct,1);
l22                     = zeros(Ntot,1);

% some initializzazions
cont_lsol               = 0;                                                     	% initialize counter of linear solutions
iteration               = 1;                                                     	% initialize iteration
minPI                   = 1e-7;                                                   	% for convergence criterion
H                       = Pserv/2+Hg(int_nodes);                                    % for PD
I_mp                    = logical(mp);                                              % this is a logical index useful for CHV
pif                     = A10*H0-sign(Rp).*Hp;                                      % for tanks and pumps in GGA

% construction of element for leaks
if sim_PD(1)
    Pserv_res               = Pserv-Pmin(:,1);                              % minimum residula pressure for a correct service (i.e. supply d_hum+d_vol)
    sqrtPserv_res           = sqrt(Pserv_res);                              % also useful later on
    C_hum                   = demand(:,1)./sqrtPserv_res;                   % coefficient of human-based demand curve (Giustolisi and Walski, 2012)
    C_hum(isnan(C_hum) | isinf(C_hum)) = 0;

    Pnodes(int_nodes,1)     = H-Hg(int_nodes);                              % nodal pressure in internal
    Pnodes(tank_nodes,1)    = H0-Hg(tank_nodes);                            % nodal pressure in tank
    Pmean                   = 0.5*(A12Gabs*Pnodes);                         % mean pressure in pipes 

    if any(C_hum)
        I                   = C_hum>0;                                      % pressure driven only inflow
        % model for indoor demands (i.e. human and volume based demands)
        [xx,zz]             = Indoor_demands(demand(I,1),C_hum(I),Pnodes(int_nodes(I)),Pserv_res(I),Pmin(I));
        qnode_act_hum(I)    = xx;
        d22(I)              = zz;
    end

    if any(betaL) && any(alfa)
        [qpipe_leak,l22]    = Background_leakages(Pmean,alfa,betaL);        % Germanopoulos model
        qnode_leak          = 0.5*(OG*qpipe_leak);                          % transform in nodal loads
        qnode_leak          = qnode_leak(int_nodes);                        % only internal useful
        l22                 = 0.5*(OG*l22);                                 % transform in nodal loads
    end
    D22                     = -sum(d22+l22(int_nodes),2);                   % derivative
end
qnodes                      = demand*(sim_PD(1)==0)+qnode_leak+qnode_act_hum;   % nodal demand

t = tic;
% initialize Q
pij                         = RL+ML;
pij                         = 1./pij;
D11_inv                     = spdiags(pij,0,Nlink,Nlink);
A                           = A21*D11_inv*A21';
DD                          = spdiags(D22,0,Njunct,Njunct);
A                           = A - DD;
pijf                        = pij.*pif;
% construction of F
F                           = A21*pijf+init_Q(1)*qnodes;
HH                          = -A(p,p)\F(p);
HH(p,1)                     = HH;
% Q solution
Q                           = -(pij.*(A21'*HH)+pijf);
iteration                   = iteration + 1;

% starting PIs
AQ                          = Q.*(ML.*abs(Q)+RL.*abs(Q).^(n-1));
AQ(I_mp)                    = AQ(I_mp) + Q(I_mp).*((Rp(I_mp).*Q(I_mp)>0).*(abs(Rp(I_mp)).*abs(Q(I_mp)).^(mp(I_mp)-1)));

PI_Continuity               = 1e30;
PI_Energy                   = 1e30;
PI                          = 1e30;

compute_HQ_new              = 1;

% this is useful to correctly update with overrelaxing in PDA
qnodes_new                  = qnodes;
D22_new                     = D22;

% for tanks
Ltanks                      = zeros(sum(tank_var_nodes),ref+1);
if not(isempty(Ltanks))
    Ltanks(:,1)             = H0(tank_var_nodes) - Hg(tank_nodes(tank_var_nodes));
end
for j=1:ref
    while (iteration<=maxiter && lambda>1e-7) && (PI>fast*minPI)
        if compute_HQ_new || lambda<=1e-7
            cont_lsol   = cont_lsol + 1;
            t1                      = tic;
            pij                     = 2*ML.*abs(Q)+n.*RL.*abs(Q).^(n-1);
            pij(I_mp)               = pij(I_mp) + (Rp(I_mp).*Q(I_mp)>0).*(mp(I_mp).*abs(Rp(I_mp)).*abs(Q(I_mp)).^(mp(I_mp)-1));
            D22                     = D22_new;
            pij                     = 1./(pij+(pij<=regolarize_pij)*regolarize_pij);
            D11_inv                 = spdiags(pij,0,Nlink,Nlink);
            DD                      = spdiags(D22,0,Njunct,Njunct);
            pijf                    = pij.*pif;
            B11                     = pijf-(Q-AQ.*pij);
            A                       = A21*D11_inv*A12;
            % construction of F
            F                       = A21*B11+qnodes_new+D22.*H;
            % Head solutions
            A                       = A-DD;
            H_new(p,1)              = -A(p,p)\F(p);
            % Q solution
            Q_new                   = -B11-pij.*(A12*H_new);
            t_solver(cont_lsol)     = toc(t1);
            H_new_store             = H_new;
            Q_new_store             = Q_new;
        else
            H_new                       = H_new_store;
            Q_new                       = Q_new_store;
        end
        
        % update with step size lambda
        Q_new                           = lambda*(Q_new - Q) + Q;
        H_new                           = lambda*(H_new - H) + H;
        
        if sim_PD(1)
            Pnodes(int_nodes)           = H_new-Hg(int_nodes);                      % nodal pressure in internal
            Pmean                       = 0.5*(A12Gabs*Pnodes);                     % mean pressure in pipes
            qnode_act_hum_new(:)        = 0;
            
            if any(C_hum)
                I                       = C_hum>0;
                % model for indoor demands (i.e. human and volume based demands)
                [xx,zz]                 = Indoor_demands(demand(I,1),C_hum(I),Pnodes(int_nodes(I)),Pserv_res(I),Pmin(I));
                qnode_act_hum_new(I)    = xx;
                d22(I)                  = zz;
            end
            
            if any(betaL) && any(alfa)
                [qpipe_leak_new,l22]    = Background_leakages(Pmean,alfa,betaL);  % Germanopoulos model
                qnode_leak_new          = 0.5*(OG*qpipe_leak_new);                          % transform in nodal loads
                qnode_leak_new          = qnode_leak_new(int_nodes);                        % only internal useful
                l22                     = 0.5*(OG*l22);                                     % transform in nodal loads
            end
            D22_new              	= -sum(d22+l22(int_nodes),2);
        end
        qnodes_new = demand*(sim_PD(1)==0)+qnode_leak_new+qnode_act_hum_new;    % nodal demand
        
        % compute some PIs
        AQ                      = Q_new.*(ML.*abs(Q_new)+RL.*abs(Q_new).^(n-1));
        AQ(I_mp)                = AQ(I_mp) + Q_new(I_mp).*((Rp(I_mp).*Q_new(I_mp)>0).*(abs(Rp(I_mp)).*abs(Q_new(I_mp)).^(mp(I_mp)-1)));
        Energy_new              = AQ+A12*H_new+pif;                     	% energy balance equation in pipes
        Continuity_new          = A21*Q_new-qnodes_new;                  	% mass balance equation at nodes
        PI_Energy_new           = max(abs(Energy_new));                    	% PI for energy
        PI_Continuity_new       = max(abs(Continuity_new));              	% PI for mass
        E_new                   = [Energy_new; Continuity_new];            	% mathematical system equations
        PI_new                  = (E_new'*E_new)/length(E_new);            	% PI of error surface for iteration meaning SSE
        
        % Update only if criterion has decreased
        if (PI_new-PI)<=0
            
            if sim_PD(1)
                qpipe_leak      = qpipe_leak_new;
                qnode_leak      = qnode_leak_new;
                qnode_act_hum   = qnode_act_hum_new;
                qnodes          = qnodes_new;
            end
            H                   = H_new;                                    % assign new values
            Q                   = Q_new;                                  	% assign new values
            Energy              = Energy_new;                               % assign new values
            Continuity          = Continuity_new;                           % assign new values
            PI_Energy           = PI_Energy_new;                            % assign new values
            PI_Continuity       = PI_Continuity_new;                        % assign new values
            PI                  = PI_new;                                   % assign new values
            lambda              = lambda*lstep(1);
            if lambda>=max_lambda
                lambda          = max_lambda;
            end
            compute_HQ_new      = 1;                                        % compute parameters again
        else
            lambda_used         = lambda_used + 1;                          % increase # of time is used lambda
            lambda              = lambda/lstep(2);                          % decrease step size
            compute_HQ_new      = 0;                                        % do not compute parameters
        end
        iteration               = iteration + 1;
    end
    
    % Mass balance at the tanks is here decoupled with respect to energy balance
    % (i.e. the increase of tank level is not considered during the steady-state simulation run)
    % In orderr to better simulate tanks it is necessary to couple mass balance equation as in WDNetXL system, see
    % Giustolisi, O., Berardi, L. & Laucelli, D. 2012 Generalizing WDN simulation models to variable tank levels. 
    % J. Hydroinf., 12(3) 562–573.
    if any(tank_var_nodes)
        disp(['Error model before tanks mass balance= ' num2str([PI PI_Energy PI_Continuity mean(abs(Energy)) mean(abs(Continuity)) iteration-1])])
        
        % new reservoir level based on mass balance in the reservoir nodes
        qnode_tank_leak             = 0.5*(OG(tank_nodes(tank_var_nodes),:)*qpipe_leak);   % nodal leakages of tanks
        H0(tank_var_nodes)          = (A10(:,tank_var_nodes)'*Q-qnode_tank_leak+d_ext)./(A_DT*ref) + H0(tank_var_nodes);
       
        % this is to reset to maximum
        H0(tank_var_nodes)          = (H0(tank_var_nodes)>Hmax).*Hmax + (H0(tank_var_nodes)<=Hmax).*H0(tank_var_nodes);
        
        % this is to reset to minimum
        H0(tank_var_nodes)          = (H0(tank_var_nodes)<Hmin).*Hmin + (H0(tank_var_nodes)>=Hmin).*H0(tank_var_nodes);
        
        % this is to store for plotting
        Ltanks(:,j+1)               = H0(tank_var_nodes) - Hg(tank_nodes(tank_var_nodes));
        pif                         = A10*H0-sign(Rp).*Hp;
        if j<ref, PI = inf; lambda = 1; end
        
        % Actual mass and energy balance after corrections
        Energy_1                = AQ+A12*H+pif;                             % energy balance equation in pipes
        Continuity_1            = A21*Q_new-qnodes;                         % mass balance equation at nodes
        PI_Energy_1             = max(abs(Energy_1));                    	% PI for energy
        PI_Continuity_1         = max(abs(Continuity_1));                   % PI for mass
        E_1                     = [Energy_1; Continuity_1];                 % mathematical system equations
        PI_1                    = (E_1'*E_1)/length(E_1);                   % PI of error surface for iteration meaning SSE
        disp(['Error model after tanks mass balance = ' num2str([PI_1 PI_Energy_1 PI_Continuity_1 mean(abs(Energy_1)) mean(abs(Continuity_1))])])
    end
end
time_req                                    = toc(t);
t_solver                                    = t_solver(t_solver>0);

% this is to consider global mass balance
q(tank_nodes,1)                             = A10'*Q;
q(int_nodes,1)                              = qnodes;
mass_balance_reservoirs                     = abs(sum(q));

if not(sim_PD(1))
    qnode_act_hum                           = demand;
end


H(int_nodes)                                = H;
H(tank_nodes)                               = H0;
err(1)                                      = PI;
err(2)                                      = PI_Energy;
if isempty(PI_Continuity), PI_Continuity = 0; end
err(3)                                      = PI_Continuity;
err(4)                                      = mean(abs(Energy));
err(5)                                      = mean(abs(Continuity));
err(6)                                      = iteration-1;
err(7)                                      = cont_lsol;
err(8)                                      = lambda;
err(9)                                      = time_req;
err(10)                                     = lambda_used;
err(11)                                     = mass_balance_reservoirs;
end

function [qleak,l22kk] = Background_leakages(P,alfa,beta)

% Germanopoulos' model for Background leakages
%
%   Germanopoulos, G. (1985). “A technical note on the inclusion of pressure dependent demand and leakage terms in 
%   water supply network models.” Civil Engineering Systems, 2(September), 171–179.
%   Giustolisi, O., Savic, D.A. & Kapelan, Z. 2008a Pressure-driven demand and leakage simulation for water distribution networks. 
%   J. Hydr. Eng., 134(5), 626–635.
%
%   Input Variables
%   P           = actual pressure
%   alfa        = exponent of Germanopoulos' model
%   beta        = background leakage flow for unitary pressure
%
%   Output variable
%   qleak       = actaul leakage flow
%   l22kk       = derivative of Germanopoulos' curve
%
%   Programmed by   Orazio Giustolisi, dICAR, Technical University of Bari
%   LastEditDate    January, 2014
%   e-mail          o.giustolisi@poliba.it
%   web page        www.hydroinformatics.it

m               = length(P);
qleak           = zeros(m,1);                       % initialize leakage flow rate
l22kk           = zeros(m,1);                       % initialize to zero elements of the derivative
I               = (P>0 & beta>0 & alfa>0);          % index to compute functions
qleak(I)        = beta(I).*P(I).^alfa(I);           % leakage flow rate
l22kk(I)        = alfa(I).*qleak(I)./P(I);          % derivative
end

function [qact_hum,d22kk] = Indoor_demands(d_hum,C_hum,P,Pserv_res,Pmin)

% (d_hum,C_hum,d_vol,bmax,P,Pserv_res,Pmin,x,a,b,P0,smooth)
% Model for indoor demands (i.e. human and volume based demands) 
%
%   Wagner, J.M., Shamir, U., and Marks, D.H. (1988). “Water distribution reliability: simulation methods.” 
%   J. Water Resour. Plan. and Manage., 114(3), 276-294.
%   Giustolisi, O. & Walski, T.M. 2012 A Demand Components in Water Distribution Network Analysis. 
%   J. Water Resour. Plan. Manage., 138(4), 356 -367.
%
%   Input Variables
%   d_hum       = statistical required human-based demand over time in normal working conditions (i.e. P>Pserv) 
%   C_hum       = coefficient of human-based demand curve, usually C=qserv/sqrt(Pserv-res)
%   P           = model pressure
%   Perv_res    = minimum residula pressure for a correct service (i.e. supply d_hum+d_vol)
%   Pmin        = minimum pressure for supplying water
%
%   Output variable
%   qact_hum    = actual supplied human-based demand
%   d22kk       = sum od the derivative of the models
%
%   Programmed by   Orazio Giustolisi, dICAR, Technical University of Bari
%   LastEditDate    January, 2014
%   e-mail          o.giustolisi@poliba.it
%   web page        www.hydroinformatics.it

Pres                = P - Pmin;                                         % residual pressure
m                   = length(Pres);                                     % length of the vector
y                   = zeros(m,1);                                       % initialize to zero
d22kk               = zeros(m,1);                                       % initialize to zero
qact_hum            = d_hum;                                            % initialize to d_hum (i.e. to normal conditions P>Pserv)
I                   = (Pres>0) & (Pres<=Pserv_res);                     % indices of 0<Pres=<Pserv_res (i.e. nodes in pressure-deficient conditions)
y(I)                = sqrt(Pres(I));                                    % DH^0.5 of function qact = C*DH^0.5

if any(d_hum>0)
    II              = I & d_hum>0;
    qact_hum(II)    = y(II).*C_hum(II);                                 % actual human-based demand in pressure-deficient conditions
    d22kk(II)       = 0.5*C_hum(II)./y(II);                             % derivative of human-based demand model
end

I                   = Pres<=0 & d_hum>0;                                % indices of P=<0
qact_hum(I)         = 0;                                                % actual human-based demand equal to exp(...) if P<=P0
end
