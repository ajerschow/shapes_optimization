% Leonard Strouk, Boris Kharkov
% Developped for Spinach 1.9
%
% In this simulation, the signal from conducitve part (metal) should be
% averaged out, while the signal from the liquid part is maximized.

% varargin is either empty or containes initial guess for the waveform

function waveform = MAIN_Minimization_Metal_and_Liquid(varargin) 

%% Global parameters for the grape optimization

POWER_LEVEL            = -1*2*pi*1e3;   % 1kHz
PULSE_DURATION         = 1.6e-3;        % 1.6 ms
N_STEPS                = 10;            % Nbr of points in the waveform
GRAPE_PLOT_MONITORING  = 'yes';         % If you want to follow how the waveform
                                        % is evolving, put 'yes' or 'y'

%% Basic Spinach Procedure

% Disable console output
sys.output = 'hush';

% Magnetic field
sys.magnet=11.7; % 500 MHz

% System specification
Nbr_Spins = 150;          % Results depend on the number of spins. For accuracy this parameter should not be less than 100
sys.isotopes = cell(1, Nbr_Spins);
sys.isotopes(:) = {'1H'};
inter.zeeman.scalar=num2cell(zeros(1,Nbr_Spins));

% basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

%% Livouillian/Hamiltonian of the system

% Get the drift Liouvillian (including relaxation, if any)
L=hamiltonian(assume(spin_system,'nmr'))+1i*relaxation(spin_system);

%% Initial state 

% Set up and normalize the initial state
a = 0; % Lx
b = 1; % Lz
c = 0; % Ly
rho=(...
    a*((state(spin_system,'L+','all')+state(spin_system,'L-','all'))/2)+...
    b*state(spin_system,'Lz','all')+...
    c*((state(spin_system,'L+','all')-state(spin_system,'L-','all'))/(2i))...
    )...
    /sqrt(a^2+b^2+c^2);

%% Skin Depth

% Take into account the skin depth

phi = linspace(0,2*pi,numel(spin_system.logs.sys.isotopes));

%% Control and State operators - Normalization factor

[Operators] =  Operator_Generator_Metal(spin_system,length(rho),phi);

%% Target State

Target_X = cell(1,Nbr_Spins);
Target_Y = cell(1,Nbr_Spins);
Target_Z = cell(1,Nbr_Spins);

for j = 1:Nbr_Spins
    Target_X{j} = Operators.Lx_State{j};
    Target_Y{j} = Operators.Ly_State{j};
    Target_Z{j} = state(spin_system,{'Lz'},{j});
end

% Set up the optimization target
target_1 = spalloc(length(rho),1,length(rho));
target_2 = spalloc(length(rho),1,length(rho));
target_3 = spalloc(length(rho),1,length(rho));

for j=1:Nbr_Spins
    target_1 = target_1 + Target_X{j};
    target_2 = target_2 + Target_Y{j};
end
target_3 = Target_X{1};

target_1 = - target_1; target_2 = - target_2; target_3 = - target_3;

%% Define control parameters

% Grape 1
ctrl_sys_1.method='grape';                                                 % Algorithm
ctrl_sys_1.control_ops={Operators.SpinX,Operators.SpinY};                  % Control operators
ctrl_sys_1.rho=rho;                                                        % Starting state
ctrl_sys_1.target=target_1;                                                % Destination state
ctrl_sys_1.pulse_duration=PULSE_DURATION;                                  % Pulse duration                                                     % Number of steps
ctrl_sys_1.nsteps=N_STEPS;                                                 % Number of steps
ctrl_sys_1.time_step=ctrl_sys_1.pulse_duration/ctrl_sys_1.nsteps;          % Time step
ctrl_sys_1.power_level = POWER_LEVEL;                                      % Set the power level                                     % Set the power level
ctrl_sys_1=control_sys(spin_system,ctrl_sys_1);                            % Make control system structure

% Grape 2
ctrl_sys_2 = ctrl_sys_1;
ctrl_sys_2.target=target_2;                                                % Destination state
ctrl_sys_2=control_sys(spin_system,ctrl_sys_2);                            % Make control system structure

% Grape 3
ctrl_sys_3 = ctrl_sys_1;
ctrl_sys_3.target=target_3;                                                % Destination state
ctrl_sys_3=control_sys(spin_system,ctrl_sys_3);                            % Make control system structure

%% Define optimisation parameters
optim.method='lbfgs';                                                      % Optimisation method
optim.max_iterations=Inf;                                                  % Maximum number of iterations
optim.extremum='minimum';                                                  % find minimum fidelity
optim.tol_linesearch_fx = 1e-4;                                            % cf. Nocedal and Wright optimization book
optim.tol_x = 1e-9;
optim.tol_gfx = 1e-9;

%% Guess for the control phase sequence
if (nargin == 1)
    guess = varargin{1};
else
    guess = zeros(numel(ctrl_sys_1.control_ops),ctrl_sys_1.nsteps);
end

%% Run the optimization
waveform = fminnewton(@cost_function,guess,optim);

% Build the error functional as a sum over power levels
    function [cost,cost_grad]=cost_function(waveform)
        
        % grape 1
        [diag_data_1,cost_1,cost_grad_1]=grape_modified_version(spin_system,ctrl_sys_1,L,waveform,'modified');            
        % grape 2
        [diag_data_2,cost_2,cost_grad_2]=grape_modified_version(spin_system,ctrl_sys_2,L,waveform,'modified');       
        % grape 3
        [diag_data_3,cost_3,cost_grad_3]=grape_modified_version(spin_system,ctrl_sys_3,L,waveform,'original');
        
        %% total cost
        
        weight_1 = 1;
        weight_2 = 1;
        weight_3 = 25;
        
        cost = (weight_1*cost_1 + weight_2*cost_2 + weight_3*cost_3)/(weight_1+weight_2+weight_3);
        cost_grad = (weight_1*cost_grad_1 + weight_2*cost_grad_2 + weight_3*cost_grad_3)/(weight_1+weight_2+weight_3);
        diag_data.trajectory = (diag_data_1.trajectory + diag_data_2.trajectory + diag_data_3.trajectory)/3;

        %% Penalize excursions outside the power envelope
        
        [pen,pen_grad]=penalty(waveform,'SNS',-ones(size(waveform)),ones(size(waveform)));
        cost=cost+pen; 
        cost_grad=cost_grad+pen_grad;
        
        %% Plot the current waveform
        
        if (strcmp(GRAPE_PLOT_MONITORING,'yes') || strcmp(GRAPE_PLOT_MONITORING,'y'))
                display_grape_results(diag_data.trajectory,ctrl_sys_1,waveform,Operators,3);
        end
       
    end

Trajectory = waveform_to_final_rho(spin_system, L, ctrl_sys_1, waveform);
[X_1,Y_1,Z_1] = display_grape_results(Trajectory,ctrl_sys_1,waveform,Operators,4);
POP_1 = [X_1;Y_1;Z_1];

%%
fprintf('\n');
fprintf('********************************************************* \n');
fprintf('************************* Final ************************* \n');
fprintf('********************************************************* \n');
fprintf('\n');

fprintf('Population of the final states \n');
fprintf('X = %f \n', POP_1(1,end));
fprintf('Y = %f \n', POP_1(2,end));
fprintf('Z = %f \n', POP_1(3,end));

fprintf('\n');

fprintf('hdot(Lx_State{1},Lx_State{1})       = %f \n', full(hdot(Operators.Lx_State{1},Operators.Lx_State{1})));
fprintf('hdot(Trajectory{end,1},Lx_State{1}) = %f \n', full(hdot(Trajectory{end,1},Operators.Lx_State{1})));

fprintf('\n');
fprintf('\n');

fprintf('FIRST SPIN \n')

fprintf('X = %f \n', full(hdot(Trajectory{end,1},Operators.Lx_State{1})));
fprintf('Y = %f \n', full(hdot(Trajectory{end,1},Operators.Ly_State{1})));
fprintf('Z = %f \n', full(hdot(Trajectory{end,1},Operators.Lz_State{1})));

VECTOR_FIRST_SPIN = [full(hdot(Trajectory{end,1},Operators.Lx_State{1}));
                     full(hdot(Trajectory{end,1},Operators.Ly_State{1}));
                     full(hdot(Trajectory{end,1},Operators.Lz_State{1}))]; 
VECTOR_FIRST_SPIN = VECTOR_FIRST_SPIN / norm(VECTOR_FIRST_SPIN);
VECTOR_FIRST_SPIN

fprintf('********************************************************* \n');
fprintf('********************************************************* \n');
fprintf('********************************************************* \n');


%% Write the shape file associated to the newly born pulse

fprintf('\n\n\n')
choice = input('Do you want to write the shape file associated to this solution (yes/no) ? \n','s');
if strcmp(choice,'yes')
    fname = input('Name of the shape file ? \n','s');
    if isempty(fname)
        fname = 'Shape_file.txt';
    end
    SPINACH_Shape_File_Writer(waveform,fname,'yes');
end


end