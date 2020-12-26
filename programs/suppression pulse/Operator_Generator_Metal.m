
function [Operators] =  Operator_Generator_Metal(spin_system,L_rho,phi)

% L_rho = length(rho) for the intizialization of the matrices

%% Number of spins in the simulations
Nbr_Spins = numel(spin_system.logs.sys.isotopes);

%% Control and State operators - Normalization factor

% Control Operators - Pre-allocation
Lx = cell(1,Nbr_Spins);
Ly = cell(1,Nbr_Spins);
Lz = cell(1,Nbr_Spins);
SpinX = spalloc(L_rho,1,L_rho);
SpinY = spalloc(L_rho,1,L_rho);
SpinZ = spalloc(L_rho,1,L_rho);

% State operators - Pre-allocation
Lp_plus_state = cell(1,Nbr_Spins);
Lp_minus_state = cell(1,Nbr_Spins);
Lx_State = cell(1,Nbr_Spins);
Ly_State = cell(1,Nbr_Spins);
Lz_State = cell(1,Nbr_Spins);
SpinX_State = zeros(L_rho,1);
SpinY_State = zeros(L_rho,1);
SpinZ_State = zeros(L_rho,1);

% Normalization - Pre-allocation
Lx_State_NORMALIZATION = cell(1,Nbr_Spins);
Ly_State_NORMALIZATION = cell(1,Nbr_Spins);
Lz_State_NORMALIZATION = cell(1,Nbr_Spins);
SpinX_State_NORMALIZATION = zeros(L_rho,1);
SpinY_State_NORMALIZATION = zeros(L_rho,1);
SpinZ_State_NORMALIZATION = zeros(L_rho,1);

% Coil
Coil_X = zeros(L_rho,1);
Coil_Y = zeros(L_rho,1);
Coil_Z = zeros(L_rho,1);

for kk=1:Nbr_Spins
    kk
    % Spinach functions to create the operators and states
    Lp = operator(spin_system,{'L+'},{kk});
    Lp_plus_state{kk} = state(spin_system,{'L+'},{kk});
    Lp_minus_state{kk} = state(spin_system,{'L-'},{kk});
    
    % control operator - individual 
    X = ((Lp+Lp')/2);
    Y = ((Lp-Lp')/2i);
    
    Lx{kk} = (X*cos(phi(kk))+Y*sin(phi(kk)))*exp(-phi(kk));
    Ly{kk} = (Y*cos(phi(kk))-X*sin(phi(kk)))*exp(-phi(kk));
    Lz{kk} = (operator(spin_system,{'Lz'},{kk})); %no attenuation
    
    % control operator - global
    SpinX = SpinX + Lx{kk};
    SpinY = SpinY + Ly{kk};
    SpinZ = SpinZ + Lz{kk};
    
    % state operator --> to work in the liouville state
    XX = ((Lp_plus_state{kk} + Lp_minus_state{kk})/2);
    YY = ((Lp_plus_state{kk} - Lp_minus_state{kk})/(2i));
    
    % state operator - individual
    Lx_State{kk} = (XX*cos(phi(kk))-YY*sin(phi(kk)))*exp(-phi(kk)); 
    Ly_State{kk} = (YY*cos(phi(kk))+XX*sin(phi(kk)))*exp(-phi(kk)); 
    Lz_State{kk} = (state(spin_system,{'Lz'},{kk}))*exp(-phi(kk));
    
    % state operator - global
    SpinX_State = SpinX_State + Lx_State{kk};
    SpinY_State = SpinY_State + Ly_State{kk};
    SpinZ_State = SpinZ_State + Lz_State{kk};
    
    % state operator - Normalization
    Lx_State_NORMALIZATION{kk} = (XX*cos(phi(kk))-YY*sin(phi(kk))); 
    Ly_State_NORMALIZATION{kk} = (YY*cos(phi(kk))+XX*sin(phi(kk))); 
    Lz_State_NORMALIZATION{kk} = (state(spin_system,{'Lz'},{kk}));
    SpinX_State_NORMALIZATION = SpinX_State_NORMALIZATION + Lx_State_NORMALIZATION{kk};
    SpinY_State_NORMALIZATION = SpinY_State_NORMALIZATION + Ly_State_NORMALIZATION{kk};
    SpinZ_State_NORMALIZATION = SpinZ_State_NORMALIZATION + Lz_State_NORMALIZATION{kk};
    
    %Coil --> not shifted !
    Coil_X = Coil_X + XX;
    Coil_Y = Coil_Y + YY;
    Coil_Z = Coil_Z + state(spin_system,{'Lz'},{kk});
    
end

%% Output

% Control operator
Operators.Lx = Lx;
Operators.Ly = Ly;
Operators.Lz = Lz;
Operators.SpinX = SpinX;
Operators.SpinY = SpinY;
Operators.SpinZ = SpinZ;

% State operators 
Operators.Lp_plus_state = Lp_plus_state;
Operators.Lp_minus_state = Lp_minus_state;
Operators.Lx_State = Lx_State;
Operators.Ly_State = Ly_State;
Operators.Lz_State = Lz_State;
Operators.SpinX_State = SpinX_State;
Operators.SpinY_State = SpinY_State;
Operators.SpinZ_State = SpinZ_State;

% Normalization - Pre-allocation
Operators.Lx_State_NORMALIZATION = Lx_State_NORMALIZATION;
Operators.Ly_State_NORMALIZATION = Ly_State_NORMALIZATION;
Operators.Lz_State_NORMALIZATION = Lz_State_NORMALIZATION;
Operators.SpinX_State_NORMALIZATION = SpinX_State_NORMALIZATION;
Operators.SpinY_State_NORMALIZATION = SpinY_State_NORMALIZATION;
Operators.SpinZ_State_NORMALIZATION = SpinZ_State_NORMALIZATION;

% Coil
Operators.Coil.Coil_X = Coil_X;
Operators.Coil.Coil_Y = Coil_Y;
Operators.Coil.Coil_Z = Coil_Z;

end