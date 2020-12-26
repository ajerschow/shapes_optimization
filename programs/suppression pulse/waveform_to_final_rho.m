% Starting from a waveform and a final state, this function returns the
% intermediate density matrixes and the final one
% Formalism taken from the grape.m function. 


function Trajectory = waveform_to_final_rho(spin_system, L, ctrl_sys, waveform)
    
    waveform_1 = waveform * ctrl_sys.power_level;
    
    Trajectory = cell(ctrl_sys.nsteps+1,1);
    Trajectory{1,1} = ctrl_sys.rho;
    
    for n = 1:ctrl_sys.nsteps
        H = L;
        for control = 1:length(ctrl_sys.control_ops)
            H = H + waveform_1(control,n)*ctrl_sys.control_ops{control};
        end
        Trajectory{n+1,1} = step(spin_system,H,Trajectory{n,1},ctrl_sys.time_step);
    end 
    
end