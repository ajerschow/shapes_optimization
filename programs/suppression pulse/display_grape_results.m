% This function has been written to lighten the main file performing the
% grape pulse optimization
%
% It allows to monitor :
%
% (*)       The shape of the waveform
% (**)      The relative population of the Ix, Iy and Iz states
%           This function can also, if requested, returns as an OUTPUT 
%           these relative populations at each timestep
%
% INPUT
% 
%           Trajectory        : density matrix at each timestep, can be a cell or
%                               an array/ a matrix
%
%           ctrl_sys          : the parameters needed by grape.m to function,
%                               basically the number of steps, the timestep ...
%
%           SpinX_ State & co : set of control operators
%
%           n                 : indicates the index of the figure in which
%                               the plot will be displayed


function [Population_Ix,Population_Iy,Population_Iz] = display_grape_results(Trajectory,ctrl_sys,waveform,...
                            Operators,n)
    
    % Initialization of the outputs
    Population_Ix = zeros(1,ctrl_sys.nsteps+1);
    Population_Iy = zeros(1,ctrl_sys.nsteps+1);
    Population_Iz = zeros(1,ctrl_sys.nsteps+1);
    
    choice = class(Trajectory);

    switch(choice)
        
        case 'cell'  % when function called *outside* of the grape function
            
            for k=1:(ctrl_sys.nsteps+1) 
                Population_Ix(1,k) = (hdot(Trajectory{k,1},Operators.SpinX_State))/...
                                      hdot(Operators.SpinX_State,Operators.SpinX_State_NORMALIZATION);
                Population_Iy(1,k) = (hdot(Trajectory{k,1},Operators.SpinY_State))/...
                                      hdot(Operators.SpinY_State,Operators.SpinY_State_NORMALIZATION);
                Population_Iz(1,k) = (hdot(Trajectory{k,1},Operators.SpinZ_State))/...
                                      hdot(Operators.SpinZ_State,Operators.SpinZ_State_NORMALIZATION);
            end
            
        case 'double' % when function called *inside* the grape function
            
            for k=1:(ctrl_sys.nsteps+1) 
                Population_Ix(1,k) = (hdot(Trajectory(:,k),Operators.SpinX_State))/...
                                      hdot(Operators.SpinX_State,Operators.SpinX_State_NORMALIZATION);
                Population_Iy(1,k) = (hdot(Trajectory(:,k),Operators.SpinY_State))/...
                                      hdot(Operators.SpinY_State,Operators.SpinY_State_NORMALIZATION);
                Population_Iz(1,k) = (hdot(Trajectory(:,k),Operators.SpinZ_State))/...
                                      hdot(Operators.SpinZ_State,Operators.SpinZ_State_NORMALIZATION);
            end 
            
        otherwise
            
            error('The first argument of the function should be either a matrix either a cell');
            
    end

    figure(n)
    subplot(1,2,1)
    plot(linspace(0,ctrl_sys.time_step*ctrl_sys.nsteps,ctrl_sys.nsteps),waveform);
    title('control amplitudes'); 
    xlabel('time, seconds');
    ylabel('control amplitude'); 
    legend('X','Y')
    axis tight;
    drawnow;
    subplot(1,2,2)
    plot(...
         1:(ctrl_sys.nsteps+1),real(Population_Ix(1,:)),'r',...
         1:(ctrl_sys.nsteps+1),real(Population_Iy(1,:)),'b',...
         1:(ctrl_sys.nsteps+1),real(Population_Iz(1,:)),'m' ...
        );
    xlabel('number of steps');
    drawnow;
    
end