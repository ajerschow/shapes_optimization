% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
%
% This is the modified version of the grape function !
%
% What are the modifications ?
%
%           * This modified version allows to minimize the population of
%             one state
%             ex: if you set the target as being <Ix>, this version allows
%             you to have <Ix> = 0 (and not +1/-1 as in the original grape
%             function)
%
%           * How ?
%             New cost function  & new gradient
%             cost_new = (cost_grape)^2
%             grad_new = 2 * cost_grape * grad_grape
%
%           * What has been modified ?
%
%                   (1) One new input in the main function. The new input
%                       is 'version'. It has to be a string. Two choices :
%                       'original' or 'modified' which allows to specify
%                       wether you want to use the original version of the 
%                       grape function or the modified one.
%
%                   (2) In the main function, I added at line 337 the
%                       choice between the two cost functions (called
%                       fidelity here)
%
%                   (3) I modified the grape_grad function so that we can
%                       compute the gradients associated with the new cost
%                       function.
%                       Note = Between the different choices available to
%                       compute the gradients, only 'auxmat' has been
%                       tested !
%
% Note: An error message is displayed if the user want to use the hessians
%       matrixes with the modified version of the grape function. Why ?
%       Simply, because I have not modified this part of the grape
%       function.
%
% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\
% /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\


% Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
% and Hessian. Propagates the system through a user-supplied shaped pulse
% from a given initial state and projects the result onto the given final
% state. The real part of the projection is returned, along with its
% gradient and Hessian with respect to amplitudes of all operators in every
% time step of the shaped pulse.
%
% <http://spindynamics.org/wiki/index.php?title=Grape.m>

function [diag_data,fidelity,total_grad,total_hess]=...
    grape_modified_version(spin_sys,ctrl_sys,drift,waveform_in,version)       %#ok<*PFBNS>

% /!\ /!\ /!\ /!\ /!\ /!\
%  NEW INPUT --> version
% /!\ /!\ /!\ /!\ /!\ /!\
    
% Allow numeric input for single source calculations
if isnumeric(ctrl_sys.rho)
    rho={ctrl_sys.rho};
else
    rho=ctrl_sys.rho;
end

% Allow numeric input for single target calculations
if isnumeric(ctrl_sys.target)
    target={ctrl_sys.target};
else
    target=ctrl_sys.target;
end

% Allow numeric input for single control calculations
if isnumeric(ctrl_sys.control_ops)
    controls={ctrl_sys.control_ops};
else
    controls=ctrl_sys.control_ops;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if the user forgets to specify which version of the grape function, he
% wants to use, choose automatically the original version #OG4LIFE

if nargin <= 4
    version = 'original';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% if offset range is set - define variables
if isfield(ctrl_sys,'offset_ops')
    if isnumeric(ctrl_sys.offset_ops)
        offset_ops={ctrl_sys.offset_ops};
    else
        offset_ops=ctrl_sys.offset_ops;
    end
else
    offset_ops=eye(size(controls{1}));
end

% Preallocate results
fidelity=0; total_grad=zeros(size(waveform_in));
total_hess=zeros(size(waveform_in,1)*size(waveform_in,2));
current_state_array=cell(1,numel(target));

% Silence spinach output for propagator diagnostics
outtemp=spin_sys.sys.output;
spin_sys.sys.output='hush';

% Loop over power level miscalibrations
for pwr_ind=1:length(ctrl_sys.power_level)
    
    % Scale the waveform
    waveform=ctrl_sys.power_level(pwr_ind)*waveform_in;
    
    % Loop over offsets
    for offset_ind=1:length(ctrl_sys.offset_vals)
        
        % Loop over sources and targets
        for contribution=1:numel(rho)
            
            % Grab the current source and target
            init_state=rho{contribution};
            target_state=target{contribution};
            
            % Preallocate trajectory arrays
            fwd_traj=zeros(size(init_state,1),ctrl_sys.nsteps+1);
            bwd_traj=zeros(size(init_state,1),ctrl_sys.nsteps+1);
            
            % Run the forward propagation
            fwd_traj(:,1)=init_state;
            for n=1:ctrl_sys.nsteps
                
                % Add current controls
                L=drift;
                for k=1:length(controls)
                    L=L+waveform(k,n)*controls{k};
                end
                
                % Add the offset contribution
                if length(ctrl_sys.offset_vals)>1
                    for m=1:length(offset_ops)
                        L=L+ctrl_sys.offset_vals(offset_ind)*offset_ops{m};
                    end
                end
                fwd_traj(:,n+1)=step(spin_sys,L,fwd_traj(:,n),ctrl_sys.time_step);
            end
            
            % Run the backward propagation
            bwd_traj(:,ctrl_sys.nsteps+1)=target_state;
            for n=ctrl_sys.nsteps:-1:1
                
                % Add current controls
                L=drift';
                for k=1:length(controls)
                    L=L+waveform(k,n)*controls{k};
                end
                
                % Add the offset contribution
                if length(ctrl_sys.offset_vals)>1
                    for m=1:length(offset_ops)
                        L=L+ctrl_sys.offset_vals(offset_ind)*offset_ops{m};
                    end
                end
                bwd_traj(:,n)=step(spin_sys,L,bwd_traj(:,n+1),-ctrl_sys.time_step);
            end
            
            % Compute the gradient
            if nargout==3
                
                % Loop over time steps in parallel
                parfor n=1:ctrl_sys.nsteps
                    
                    % Initialize the Liouvillian
                    L=drift;
                    
                    % Add current controls
                    for k=1:length(controls)
                        L=L+waveform(k,n)*controls{k};
                    end
                    
                    % Add the offset contribution
                    if length(ctrl_sys.offset_vals)>1
                        for m=1:length(offset_ops)
                            L=L+ctrl_sys.offset_vals(offset_ind)*offset_ops{m};
                        end
                    end
                    
                    % Calculate gradient at this timestep
                    grad(:,n)=grape_grad(spin_sys,L,controls,waveform,fwd_traj,bwd_traj,ctrl_sys.time_step,n,target_state,version);
                    
                end
                
                % Add the normalised gradient to the total
                total_grad=total_grad+ctrl_sys.power_level(pwr_ind)*real(grad);
                
                % Compute the gradient and hessian
            elseif nargout>3
                
                if strcmp(version,'modified')
                   error('The modified version of the grape function does not handle the hessian yet ...\nCome back later!\n\n'); 
                end
                
                % initialise cell structure for propagators and their derivatives
                props=cell(1,ctrl_sys.nsteps);
                props_D=cell(length(controls),ctrl_sys.nsteps);
                props_DD=cell(length(controls),length(controls),ctrl_sys.nsteps);
                
                % propagator and derivative calculations
                parfor n=1:ctrl_sys.nsteps
                    
                    % Initialize the Liouvillian at this timestep
                    L=drift;
                    
                    % Add all current controls to Liouvillian for the current timestep
                    for k=1:length(controls)
                        L=L+waveform(k,n)*controls{k};
                    end
                    
                    % Add the offset contribution
                    if length(ctrl_sys.offset_vals)>1
                        for m=1:length(offset_ops)
                            L=L+ctrl_sys.offset_vals(offset_ind)*offset_ops{m};
                        end
                    end
                    
                    % Compute and store propagators and their derivatives
                    [props(1,n),props_D(:,n),props_DD(:,:,n)]=grape_hess(spin_sys,L,controls,ctrl_sys.time_step,fwd_traj(:,n));
                    
                end
                
                %initialise gradient vector and Hessian matrix
                grad=zeros(length(controls),ctrl_sys.nsteps);
                hess=cell(ctrl_sys.nsteps,ctrl_sys.nsteps);                   hess(:,:)={zeros(length(controls))};
                
                % split the hessian into two, to aid parfor over ctrl_system.nsteps/2
                hess_upper=cell(ceil(ctrl_sys.nsteps/2),ctrl_sys.nsteps);     hess_upper(:,:)={zeros(length(controls))};
                hess_lower=cell(floor(ctrl_sys.nsteps/2),ctrl_sys.nsteps);	hess_lower(:,:)={zeros(length(controls))};
                
                parfor n_prime=1:ceil(ctrl_sys.nsteps/2)
                    
                    % Preallocate space for parfor use
                    grad_k=zeros(length(controls),ctrl_sys.nsteps);
                    hess_k=zeros(length(controls));
                    
                    % create a hessian block vector to aid parfor loop
                    hess_vec_l=cell(1,ctrl_sys.nsteps);      hess_vec_l(:,:)={zeros(length(controls))};
                    hess_vec_u=cell(1,ctrl_sys.nsteps);      hess_vec_u(:,:)={zeros(length(controls))};
                    
                    for m_prime=1:ctrl_sys.nsteps+1
                        
                        % calculate n and m proper from primed coordinates
                        m=m_prime+n_prime-1;
                        if m>ctrl_sys.nsteps
                            m=1+ctrl_sys.nsteps-mod(m,ctrl_sys.nsteps);
                            n=ctrl_sys.nsteps-n_prime+1;
                        else
                            n=n_prime;
                        end
                        
                        for k=1:length(controls)
                            for j=1:length(controls)
                                
                                % <sigma|UN...Un+1 D^2(n) Un...U1|rho>
                                if m==n && j==k
                                    
                                    % Calculate gradient at this timestep
                                    grad_k(k,n)=bwd_traj(:,n+1)'*(props_D{k,n}*fwd_traj(:,n));
                                    
                                    % calculate KxK hessian diagonal
                                    hess_k(k,j)=bwd_traj(:,n+1)'*props_DD{k,j,n}*fwd_traj(:,n);
                                    
                                elseif m==n
                                    % calculate KxK hessian off-diagonal
                                    hess_k(k,j)=bwd_traj(:,n+1)'*props_DD{k,j,n};
                                    
                                elseif n<m, N=[n m]; K=[k j]
                                    
                                    % Calculate and store initial derivative(1)*forward_traj
                                    forward_store=props_D{K(1),N(1)}*fwd_traj(:,N(1));
                                    
                                    % loop over timesteps between derivatives
                                    for n_ind=N(1)+1:N(2)-1
                                        
                                        % update forward_traj between derivatives
                                        forward_store=props{1,n_ind}*forward_store;
                                    end
                                    
                                    % Calculate and derivative(2)*forward_traj
                                    forward_store=props_D{K(2),N(2)}*forward_store;
                                    
                                    % Calculate Hessian element
                                    hess_k(k,j)=bwd_traj(:,N(2)+1)'*forward_store;
                                end
                            end
                        end
                        
                        % Allocate the Hessian elements to parfor vector
                        if n>ceil(ctrl_sys.nsteps/2)
                            hess_vec_l{1,m}=hess_k;
                        else
                            hess_vec_u{1,m}=hess_k;
                        end
                    end
                    
                    % allocate hessian parfor vector to hessian blocks
                    hess_lower(n_prime,:)=hess_vec_l;
                    hess_upper(n_prime,:)=hess_vec_u;
                    
                    % Allocate the gradient over controls at this timestep
                    grad=grad+grad_k;
                    
                end
                
                % construct hessian proper from upper and lower blocks
                hess(1:ceil(ctrl_sys.nsteps/2),:)=hess_upper;
                hess(1+ceil(ctrl_sys.nsteps/2):end,:)=hess_lower(end:-1:1,:);
                
                % Add the hessian to the total
                total_hess=total_hess+(ctrl_sys.power_level(pwr_ind)^2)*real(cell2mat(hess));
                
                % Add the gradient to the total
                total_grad=total_grad+ctrl_sys.power_level(pwr_ind)*real(grad);
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the cost function and add to the total
            
            if strcmp(version,'original')
                fidelity=fidelity+real(target_state'*fwd_traj(:,end));
            elseif strcmp(version,'modified')
                fidelity=fidelity+((real(target_state'*fwd_traj(:,end)))^2);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the cost function and add to the total
            % fidelity=fidelity+real(target_state'*fwd_traj(:,end));
            
            % store the current state for zero offset
            if ctrl_sys.offset_vals(offset_ind)==0 && pwr_ind==(1+length(ctrl_sys.power_level))/2
                current_state_array(contribution)={fwd_traj(:,end)};
            end
        end
        
    end
end

% normalisation constant
norm_const=numel(rho)*length(ctrl_sys.offset_vals)*length(ctrl_sys.power_level);

% Normalize the total objective, total gradient, and total Hessian
fidelity=fidelity/(norm_const);
total_grad=total_grad/(norm_const);
total_hess=total_hess/(norm_const);

% create block-diagonal blanking matrix
offdiag_blank=kron(eye(size(waveform,2)),ones(size(waveform,1)));
diag_blank=~offdiag_blank;

% Force Hessian symmetry and fill empty Hessian entries
total_hess=(total_hess+total_hess').*diag_blank+...
    (total_hess+total_hess').*offdiag_blank/2;

% re-enable spinach console output
spin_sys.sys.output=outtemp;

% Compile diagnostics data
diag_data.current_state=current_state_array;
diag_data.rho=rho;
diag_data.target=target;
diag_data.total_objective=fidelity;
diag_data.spin_system=spin_sys;
diag_data.ctrl_system.power_level=ctrl_sys.power_level;
diag_data.trajectory=fwd_traj;
diag_data.total_grad=total_grad;
diag_data.total_hess=total_hess;
diag_data.dt=ctrl_sys.time_step;
diag_data.controls=controls;
diag_data.waveform=waveform;
diag_data.ctrl_system.nsteps=ctrl_sys.nsteps;
diag_data.drift=drift;

end

function [props,props_d,props_dd]=grape_hess(spin_system,L,ctrls,dt,fwd_traj)

% Preallocate propagator derivatives
props_d=cell(1,length(ctrls));
props_dd=cell(length(ctrls),length(ctrls));

switch spin_system.tols.dP_method
    
    case 'auxmat'
        
        % loop over K*K controls
        for k=1:length(ctrls)
            for j=1:length(ctrls)
                
                if k==j
                    %find explicit propagator to reuse for off-diagonal
                    prop_dirdiff=dirdiff(spin_system,L,{ctrls{k},ctrls{k}},dt,3);
                    
                    % store 1st and 2nd order propagator derivatives
                    props_d(k)=prop_dirdiff(2); props_dd(k,k)=prop_dirdiff(3);
                else
                    % Create the auxiliary matrix
                    aux_matrix=[L, ctrls{k}, 0*L; 0*L, L, ctrls{j}; 0*L, 0*L, L];
                    
                    % Propagate the auxiliary vector
                    aux_vector=step(spin_system,aux_matrix,[zeros(size(fwd_traj)); zeros(size(fwd_traj)); fwd_traj],dt);
                    
                    % Store propagated second order derivative vector
                    props_dd{k,j}=2*aux_vector(1:(end/3));
                end
            end
        end
        
        % Store propagator of Liouvillian, L.
        props=prop_dirdiff(1);
        
    otherwise
        
        error('grape: unknown second order differentiation method.');
end

end

function control_derivatives=grape_grad(spin_system,L,ctrls,wf,fwd_traj,bwd_traj,dt,n,target_state,version)

% Preallocate derivatives
control_derivatives=zeros(length(ctrls),1);

% Generate control derivatives
switch spin_system.tols.dP_method
    
    case 'auxmat'
        
        % Loop over controls
        for k=1:length(ctrls)
            
            % Create the auxiliary matrix
            aux_matrix=[L, ctrls{k}; 0*L, L];
            
            % Propagate the auxiliary vector
            aux_vector=step(spin_system,aux_matrix,[zeros(size(fwd_traj(:,n))); fwd_traj(:,n)],dt);
            
            % Compute the derivative
            control_derivatives(k)=bwd_traj(:,n+1)'*aux_vector(1:(end/2));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(version,'modified')
                % Modified gradients
                control_derivatives(k) = 2*hdot(target_state,fwd_traj(:,end))*control_derivatives(k);  
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
        
    case 'fd_O(h^2)'
        
        % Get the finite difference step
        delta=max(abs(wf(:)))*sqrt(eps);
        
        % Loop over controls
        for k=1:length(ctrls)
            
            % Compute control derivatives with central finite differences
            control_derivatives(k)=bwd_traj(:,n+1)'*(step(spin_system,L+delta*ctrls{k},fwd_traj(:,n),dt)-...
                step(spin_system,L-delta*ctrls{k},fwd_traj(:,n),dt))/(2*delta);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if strcmp(version,'modified')  
            error('MODIFIED GRAPE : spin_system.tols.dP_method = fd_0 (h^2) - This part has not been modified yet')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    case 'fd_O(h^4)'
        
        % Get the finite difference step
        delta=max(abs(wf(:)))*sqrt(eps);
        
        % Loop over controls
        for k=1:length(ctrls)
            
            % Compute control derivatives with central finite differences
            control_derivatives(k)=bwd_traj(:,n+1)'*(-1*step(spin_system,L+2*delta*ctrls{k},fwd_traj(:,n),dt)...
                +8*step(spin_system,L+delta*ctrls{k},fwd_traj(:,n),dt)...
                -8*step(spin_system,L-delta*ctrls{k},fwd_traj(:,n),dt)...
                +1*step(spin_system,L-2*delta*ctrls{k},fwd_traj(:,n),dt))/(12*delta);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        if strcmp(version,'modified')  
            error('MODIFIED GRAPE : spin_system.tols.dP_method = fd_0 (h^4) This part has not been modified yet')          
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    otherwise
        
        % Complain and bomb out
        error('unknown differentiation method.');
        
end


end

% In any culture, subculture, or family in which belief is valued above
% thought, self-surrender is valued above self-expression, and conformity
% is valued above integrity, those who preserve their self-esteem are
% likely to be heroic exceptions.
%
% Nathaniel Branden


