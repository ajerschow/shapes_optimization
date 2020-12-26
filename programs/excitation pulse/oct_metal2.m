% Boris Kharkov
% Developped for Spinach 1.6
%

% % In this simulation we find optimal pulse for transfering z magmetization
% into x magnetization in a metal.
% 

% nsteps is the number of steps in the waveform;
% flip_initial defines the pulse duration in connection to the maximal
% pulse amplitude
% for example, run oct_metal2(50, 0.8)

function oct_metal2(nsteps, flip_initial)

% Magnetic field
sys.magnet=11.7434;

% Isotopes
nOfSpins = 100;
sys.isotopes = cell(1, nOfSpins);
sys.isotopes(:) = {'1H'};

%% Spinach housekeeping
% Relaxation theory
inter.relaxation='none';

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=1; 

% Run Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

H=hamiltonian(assume(spin_system,'nmr'));

%% States and operators specification
% Set up relevant states
Lp_ = cell(1,nOfSpins); Lm_ = cell(1,nOfSpins); Lz_ = cell(1,nOfSpins);
Lx_ = cell(1,nOfSpins); Ly_ = cell(1,nOfSpins);
Lp = cell(1,nOfSpins); Lz = cell(1,nOfSpins);
Lx = cell(1,nOfSpins); Ly = cell(1,nOfSpins);
for k = 1:nOfSpins
    Lp_{k} = state(spin_system,{'L+'},{k});
    Lm_{k} = state(spin_system,{'L-'},{k});
    Lz_{k} = state(spin_system,{'Lz'},{k});
    Lx_{k} = (Lp_{k} + Lm_{k})/2;
    Ly_{k} = (Lp_{k} - Lm_{k})/2i;
    
    Lz{k} = operator(spin_system,{'Lz'},{k});
    Lp{k} = operator(spin_system,{'L+'},{k});
    Lx{k} = (Lp{k} + Lp{k}') / 2; Ly{k} = (Lp{k} - Lp{k}') / 2i;
end
Lp_total_ = state(spin_system,'L+','all');Lm_total_ = state(spin_system,'L-','all');
Lz_total_ = state(spin_system,'Lz','all');
Lx_total_ = (Lp_total_ + Lm_total_)/2;
Ly_total_ = (Lp_total_ - Lm_total_)/2i;

% Get the control operators
Lp_total=operator(spin_system,'L+','1H');
Lx_total=(Lp_total+Lp_total')/2; Ly_total=(Lp_total-Lp_total')/2i;

%% Set up source and target states
attenuation = 1;
range = 2*pi;
x = linspace(0, range, nOfSpins);
amps = exp(- x / attenuation);

targetX_ = zeros(size(Lz_total_));
idealX_ = zeros(size(Lz_total_));
targetY_ = zeros(size(Lz_total_));
Lx_control = zeros(size(Lx_total));
Ly_control = zeros(size(Lx_total));
for k = 1:nOfSpins
    targetX_ = targetX_ + amps(k) * (Lx_{k} * cos(x(k)) - Ly_{k} * sin(x(k)));
    targetY_ = targetY_ + amps(k) * (Ly_{k} * cos(x(k)) + Lx_{k} * sin(x(k)));
    Lx_control = Lx_control + amps(k) * (Lx{k} * cos(x(k)) + Ly{k} * sin(x(k)));
    Ly_control = Ly_control + amps(k) * (Ly{k} * cos(x(k)) - Lx{k} * sin(x(k)));
    idealX_ = idealX_ + (Lx_{k} * cos(x(k)) - Ly_{k} * sin(x(k)));
end

sources = {Lz_total_};
targets = -1 * {targetX_};  % Here, -1 because in spinach oct searches the minimum of  <-Lx_, rho>, which is equivalent to finding a maximum of <Lx_, rho>. And thats what we need.

%% Simulation
pl_Hz = 1000;
% nsteps = 20;
pulse_duration = flip_initial / pl_Hz;
time_step = pulse_duration / nsteps;
power_level = -1*2*pi * pl_Hz;   % -1 because effective nutation frecuency will be negative for 1H: for instance, w0 = -1 * gamma * Iz
guess = [zeros(1,nsteps);...
         zeros(1,nsteps)];
% guess = rand(2,nsteps);
controls = {Lx_control, Ly_control};

pulseq.power_level = power_level;
pulseq.controls = {Lx_control, Ly_control};
pulseq.time_grid = time_step * ones(1,nsteps);


% Run the optimization
options = struct('method','lbfgs','max_iterations',1000);
resultingWaveform=fminnewton(spin_system,@cost_function,guess,options);

pulseq.waveform = resultingWaveform;

waveform_ampph = ampsPhases(resultingWaveform);

rhos = applyPulSeq(spin_system, H, Lz_total_, pulseq, 'trajectory');
fid = calcFID(rhos,{targetX_,targetY_})  * hdot(targetX_, targetX_) / hdot(idealX_, targetX_);

% calculate and plot spin orientation in the sample
rho_final = rhos{size(rhos, 2)};
spaceDistr = zeros(1, nOfSpins);
for k = 1:nOfSpins
    spinDetectOp = (Lx_{k} * cos(x(k)) - Ly_{k} * sin(x(k)));
    spaceDistr(k) = hdot(spinDetectOp, rho_final) / hdot(spinDetectOp, spinDetectOp);
end

%% plotting section

figure(2);
plot(real(spaceDistr)); hold on;
plot(real(spaceDistr) .* amps);
plot(real(amps));

% Plot the fid
figure(100)
subplot(1,2,1); plot(real(fid(1,:))); 
title('FID X component'); axis tight;
subplot(1,2,2); plot(real(fid(2,:)));
title('FID Y component'); axis tight;

% Plot the pulse amplitudes
figure(300)
subplot(1,2,1); plot(resultingWaveform(1,:)); 
title('Pulse X component'); axis tight;
subplot(1,2,2); plot(resultingWaveform(2,:));
title('Pulse Y component'); axis tight;

% Plot the pulse amplitudes and phase
figure(301)
subplot(1,2,1); plot(waveform_ampph(1,:)); 
title('Pulse amplitude'); axis([0 nsteps 0 2]);% axis tight;
subplot(1,2,2); plot(waveform_ampph(2,:));
title('Pulse phase'); axis([0 nsteps -180 180]);

%% save to file
% filename = ['results/octMetal2' '_' num2str(nsteps) '_' num2str(pulse_duration*pl_Hz) '.mat'];
% save(filename);

%% cost functional is here
    % Set the cost functional
    function [cost,cost_grad]=cost_function(waveform)
        
        % Preallocate the outputs
        cost=0; cost_grad=zeros(size(waveform));
            
            % Calculate the objective and its gradient
            [diag_data,current_cost,current_cost_grad]=grape(spin_system,H,controls,waveform,time_step,nsteps,sources,targets,power_level);
            
            %plot(current_cost_grad(1,:))
            
            % Penalize excursions outside the power envelope
            [current_pen,current_pen_grad]=penalty(waveform,'mean_square_spillout',-ones(size(waveform)),ones(size(waveform)));
            
            cost=cost+current_cost+current_pen;
            cost_grad=current_cost_grad+current_pen_grad;
            
            current = pulseq;
            current.waveform = waveform;
            rhos1 = applyPulSeq(spin_system, H, Lz_total_, current, 'trajectory');
            fid1 = calcFID(rhos1,{targetX_,targetY_})  * hdot(targetX_, targetX_) / hdot(idealX_, targetX_);

        % Plot the current trajectory
        subplot(1,2,1); plot(linspace(0,time_step*nsteps,nsteps+1),real(fid1));
        title('trajectories'); xlabel('time, seconds');
        ylabel('amplitude'); axis tight; drawnow;

        % Plot the current waveform
        subplot(1,2,2); plot(linspace(0,time_step*nsteps,nsteps),waveform);
        title('control amplitudes'); xlabel('time, seconds');
        ylabel('control amplitude'); axis tight; drawnow;
        
    end

end

function ret = ampsPhases(fid)
    ret = zeros(size(fid));
    ret(1,:) = (fid(1,:).^2 + fid(2,:).^2).^0.5;
    ret(2,:) = (atan(fid(2,:) ./ fid(1,:))/pi*180) + 180 * (fid(1,:) < 0).*(sign(fid(2,:)));
end

