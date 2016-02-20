clearvars;
syn_params = struct(...
    'Tampa', 3.5,'Eampa', 0,...
    'Tnmda1', 80,'Tnmda2', 20,'Enmda', 0,...
    'Tgaba', 5.5,'Egaba', -65);
params = struct(...
    'mem', struct('Er',-65,'Et',-40,'tau', 20,'Vr',-70,'Tref', 3),...
    'syn', syn_params);
rate_e = 0.40;
rate_i = 0.40;

M = 10;
neur = cell(M, 1);
for i = 1:M
    neur{i} = make_neuron(params);
end
neur = [neur{:}];
N = 1000;
t0 = 0:N;
ts = [];
Y0 = [-50, 0.05*rand(1,4)]';
SOL = zeros(N+1, length(Y0));
SOL(1,:) = Y0;
opts = odeset('Events',neur.events,'RelTol',1e-3,'AbsTol',1e-3,...
    'InitialStep', 1,'MaxStep', 1);
exit = false;
Tend = 1;
for i = 1:N
    [T, Y, TE, YE] = ode45(neur.odefun, [t0(i), t0(i+1)], Y0, opts);
    if ~isempty(TE)
        ts = [ts; TE];
        Y = YE(end,:)';
        Y(1) = params.mem.Vr;
        neur.set_spike(TE);
        [T, Y] = ode15s(neur.odefun, [TE, t0(i+1)], Y, opts);
    end
    SOL(i+1,:) = Y(end,:);
    % update for next iteration
    Y0 = Y(end,:);
    spike_e = poissrnd(rate_e);
    Y0(2) = Y0(2) + 0.02 * spike_e;
    Y0(5) = Y0(5) + 0.01 * spike_e;
    spike_i = poissrnd(rate_i);
    Y0(3) = Y0(3) + 0.04*spike_i;
    
    Y0(2:end) = min(Y0(2:end), repmat(0.2, 1, length(Y0)-1));
end

time = t0;
subplot(211); hold on
plot(time, SOL(:,1),'b')
stem(ts, 20 *ones(size(ts)),'-b');
subplot(212);
plot(time, SOL(:,2:end));