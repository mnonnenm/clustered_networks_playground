function [spikes, setup] = run_network(pars, sf)

    if nargin < 2
        sf = ''; % savefile, with absolute or relative parth to destination
    end
    
    if strcmp(sf, '')
        if_save = false;
    else 
        if_save = true;
    end
    
    neur = cell(pars.N, 1);
    for i = 1:pars.N
        if i <= pars.Ne
            neur{i} = make_neuron(pars.params_e);
        else
            neur{i} = make_neuron(pars.params_i);
        end
    end
    neur = [neur{:}];


    t0 = 0:pars.dt:pars.T;
    Y0 = arrayfun(@(~) {[-65 + rand(1) * 10,0.05*rand(1,4),20]'},1:pars.N);
    

    opts = cell(pars.N,1);
    for i = 1:pars.N
        opts{i} = odeset(...
            'Events',neur(i).events,...
            'RelTol',1e-3,'AbsTol',1e-3,...
            'InitialStep', t0(2) - t0(1),'MaxStep', t0(2) - t0(1));
    end
    
    if if_save
        spikes = zeros(pars.chunksize, pars.N);        
        SOL = cellfun(@(y0) {zeros(pars.chunksize+1, length(y0))}, Y0);
    else
        spikes = zeros(pars.T/pars.dt, pars.N); 
        SOL = cellfun(@(y0) {zeros(length(t0), length(y0))}, Y0);
    end
    for i = 1:pars.N
        SOL{i}(1,:) = Y0{i};
    end
   
    [c_ee, c_ei, c_ie, c_ii] = setup_conns(pars.Ne, pars.Ni, ...
                                           pars.K, pars.ps);

    clsz_e = pars.Ne/pars.K;
    clsz_i = pars.Ni/pars.K;
                                       
    Ge = [125 * c_ee, 60 * c_ei] / clsz_e;
    Gi = [90 * c_ie,  20 * c_ii] / clsz_i;
    ampa_ratio = pars.ampa_ratio;
    nmda_ratio = 1 - ampa_ratio;

    %%
    vr = [repmat(pars.params_e.mem.c, pars.Ne, 1); 
          repmat(pars.params_i.mem.c, pars.Ni, 1)];
    ur = [repmat(pars.params_e.mem.d, pars.Ne, 1); 
          repmat(pars.params_i.mem.d, pars.Ni, 1)];
      
    idx_t = 1;
    for i = 1:length(t0)-1

        if mod(i, 100) == 0
            disp([num2str(i), '/', num2str(pars.T/pars.dt)])
        end
        % update membrane potentials
        tic
        t_1 = t0(i);
        t_2 = t0(i+1);
        for j = 1:pars.N
            [T, Y, TE, YE] = ode23(neur(j).odefun,[t_1,t_2],Y0{j},opts{j});
            if ~isempty(TE)
                spikes(idx_t,j) = 1;
                Y = YE(end,:)';
                Y(1) = vr(j);
                Y(6) = Y(6) + ur(j);
                neur(j).set_spike(TE);
                [T, Y] = ode23(neur(j).odefun, [TE, t_2], Y, opts{j});
            end
            SOL{j}(idx_t+1,:) = Y(end,:);
        end
        toc
        % propagate effect of spikes
        
        s = spikes(idx_t,:);
        
        i_e = s(1:pars.Ne) * Ge;
        i_i = s(pars.Ne+1:end) * Gi;
        for j = 1:pars.N
            % update for next iteration
            Y0{j} = SOL{j}(idx_t+1,:);
            Y0{j}(3) = Y0{j}(3) + i_i(j);
            if j <= pars.Ne
                Y0{j}(2) = Y0{j}(2) + ampa_ratio * i_e(j);
            else
                Y0{j}(2) = Y0{j}(2) + i_e(j);
            end

            if j <= pars.Ne
                spike_e = poissrnd(pars.rate_e);

                Y0{j}(2) = Y0{j}(2) + pars.jXE2 * spike_e;
                Y0{j}(5) = Y0{j}(5) + pars.jXE5 * spike_e;
                spike_i = poissrnd(pars.rate_i);
                Y0{j}(3) = Y0{j}(3) + pars.jXI * spike_i;
                % nmda only on excitatory neurons
                Y0{j}(5) = Y0{j}(5) + nmda_ratio * i_e(j);
            end
            % saturate conductances
            Y0{j}(2:end-2) = min(Y0{j}(2:end-2), [100, 100, 50]);
        end
        
        
        if mod(i,pars.chunksize)==0 && if_save && i > 1
            
            spikes = sparse(spikes);
            
            idx_sf = (i - mod(i, pars.chunksize))/pars.chunksize; % slightly retarted
            save([sf,'_ck',num2str(idx_sf)], 'spikes', 'SOL', 'pars'); 
            
            if i <  length(t0)-1 % i.e. there is still some time left to simulate
                spikes = zeros(pars.chunksize, pars.N);    
                SOL_old = SOL;
                SOL = cellfun(@(y0) {zeros(pars.chunksize, length(y0))}, Y0);
                for j = 1:pars.N
                    SOL{j}(1,:) = SOL_old{j}(idx_t+1,:);
                end            
            else % all data is stored, no left-overs to return
                spikes = [];
                SOL = cell(1, pars.N);
            end
            
            idx_t = 0;
        end
        
        idx_t = idx_t + 1;
    end
    
    %setup.t0 = t0;
    setup.Y0 = Y0;
    setup.Ge = Ge;
    setup.Gi = Gi;
    setup.SOL = SOL;
    %if if_save
    %    spikes = [];
    %end
    
end
