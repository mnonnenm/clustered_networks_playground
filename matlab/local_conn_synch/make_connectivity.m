function [Ge, Gi] = make_connectivity(Me, Mi)
    K = 10;
    M = Me + Mi;
    %% Define cells' position
    pos_e = rand(Me, 2);
    pos_i = rand(Mi, 2);
    
    %%
    figure; hold on
    xy = reshape(pos_e, [], 2, K);
    for k = 1:K
        plot(xy(:,1,k), xy(:,2,k),'^');
    end
    plot(pos_i(:,1), pos_i(:,2),'ok');
    pause(0.1);
    
    %% Build connectivity matrix: E -> E
    c_ee = cell(K);
    for i = 1:K
        for j = 1:K
            if i == j
                pc = 0.30;
            else
                pc = 0.01;
            end
            c_ee{i,j} = sparse(rand(Me/K, Me/K) < pc);
        end
    end
    c_ee = cell2mat(c_ee);
    
    %% Build connectivity matrix: I -> I
    c_ii = zeros(Mi);
    for i = 1:Mi-1
        for j = i+1:Mi
            dp = pos_i(i,:) - pos_i(j,:);
            r = sqrt(sum(dp.^2));
            if r < 0.2
                c_ii(i,j) = 1;
                c_ii(j,i) = 1;
            end
        end
    end
    c_ii = sparse(c_ii);
    
    %% Build connectivity matrix: E -> I, I -> E
    dist = 0.2;
    for i = 1:Mi
        dx = pos_e(:,1) - pos_i(i,1);
        dy = pos_e(:,2) - pos_i(i,2);
        r = sqrt(dx.^2 + dy.^2);
        % inhibitory to excitatory
        pc_ie = r < dist;
        c_ie(i,:) = sparse(rand(Me, 1) < pc_ie);
        % excitatory to inhibitory
        pc_ei = r < dist;
        c_ei(:,i) = sparse(rand(Me, 1) < pc_ei);
    end
    
    %% Combine the matrices
    Ge = [K * 125 * c_ee, 250 * c_ei] / M;
    Gi = [450 * c_ie, 250 * c_ii] / M;
    
    %%
    figure
    spy([Ge; Gi])
    
    save(fullfile(pwd,'connectivity'),'Ge','Gi');
end