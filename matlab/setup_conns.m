function [c_ee, c_ei, c_ie, c_ii] = setup_conns(Me, Mi, K, ps)
% setting up connectivity matrices for several E-I modules. 

% Note this function uses a naming convention that is based on left-(!)-
% multiplying state and connectivity: 
% pEI is the connection probability from E to I etc.

% Cross-talk accross modules is only allowed for E-to-E and E-to-I
% connections. The code however is easily extendable to I-to-E and I-to-I.

pEE   = ps.pEE; 
pEEcr = ps.pEEcr;

pEI   = ps.pEI;
pEIcr = ps.pEIcr;

pIE   = ps.pIE;
pII   = ps.pII;

clsz_e = Me/K; % will break below if not integer
clsz_i = Mi/K;


c_ee = false(Me, Me);
for i = 1:K
    for j = 1:K
        if i == j
            c_ee((i-1)*clsz_e+(1:clsz_e),(j-1)*clsz_e+(1:clsz_e)) = ...
                                               rand(clsz_e, clsz_e) < pEE;
        else
            c_ee((i-1)*clsz_e+(1:clsz_e),(j-1)*clsz_e+(1:clsz_e)) = ...
                                             rand(clsz_e, clsz_e) < pEEcr;
        end
    end
end
c_ee = sparse(c_ee);

c_ei = false(Me, Mi);
for i = 1:K
    for j = 1:K
        if i == j
            c_ei((i-1)*clsz_e+(1:clsz_e),(j-1)*clsz_i+(1:clsz_i)) = ...
                                               rand(clsz_e, clsz_i) < pEI;
        else
            c_ei((i-1)*clsz_e+(1:clsz_e),(j-1)*clsz_i+(1:clsz_i)) = ...
                                               rand(clsz_e, clsz_i) < pEIcr;
        end
    end
end
c_ei = sparse(c_ei);

c_ie = false(Mi, Me);
for i = 1:K
    for j = 1:K
        if i == j
            c_ie((i-1)*clsz_i+(1:clsz_i),(j-1)*clsz_e+(1:clsz_e)) = ...
                                               rand(clsz_i, clsz_e) < pIE;
        else
            % we just leave this zero!
        end
    end
end
c_ie = sparse(c_ie);


c_ii = false(Mi, Mi);
for i = 1:K
    for j = 1:K
        if i == j
            c_ii((i-1)*clsz_i+(1:clsz_i),(j-1)*clsz_i+(1:clsz_i)) = ...
                                               rand(clsz_i, clsz_i) < pII;
        else
            % we just leave this zero!
        end
    end
end
c_ii = sparse(c_ii);

end