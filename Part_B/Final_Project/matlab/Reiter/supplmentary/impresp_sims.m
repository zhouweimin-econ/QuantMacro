%% Computes impulse responses for solution of linear RE models,
% assumes i.i.d. shock process (z)
% Inputs:
%   G1 and Impact: Output from Sims' package
%   HORIZON:       number of periods for which to compute IR
%   scale_shocks:  scaling factor of shocks
%   scale_y:       scaling factor for variables (i.e. for annualized variables, etc.)

% Ouput: Cell array of impulse responses, one element for each shock

function Resp_mat = impresp_sims(G1,Impact,HORIZON,scale_shocks,scale_y)

nz = size(Impact,2);
NN = zeros(nz,nz);
[m_states,k_exog] = size(Impact);
Resp_mat = zeros(m_states+k_exog,HORIZON,k_exog);

for shock_counter = 1 : k_exog
    Response = zeros(m_states+k_exog,HORIZON);
    iVar = m_states+shock_counter;
    Response(iVar,1) = 1;
    II_lag = [G1,zeros(m_states,k_exog)
        zeros(k_exog,(m_states)), NN];
    aux1 = [ zeros(m_states,m_states), Impact];
    aux3 = [ zeros(k_exog,  m_states), zeros(k_exog,k_exog) ];
    II_contemp = eye(m_states+k_exog) + [aux1;aux3];

    Response(:,1) = II_contemp*Response(:,1);
    for time_counter = 2 : HORIZON
        Response(:,time_counter) = II_contemp*II_lag*Response(:,time_counter-1);
    end
    Response=Response.*repmat(scale_shocks(shock_counter)*[scale_y;zeros(k_exog,1)],1,HORIZON);

    Resp_mat(:,:,shock_counter) =  Response;
end
    
    
