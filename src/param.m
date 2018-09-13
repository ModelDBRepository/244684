% This code was used in: Masquelier & Kheradpisheh (2018) Optimal localist and distributed coding of spatiotemporal spike patterns through STDP and coincidence detection. Frontiers in Computational Neuroscience.
% with Matlab R2016b
% Aug 2018
% timothee.masquelier@cnrs.fr
%__________________________________________
% PARAMETERS

% Note: a period contains some Poisson activity followed by the
% presentation of (one of) the pattern(s).
%

% stimulation
f = 3.2; % mean rate in Hz
n_pattern = 5; % number of patterns
n_period_record = 5; % number of periods to record potential (for plotting)
n_period_record_spike = 100*n_pattern; % number of periods to estimate performance (hit rate and false alarm rate, ...)
n_period = 12000; % total number of periods
period = 4*2*50e-3; % period duration
pattern_duration = 2*50e-3; % pattern duration
n_pre = 1e4; % nb of presynaptic neurons
n_involved = n_pre; % number of afferents involved in the pattern(s)
optim_delta_t = [ 1 23e-3; 2 17e-3 ; 5 11e-3; 10 8.1e-3; 20 5.7e-3; 40  3.7e-3 ]; % see Table 1 in the paper
delta_t = optim_delta_t(optim_delta_t(:,1)==n_pattern,2); % expected duration of the reinforced subsequence
jitter = 3.2e-3; % each pattern spike is shifted by a random lag uniformly distributed in [-jitter,jitter]

n_thr = 1; % Nb of different threshold values (for parameter search)
n_dw_post = 1; % Nb of different dw_post values (for parameter search)
n_post = n_thr*n_dw_post; % nb of postynaptic neurons

% Spike timing-dependent plasticity (STDP). See Song, Miller & Abbott 2000 Nat Neurosc
tau_pre = 20e-3; % LTP time constant
% tau_post = 10e-3; % LTD time constant

da_pre = 0.1; % LTP variable increase (at each presynatpic spike)

% % Uncomment this for increasing learning rate:
% n = n_period/100;
% ratio = 20;
% alpha = ratio^(1/(n-1));
% mean_coeff = 1/n*(1-alpha^n)/(1-alpha);
% da_pre = da_pre / mean_coeff;
% disp(['da_pre initial ' num2str(da_pre)])
% disp(['da_pre final ' num2str(da_pre*alpha^(n-1))])
% da_min = 0*da_pre;
% alpha = 2/(n+1)*(da_pre-da_min);
% da_pre = da_min+alpha;
% disp(['da_pre initial ' num2str(da_pre)])
% disp(['da_pre final ' num2str(da_pre+alpha*(n-1))])

%da_post = 0*1.0*exp(log(1.1)*floor((0:n_post-1)/n_thr))'*da_pre; % LTD variable increase (at each postsynaptic spike)
da_post = zeros(n_post,1); % LTD variable increase (at each postsynaptic spike)

% Homeostatic (at each postsynaptic spike, all the synaptic weights are decreased by this value. Use an array of different values for parameter seach)
dw_post = 0.062*da_pre*exp(log(1.025)*floor((0:n_thr*n_dw_post-1)/n_thr-(n_dw_post-1)/2))'; % see Table 1 in the paper

% Neurons
optim_tau  = [ 1 18e-3 ; 2 13e-3 ; 5 8.9e-3 ; 10 6.8e-3 ; 20 5.6e-3 ; 40 5.1e-3 ]; % see Table 1 in the paper
tau_m = optim_tau(optim_tau(:,1)==n_pattern,2); % Membrane time constant
tau_s = 0*5e-3; % Synapse time constant. Put 0 for instantaneous synapse.
tau_thr = 80e-3; % Adaptive threshold time constant
d_thr = 1.8; % Adaptive threshold increase at each postsynaptic spikes (divided by the resting threshold)
dt = .1e-3; % Time step
if jitter==0
    v_max = (1-exp(-delta_t/tau_m));
else
    v_max = min(1,delta_t/(2*jitter))-tau_m/(2*jitter)*log(1-exp(-max(2*jitter,delta_t)/tau_m)+exp(-abs(2*jitter-delta_t)/tau_m));
end

% Expected number of selected afferents
M = n_involved * (1 - exp(-n_pattern * f * delta_t));

% Threshold. Use an array of different values for parameter seach.
% thr =  tau_m*f*(M+v_max*(n_involved-M))*repmat( exp(log(1.025)*(-(n_thr-1):0)'),n_dw_post,1); 
thr =  178*repmat( exp(log(1.025)*(-(n_thr-1):0)'),n_dw_post,1); % see Table 1 in the paper

initial_distance_to_threshold = -1; % distance between E(V) and thr divided by std(V) (negative <=> mean driven regime; positive <=> fluctuation driven regime)

disp(['E(n_selected)=' num2str(M) ])
disp(['E(V_noise)=' num2str(tau_m*M*f) ])
disp(['thr=' num2str(thr(round((length(thr)+1)/2)))])
%disp(['thr=' num2str(thr(end))])


% parameters for the pre-calculated poisson input spikes
poisson_spike.m = n_pre;
poisson_spike.n = 2e6;
poisson_spike.p = f*dt;
poisson_spike.cursor = 1;


% % Rapidly adapting thr (Fontaine et al PCB 2014)
% tau_thr = Inf; % put Inf not to use
% V_thr = 0;
% if tau_thr < Inf
%     exp_filter = exp(-(0:dt:5*tau_thr)/tau_thr);
%     epsp = exp(-(0:dt:5*tau_m)/tau_m)-exp(-(0:dt:5*tau_m)/tau_s);
%     epsp = epsp/max(epsp);
%     effective_epsp = epsp-filter(exp_filter,sum(exp_filter),epsp);
%     coef = mean(effective_epsp)/mean(epsp);
%     coef = .019;
% else
%     coef = 1;
% end
