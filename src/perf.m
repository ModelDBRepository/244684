miss = zeros(n_post,n_pattern);
hit = zeros(n_post,n_pattern);
false_alarm = zeros(1,n_post);

for n=1:n_period_record_spike
    p = mod(n-1,n_pattern)+1;

    %check for false alarms
%     range = round((n-1)*period/dt)+2:round((n*period-pattern_duration)/dt);
%     range_prec = round((n-1)*period/dt)+1:round((n*period-pattern_duration)/dt)-1;
%     [neuron, time ] = find(V_post(:,range)==0 & V_post(:,range-1)~=0); 
    neuron = spike_list( spike_list(:,1)>(n_period-n)*period & spike_list(:,1)<=(n_period-n+1)*period-pattern_duration-2*jitter , 2);
    false_alarm(ismember(1:n_post,neuron)) = false_alarm(ismember(1:n_post,neuron))+1;
    
    %check for misses
%     range = round((n*period-pattern_duration)/dt)+2:round(n*period/dt);
%     range_pred = round((n*period-pattern_duration)/dt)+1:round(n*period/dt)-1;
%     [neuron, time ] = find(V_post(:,range)==0 & V_post(:,range-1)~=0);
    neuron = spike_list( spike_list(:,1)>(n_period-n+1)*period-pattern_duration-2*jitter & spike_list(:,1)<=(n_period-n+1)*period , 2);
    miss(~ismember(1:n_post,neuron),p) = miss(~ismember(1:n_post,neuron),p)+1;
    hit(:,p) = hit(:,p) + histcounts(neuron,1:n_post+1)';
end

n_w = sum(w>.5,2);

n_w_pattern = zeros(n_pattern,n_post);

for p = 1 : n_pattern
    spike_count = sum(pattern{p},2);
    n_w_pattern(p,:) = sum( (w>.5)  .* repmat(spike_count',[n_post 1]) , 2 );
end

file_name = '../data/perf_';
if exist('seed','var')
    file_name = [ file_name sprintf('%03d',seed) ];
else
  c = clock;
  file_name = [ file_name sprintf('%04.0f',c(1)) '-' sprintf('%02.0f',c(2)) '-' sprintf('%02.0f',c(3)) '-' sprintf('%02.0f',c(4)) '-' sprintf('%02.0f',c(5)) ];
end

file_name = [ file_name '.mat' ];

save(file_name,'miss','false_alarm', 'hit', 'n_w', 'n_w_pattern','n_thr','n_dw_post','n_period_record_spike','n_pattern')
