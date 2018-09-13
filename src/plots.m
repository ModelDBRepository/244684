%__________________________________________
% PLOTTING
T = min(20*n_pattern*period,n_period_record*period)-dt;
%T = (n_period_record-1)*period;
n = round(T/dt)-1;

if ~exist('neuron_subset','var')
    neuron_subset = 1:n_post;
    % neuron_subset = 3*n_thr+1;
end

stop = size(V_post,2);

figure
co = brighten(get(gcf,'DefaultAxesColorOrder'),.8);

%subplot(2,1,1)
last_pattern_idx = mod(n_period-1,n_pattern)+1;
for i=1:round(T/period)
    rectangle('Position',[(n_period_record+1-i)*period-pattern_duration-jitter,0,pattern_duration,1.5*max(thr)],'FaceColor',co(mod(last_pattern_idx-i,min(n_pattern,7))+1,:),'EdgeColor','none')
    hold on
end
plot(dt*(max(stop-n,1):stop),V_post(neuron_subset,max(stop-n,1):stop)')
xlabel('t (s)')
ylabel('Potential')
xlim(dt*[max(stop-n,1) stop])

% legend([ int2str(thr(neuron_subset)) repmat('-',length(thr(neuron_subset)),1) num2str(dw_post(neuron_subset)/da_pre,2) ])

% % recompute adaptive threshold
% thr_rec = zeros(size(V_post));
% thr_rec(:,1)=thr;
% for i=1:size(thr_rec,2)-1;
%     for j=1:size(thr_rec,1)
%         thr_rec(j,i+1) = thr(j) + (thr_rec(j,i)-thr(j)) * (1-dt/tau_thr);
%         if V_post(j,i+1) == 0
%             thr_rec(j,i+1) = thr_rec(j,i+1) + d_thr * thr(j);
%         end
%     end
% end
% plot(dt*(max(stop-n,1):stop),thr_rec(neuron_subset,max(stop-n,1):stop)',':')



% subplot(2,1,2)
% for i=1:n_post
%     spike = dt*(max(stop-n,1)-1+find(V_post(i,max(stop-n,1):stop)==0 & V_post(i,(max(stop-n,1):stop)-1)));
%     plot(spike,i*ones(size(spike)),'.')
%     hold on
% end
% xlim(dt*[max(stop-n,1) stop])



figure
hist(w(neuron_subset,:)',20)
xlim([0 1])
ylabel('#')
xlabel('Normalized synaptic weights')

%return

n_pattern_disp = min(n_pattern,5);
if length(neuron_subset)<=3
    for i=1:length(neuron_subset)
        figure('Name',['Neuron #' int2str(neuron_subset(i))],'Position',[1 1 n_pattern_disp+2 4]*200,'Units','centimeters')
        for p=1:n_pattern_disp
            subplot(1,n_pattern_disp,p)
            %im = pattern{p} .* repmat(w',1,size(pattern{p},2));
            %imagesc(im);
            for a=1:n_pre
                spike_time = find(pattern{p}(a,:))*dt;
                plot(spike_time,a*ones(size(spike_time)),'.','Color',w(neuron_subset(i),a)*[1 0 0]+(1-w(neuron_subset(i),a))*[0 0 1],'MarkerSize',6)
                hold on
            end
            
        end
    end
end
