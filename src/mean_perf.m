% Computes the mean performance of several simulations (with different
% random seeds), for example produced using batch.py


n_involved = 1e4;
f = 3.2;
n_pattern = 40;
pattern_duration = 100e-3;
epsilon = .05;

path = '../data/';
%path = ['../data_' sprintf('%02d',n_pattern) '_pattern_paper/'];


optim_delta_t = [1 23e-3; 2 17e-3 ; 5 11e-3; 10 8.1e-3; 20 5.7e-3; 40  3.7e-3 ];
delta_t = optim_delta_t(optim_delta_t(:,1)==n_pattern,2); % expected duration of the reinforced subsequence



M = 0;
for i=1:n_pattern
    M = M + nchoosek(n_pattern,i)*(-1)^(i+1)*(1-exp(-f*delta_t))^i;
end
M = M * n_involved;


%optimal_n_w = 1e4 * (1-exp(-11e-3*3.2));
optimal_n_w = n_involved*f*delta_t+M*f*(pattern_duration-delta_t);

% optimal_n_w = 430;

list = dir([path 'perf_*.mat']);
disp(['Displaying mean of ' int2str(length(list)) ' perf files' ])
for l=1:length(list)
    load([ path list(l).name ])
    if l==1
        miss_ = miss;
        false_alarm_ = false_alarm;
        hit_ = ( hit>=90 & hit<=110 );
        %hit_ = hit;
        w_ = abs((n_w-M)/M)<epsilon;
        w_pattern = sum(abs((n_w_pattern-optimal_n_w)/optimal_n_w)<epsilon ,1 );
        %w_pattern = mean(n_w_pattern , 1 );
        %w_ = n_w;
        perfect = abs((n_w-M)/M)<epsilon & sum(hit>=1,2)==n_pattern & false_alarm'==0;
        
        n_completely_missed =  sum(miss==100,2);
        
        n_completely_missed_ = n_completely_missed;
        hit_learned = sum(100-miss,2) ./ (n_pattern-n_completely_missed);
        
        
    else
        miss_ = miss_ + miss;
        false_alarm_ = false_alarm_ + false_alarm;
        hit_ = hit_ + (hit>=90 & hit<=110);
        %hit_ = hit_ + hit;
        w_ = w_ + ( abs((n_w-M)/M)<epsilon );
        w_pattern = w_pattern + sum(abs((n_w_pattern-optimal_n_w)/optimal_n_w)<epsilon ,1 );
        %w_pattern = w_pattern + mean(n_w_pattern , 1 );
        %w_ = w_ + n_w;
        
        perfect = perfect + ( abs((n_w-M)/M)<epsilon & sum(hit>=1,2)==n_pattern & false_alarm'==0 );

        n_completely_missed =  sum(miss==100,2);
        
        n_completely_missed_ = n_completely_missed_ + n_completely_missed;
        hit_learned = hit_learned + sum(100-miss,2) ./ (n_pattern-n_completely_missed);
        
    end
    if l<=0
        figure('Name',list(l).name)
        imagesc(reshape(1-( (hit>=95 & hit<=105) ),n_thr,n_dw_post),[0 1])
    end
end
miss_ = miss_/length(list);
false_alarm_ = false_alarm_/length(list);
hit_ = hit_/length(list);
w_ = w_/length(list);
w_pattern = w_pattern/length(list);
perfect = perfect/length(list);

n_completely_missed_ = n_completely_missed_/length(list);
hit_learned = hit_learned/length(list);

max_miss = n_period_record_spike/n_pattern;
max_fa = n_period_record_spike;

figure('Position',[10 10 60 25]*20,'Units','centimeters')

subplot(2,4,1)
imagesc(reshape(false_alarm_,n_thr,n_dw_post),[0 max_fa])
colorbar
ylabel('thr')
xlabel('LTD')
title('False alarms')

subplot(2,4,2)
imagesc(reshape(sum(miss_,2),n_thr,n_dw_post))
colorbar
ylabel('thr')
xlabel('LTD')
title('Misses')

% subplot(2,2,3)
% imagesc(reshape(min(max_miss,sum(miss_,2))/max_miss+min(max_fa,false_alarm_')/max_fa,n_thr,n_dw_post),[0 1])
% colorbar
% ylabel('thr')
% xlabel('LTD')
% title('Combined score')

subplot(2,4,3)
%imagesc(reshape(1-sum(hit_,2),n_thr,n_dw_post),[0 1])
imagesc(reshape(sum(hit_,2),n_thr,n_dw_post),[0 n_pattern])
colorbar
ylabel('thr')
xlabel('LTD')
title('Hits')


subplot(2,4,4)
imagesc(reshape(w_,n_thr,n_dw_post),[0 1])
%imagesc(reshape(w_,n_thr,n_dw_post),[300 400])
colorbar
ylabel('thr')
xlabel('LTD')
title('Optimal n w')

subplot(2,4,5)
imagesc(reshape(w_pattern,n_thr,n_dw_post),[0 n_pattern])
%imagesc(reshape(w_,n_thr,n_dw_post),[300 400])
colorbar
ylabel('thr')
xlabel('LTD')
title('Optimal n w pattern')


subplot(2,4,6)
imagesc(reshape(perfect,n_thr,n_dw_post),[0 1])
colorbar
ylabel('thr')
xlabel('LTD')
title('Perfect')

subplot(2,4,7)
imagesc(reshape(n_completely_missed_,n_thr,n_dw_post),[0 n_pattern])
colorbar
ylabel('thr')
xlabel('LTD')
title('n completely missed')

subplot(2,4,8)
imagesc(reshape(hit_learned,n_thr,n_dw_post))
colorbar
ylabel('thr')
xlabel('LTD')
title('HR for learned pat.')


colormap jet
