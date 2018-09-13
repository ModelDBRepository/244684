function jittered_pattern = jitter_pattern(pattern,jitter,f,dt)

n = size(pattern,1);
m = 2*round(2*jitter/dt);

n_spike = poissrnd(m*n*f*dt);            
position = randperm(m*n,n_spike);
i_ = mod(position-1,n)+1;
j_ = floor((position-1)/n)+1;
poisson_spike = sparse(i_,j_,ones(1,n_spike),n,m);
poisson_spike = poisson_spike>0;

extended_pattern = [ poisson_spike(:,1:m/2) pattern poisson_spike(:,m/2+1:end) ];
%extended_pattern = [ sparse( rand(n,2*jitter/dt)<dt*f ) pattern sparse( rand(n,2*jitter/dt)<dt*f ) ];
m = size(extended_pattern,2);
jittered_pattern = logical(sparse(n,m));

for s = find(extended_pattern)'
    [current_line, current_col] = ind2sub([n,m],s);
    new_col = current_col + round(2*(rand-.5)*jitter/dt);
    if new_col>0 && new_col<=m
%         if jittered_pattern(current_line,new_col)
%             warning('Overwriting a spike')
%         end
        jittered_pattern(current_line,new_col) = true;
    end
end

jittered_pattern = jittered_pattern(:,jitter/dt+1:end-jitter/dt);
