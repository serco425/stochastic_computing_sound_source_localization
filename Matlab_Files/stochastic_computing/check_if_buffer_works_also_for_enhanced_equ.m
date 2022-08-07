% This script tests the idea of using a cyclic implementation not only for
% signle delays but also for enhanced equation (with intermediate shifts)

close all; clear all;

%% Configure
%for n=16 v=5 we should be able to add 2

n_simulate_vector = 32;%8:16;

%N_Major_Shifts = 1;
%N_Minor_Shifts = N_Major_Shifts+1;
%summands = (N_Major_Shifts+1)*(N_Minor_Shifts+1);
additional_buffer_length = 0;%floor(-30):floor(30);
overlap_results = zeros(1,length(additional_buffer_length));

for n_ind=1:length(n_simulate_vector)

n = n_simulate_vector(n_ind);
Prime_select = [n,n-1]; %relative Prime
k = Prime_select(2);
%n = Prime_select(1);
nk = Prime_select(1)*Prime_select(2);

for v = 10%1:floor(n/2)

%v_with_one_major = floor(n/(ceil((sqrt(4*6+1)-1)/2)+1));
%v_with_two_major = floor(n/(ceil((sqrt(4*12+1)-1)/2)+1));
%v = v_with_one_major;

N_Major_Shifts = floor((n-2*v)/(v));
N_Minor_Shifts = N_Major_Shifts+1;
N = (N_Major_Shifts+1)*(N_Minor_Shifts+1);

x_sc = zeros(N,Prime_select(1));
y_sc = zeros(N,Prime_select(2));
mul = zeros(N,nk);

for ind = 1:N
    x_sc(ind,:) = pwmVec((v)/Prime_select(1),Prime_select(1));
    y_sc(ind,:) = pwmVec((v)/Prime_select(2),Prime_select(2));
    mul(ind,:) = scAndMul(x_sc(ind,:),y_sc(ind,:),nk);
end

shift_vec = getShifts(n,v);


mul_shifted = scDoRelativeShift(mul,shift_vec);
sc_result = scMultiOrAdd(mul_shifted,size(mul_shifted,2));
%quant(sum(sc_result)/(n*(n-1)),1000*eps) == quant((N)*(v/(n-1)*v/(n)),1000*eps)
if find((sum(mul_shifted,1))>1)
    disp("Wrong Simulation Contraints");
end


%Start Buffer Implementation
M_buf = nk; % Base Buffer Length
S_Buf = N_Major_Shifts*v*n+N_Minor_Shifts*v;

buffer_size = S_Buf+M_buf;


compensate_shifts = getShifts(n,v);
compensate_prev_shifts = [0, compensate_shifts(1:end-1)];
add_delay = [0 repmat(S_Buf,1,length(compensate_shifts)-1)]+getShifts(n,v)-compensate_prev_shifts;


mul_queue = [];

for ind = 1:(size(mul,1))
    lengths(ind) = length(mul_queue);
    restarts(ind) = length(mul_queue)+length(zeros(1,add_delay(ind)));
    mul_queue = [mul_queue, zeros(1,add_delay(ind)), mul(ind,:)];
    
end


row = 1;
col = 1;
mul_stacked = zeros(ceil(length(mul_queue)/buffer_size),buffer_size);%zeros(size(mul_shifted,1),size(mul_shifted,2))*10000;

for ind = 1:length(mul_queue)  
    mul_stacked(row,col) = mul_queue(ind);
    if col == buffer_size
        %check = sum(mul_stacked,1);
        %stairs(check(1,:));
        %ylim([-0.1,1.1]);
        row = row+1;
        col = mod(col,buffer_size);
    end
    
    if mod(ind,n*k)==0
        %figure();
        %check = sum(mul_stacked,1);
        %stairs((check(1,:)));
        %ylim([-0.1,1.1]);
        a=3;
    end
    col = col+1;
end

check = sum(mul_stacked,1);
figure();plot(check);
ylim([-0.1,1.1]);
        
if length(find(check>1))>0
    disp("Overlapp occured");
        X = ["N ", string(N), " N_Major_Shifts " , string(N_Major_Shifts), ...
        " N_Minor_Shifts ", string(N_Minor_Shifts), " n ", string(n)];
    disp(X)
%    overlap_results(modify_buffer_ind) = 1+overlap_results(modify_buffer_ind);
end

%end
close all;
end
end

%figure();plot(overlap_results)

%for disp_ind = find(overlap_results==0)
%    X = ["N ", string(N), " N_major " , string(Large_Shifts), ...
%        " N_minor ", string(Small_Shifts), " Factor ", additional_buffer_length(disp_ind)...
%        " Bufs" ];
%    disp(X)
%end
%bis jetzt - bester favourite: N+1 = 6 smallshifts