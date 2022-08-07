close all;clear all;
% The test in this file has nothing to do with direction estimation
% Tests the amount of overlap when c additional bits are true


for n = 16%3:32
    Prime_select = [n,n-1]; %relative Prime

    for v = 4%1:min(floor(n/2),n-2)  

    k = Prime_select(2);
    nk = Prime_select(1)*Prime_select(2);    

    N_Periode = (nk+(1*v*n-k));

    c_max = min([floor(N_Periode/(2)-v),  v, k-v-1, floor((k^2+4*k)/(3*k+2))]); %bits required  letzte min ist wegen länge von I_minor
    c_additional_bits_vec = 0:c_max;

    x_sc = zeros(length(c_additional_bits_vec),Prime_select(1));
    y_sc = zeros(length(c_additional_bits_vec),Prime_select(2));
    mul = zeros(length(c_additional_bits_vec),nk);

        for ind = 1:length(c_additional_bits_vec)
            c = c_additional_bits_vec(ind);

            x_sc(ind,:) = pwmVec((v+c)/Prime_select(1),Prime_select(1));
            y_sc(ind,:) = pwmVec((v+c)/Prime_select(2),Prime_select(2));
            mul(ind,:) = scAndMul(x_sc(ind,:),y_sc(ind,:),nk);
        end
        

        
    end
    
            %check increase of c in an interval
        %I'_F,major
        % versuch es mit cells generischer zu machen.
        I_tick_F_major = n*v-k:k*(k-v+2)-1;
        I_begin = 0:n*v-k-1;
        I_end = k*(k-v+2):nk;
        
        x_tick_intervals_cell = {};
        x_tick_cell = {};
        
        x_tick_cell(1) = {0:v-2};
        x_tick_cell(2) = {v-1};
        x_tick_cell(3) = {v:v+c-1};
        x_tick_cell(4) = {v+c:k-(v+c)};
        x_tick_cell(5) = {k-(v+c)+1:k-v};
        x_tick_cell(6) = {k-v+1};
        x_tick_cell(7) = {k-v+2:k-1};
        
        x_tick_0 = 0:v-2;%x'_0
        x_tick_1 = v-1;
%         I_x_tick_0 = x_tick_0(1)*n:x_tick_0(end)*k+(v+c)-1;
%         I_x_tick_1 = x_tick_1(1)*n:x_tick_1(end)*k+(v+c)-1;
        
        for ind = 1:7
            cur_x_tick = cell2mat(x_tick_cell(ind));

            if ~isempty(cur_x_tick)
                if ind <= 3
                    x_tick_intervals_cell(ind) = {cur_x_tick(1)*n:cur_x_tick(end)*k+(v+c)-1};
                else
                    x_tick_intervals_cell(ind) = {(cur_x_tick(1)+1)*k : n*cur_x_tick(end)+(v+c)-1};
                end
            else
                x_tick_intervals_cell(ind) = {cur_x_tick};
            end
        end
        
        %x'_0 -> no overlap?    even with highest c
        for interval_ind = 1:7
            if max(cell2mat(x_tick_intervals_cell(interval_ind)))>=min(I_tick_F_major)
                disp(['Interval ', num2str(interval_ind-1), ' overlap']);
            end
        end
%       why overlap in interval 4?
        
        
        %x'1 -> some overlap depending on c
        sum_x_tick = zeros(7,length(c_additional_bits_vec));
        for interval_ind = 1:7
            cur_x_interval = cell2mat(x_tick_intervals_cell(interval_ind));
            for c_ind = 1:length(c_additional_bits_vec)
                tmp = mul(c_ind,:);
                tmp(1+I_begin)=0;   % Only bits in I_F,major are problematic
                tmp(1+I_end) = 0;   
                sum_x_tick(interval_ind,c_ind) = sum(tmp(1+cur_x_interval)); % sums bits in I_F,major  
            end
        end    
end