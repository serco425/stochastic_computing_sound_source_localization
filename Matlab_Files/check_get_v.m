% Thi script compares the new vs old formulas

addpath('./stochastic_computing') % add getShifts function

get_v = @(N,n)  n/((sqrt(4*N+1)-1)/2+1);

n_max = 256;
v_max = 256;


for n = 4:2:n_max

    for v= 1:v_max

        shift_vec = getShifts(n,v);
        N= length(shift_vec);


        %compare with derived formular

        v_calculate = floor(get_v(N,n)); %should return the maximum v

        if v~=v_calculate

            %calculates shift vector with new calculated v and check if
            %addends stay the same.

            N_compare = length(getShifts(n,v_calculate));

            if N_compare ~= N
                % problem
                disp('Not matching');
            end


        end


    end
end


