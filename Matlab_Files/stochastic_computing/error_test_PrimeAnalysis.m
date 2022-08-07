close all;
% The test in this file has nothing to do with direction estimation
% It just shows the error of stochastic computation

% How primes should be chosen to best add 256 Values
%% Configure
N = 256; %length of Vector (xcorrelation size)
N_repeat_for = 1000; % for averaging out the result 10000
Magnitudes_output = linspace(0.01,0.6,20); %linspace(0.1,0.2,N_output_range);
N_output_range = length(Magnitudes_output);

% Prime_mat = [];

Prime_mat = [[16,15];[16,13];[16,11];[16,9]];
N_Primes = size(Prime_mat,1);


%% Storage Variables
%absolute Values stored here
sum_reference = zeros(1,N_repeat_for);
sum_pwm_rotation_or = zeros(1,N_repeat_for);
sum_sobol = zeros(1,N_repeat_for);

dot_magnitude = zeros(1,N_repeat_for);
mean_dot_magnitude = zeros(N_Primes,N_output_range);

mean_reference = zeros(N_Primes,N_output_range);
mean_pwm_rotation_or = zeros(N_Primes,N_output_range);
mean_sobol = zeros(1,N_output_range);

%relative Values stored here
mean_rel_pwm_rotation_or = zeros(N_Primes,N_output_range);
mean_rel_sobol = zeros(N_Primes,N_output_range);

for prim_index = 1:N_Primes 
    Prime_select = Prime_mat(prim_index,:);


%% Calculation Loop
for ind_magnitude = 1:N_output_range

target_magnitude = Magnitudes_output(ind_magnitude);

%LSB
required_scale = sqrt(target_magnitude/(N*(0.5)^2));

for ind_repeat = 1:N_repeat_for
x = rand(1,N)*required_scale-required_scale/2;
y = rand(1,N)*required_scale-required_scale/2;

%processing and filtering of random numbers so it behaves more like an audio signal
%in ciurcuit half of the signal will be zero due to uniform processing
x = lowpass(x,6000,16000);
y = lowpass(y,6000,16000);
x(x<0) = 0;
y(y<0) = 0;
x = x*2;
y = y*2;

% calculate
mul_len = Prime_select(1)*Prime_select(2);
x_sc = zeros(N,Prime_select(1));
y_sc = zeros(N,Prime_select(2));
mul = zeros(N,mul_len);
x_sobol = zeros(N,mul_len);
y_sobol = zeros(N,mul_len);
mul_sobol = zeros(N,mul_len);
optional.curPrime = Prime_select;

%% SOBOL
p = sobolset(N*2);
sobol_rand = net(p,mul_len);

for ind = 1:N
    x_sobol(ind,:) = sobolVec(x(ind),sobol_rand(:,2*ind));
    y_sobol(ind,:) = sobolVec(x(ind),sobol_rand(:,2*ind-1));
    mul_sobol(ind,:) = scAndMul(x_sobol(ind,:),y_sobol(ind,:),mul_len);
end
optional.type = 'none';
sum_sc = scMultiOrAdd(mul_sobol,mul_len);
sum_sobol(ind_repeat) = Unary2Binary(sum_sc);


%% PWM
for ind = 1:N
    x_sc(ind,:) = pwmVec(x(ind),Prime_select(1));
    y_sc(ind,:) = pwmVec(x(ind),Prime_select(2));
    mul(ind,:) = scAndMul(x_sc(ind,:),y_sc(ind,:),mul_len);
end

% figure();plot(sum(mul_sobol,2)-sum(mul,2))

% do it with Shift
optional.type = 'RotateSum';
sum_sc = scMultiOrAdd(mul,mul_len,optional);
sum_pwm_rotation_or(ind_repeat) = Unary2Binary(sum_sc);

%% Refernces
% simple addition
mul_binary = Unary2Binary(mul); 
sum_reference(ind_repeat) = sum(mul_binary)+eps; %eps avoid division by zero

%Rerference
dot_magnitude(ind_repeat) = dot(x,y); 
end

mean_rel_pwm_rotation_or(prim_index,ind_magnitude) = mean(abs(sum_reference-sum_pwm_rotation_or)./sum_reference);
mean_rel_sobol(prim_index,ind_magnitude) = mean(abs(sum_reference-sum_sobol)./sum_reference);

mean_dot_magnitude(prim_index,ind_magnitude) = mean(dot_magnitude);

mean_reference(prim_index,ind_magnitude) = 		mean(sum_reference)		;
mean_pwm_rotation_or(prim_index,ind_magnitude) = 	mean(sum_pwm_rotation_or)	;
mean_sobol(prim_index, ind_magnitude) = mean(sum_sobol) ;

end
end

%% Plot Section

absolute_Add = [abs(mean_reference'-mean_pwm_rotation_or')];
relative_Add = [mean_rel_pwm_rotation_or]';
magnitudes = {num2str(mean_dot_magnitude,3)};
magnitudes_mat = mean_dot_magnitude';

% Create Legend


% load('remote_magnitude.mat')
% 
legendCell = cellstr(num2str(Prime_mat, '%-d & %-d'))';

% only comparison to Add=a+b is required
figure();
plot(magnitudes_mat,relative_Add);
xlabel('Magnitude of Result')
ylabel('Error')
title('Relative Error |x-x`|/x`');
legend(legendCell);

figure();
p = plot(magnitudes_mat,absolute_Add);
xlabel('Magnitude of Result')
ylabel('Error')
title('Absolute Error |x-x`|');
legend(legendCell);
% 
% p(8).Linewidth = 3;

%% Make plots more distuinguishable
% p(1).Marker = 'o';
% p(2).Marker = '+';
% p(3).Marker = '*';
% p(4).Marker = '.';
% p(5).Marker = 'x';
% p(6).Marker = 's';
% p(7).Marker = 'd';
% p(8).Marker = '^';
% p(9).Marker = 'v';
% p(10).Marker = '>';
% p(11).Marker = '<';
% p(12).Marker = 'p';
% p(13).Marker = 'h';