% close all;
close all;
clear all;
load('X_loop_4bit.mat')
optional.offset = 35;
degrees = linspace(0,20);

filter = true;

if filter
    phi_binary_quant = filterMat(phi_binary_quant,-180);
    phi_stochastic_multiOrApproxBinary = filterMat(phi_stochastic_multiOrApproxBinary,-180);
    phi_stochastic_half_wave = filterMat(phi_stochastic_half_wave,-180);
    
end

legendCell = cellstr(num2str((ScaleValues*max_input)', 'Scale*max = %-2.2d'))';
figure();
%% Plots
first = phi_binary_quant; %4times equal result -> not quantized
er_first = generate_Error_Rates(first,degrees,optional);
subplot(131)
p= plot(repmat(degrees',1,size(er_first,2)),er_first(:,:));
xlabel('Error Threshold in degrees')
ylabel('Azimuth Error')
title('Time Domain Double');
legend(legendCell);
ylim([0 1])
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks);
p(1).Marker = 'o';
p(2).Marker = '+';
p(3).Marker = '*';
p(4).Marker = '.';
p(5).Marker = 'x';

second = phi_stochastic_multiOrApproxBinary;
subplot(132);
er_second = generate_Error_Rates(second,degrees,optional);
p = plot(repmat(degrees',1,size(er_second,2)),er_second(:,:));
xlabel('Error Threshold in degrees')
ylabel('Azimuth Error')
title('Mixed: AND + Binary Adder');
legend(legendCell);
ylim([0 1])
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks);
p(1).Marker = 'o';
p(2).Marker = '+';
p(3).Marker = '*';
p(4).Marker = '.';
p(5).Marker = 'x';


second = phi_stochastic_half_wave;
subplot(133);
er_second = generate_Error_Rates(second,degrees,optional);
p=plot(repmat(degrees',1,size(er_second,2)),er_second(:,:));
xlabel('Error Threshold in degrees')
ylabel('Azimuth Error')
title('Full Stochastic');
legend(legendCell);
ylim([0 1])
a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
pct = char(ones(size(a,1),1)*'%'); 
new_yticks = [char(a),pct];
set(gca,'yticklabel',new_yticks);
p(1).Marker = 'o';
p(2).Marker = '+';
p(3).Marker = '*';
p(4).Marker = '.';
p(5).Marker = 'x';

figure();
p= plot(max_sum_result(:,1:3));title('Cross Correlation Maximum of Reference');
legend(legendCell);
ylim([0,10]);
p(1).Marker = 'o';
p(2).Marker = '+';
p(3).Marker = '*';
p(4).Marker = '.';
p(5).Marker = 'x';