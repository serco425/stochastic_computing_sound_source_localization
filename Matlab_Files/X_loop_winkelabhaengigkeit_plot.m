% This file takes the location estimation for different azimuth angels
% and wants to compare the azimuth dependent effect of scaleing
% it should show the best magnitude to feed to the pwm module.
% input is 16 and 15 bit

% >> load('X_loop_winkel_extensive.mat','phi_binary_quant','phi_stochastic_multiOrApproxBinary','phi_stochastic_half_wave','ScaleValues','AzimuthValues','max_input')

% Variables
close all;
clear all;
load('X_loop_winkel_select_remote.mat')
ideal = AzimuthValues;

% Format of Input Vectors: % Scale x blocks x Azimuth

filter = true;
plotOnly = 2:5; % allows to remove some curves from plotting
ScaleValues = ScaleValues(plotOnly);

if filter    %all zeroy only values are set to -180 >> this allows filtering those and reducing the number of total_blocks
    phi_binary_quant = filterMat3D(real(phi_binary_quant(:,plotOnly,:)),-180);
    phi_stochastic_multiOrApproxBinary = filterMat3D(real(phi_stochastic_multiOrApproxBinary(:,plotOnly,:)),-180);
    phi_stochastic_half_wave = filterMat3D(real(phi_stochastic_half_wave(:,plotOnly,:)),-180);
end

ideal_rep = repmat(ideal,length(ScaleValues),1); %might change total_blocks to size(phi,1)
legendCell = cellstr(num2str((ScaleValues*max_input)', 'Scale*max(Input) = %-2.2d'))';

for i = 1:3
    
    if i==1
       if filter
           data = phi_binary_quant;
       else
           data = phi_binary_quant(:,plotOnly,:);
       end
       title_str = ['Reference with quantization'];
    elseif i==2
       if filter
           data = phi_stochastic_multiOrApproxBinary;
       else
           data = phi_stochastic_multiOrApproxBinary(:,plotOnly,:);
       end
       title_str = ['Stochastic Products, Binary Summing'];
    elseif i ==3
       if filter
           data = phi_stochastic_half_wave;
       else
           data = phi_stochastic_half_wave(:,plotOnly,:);
       end
       title_str = ['Stochastic Products, Stochastic Summing'];
    end
data_permute = permute(data,[2,1,3]);
% data_plot = squeeze(var(data_permute,0,2)); %calculate Variance
% ylabel_str = 'Variance of Estimation in degree';

data_plot = squeeze(mean(data_permute,2)); %calculate Variance
ylabel_str = 'Estimation in degree';

figure()
plot(ideal_rep',data_plot')
legend(legendCell,'Interpreter','latex');
title(title_str);
xlabel('Azimuth of Source in degree')
ylabel(ylabel_str)
end