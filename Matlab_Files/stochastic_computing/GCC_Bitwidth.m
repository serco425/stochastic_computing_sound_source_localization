close all;
load('simulation_data.mat','X','X_Frame','Frame_FrameIndex_Channel_Snr','stdvalues');

%Die Frage war Range und precision

scale = 100; % Bring it to [-1, 1] %will be done before pwm mudule
maxLag = 4;
Bitwidth = [8,10,12,16];
Tests_Bitwidth = length(Bitwidth);
total_blocks = floor(length(X)/stdvalues.windowlength);

Bitwidth_Frame_FrameIndex = zeros(Tests_Bitwidth,total_blocks,maxLag*2+1);

%% Bitwidht Loop
for bit_ind = 1:Tests_Bitwidth
    
    for frame_ind = 1:total_blocks
        
    X_Frame = squeeze(X(1+(frame_ind-1)*stdvalues.windowlength:(frame_ind)*stdvalues.windowlength,:)); %1.
    X_Frame = X_Frame*scale;
    
    Bitwidth_Frame_FrameIndex(bit_ind,frame_ind,:) = double(Corr_fi(X_Frame(:,1),X_Frame(:,2),maxLag,Bitwidth(bit_ind)));

    end
end

%% reference
Ref_Corr = zeros(total_blocks,maxLag*2+1);
for frame_ind = 1:total_blocks
    X_Frame = squeeze(Frame_FrameIndex_Channel_Snr(:,frame_ind,:)).*scale;

    Ref_Corr(frame_ind,:) = xcorr(X_Frame(:,1), X_Frame(:,2),maxLag);
end



%% Refernce vs Bit implementation
figure('Name','Fixed Point Test','NumberTitle','off');
for bit_ind = 1:Tests_Bitwidth %% plot Frames for Bitwidths
    
    subplot(Tests_Bitwidth+1,1,bit_ind);
    plot(squeeze(Bitwidth_Frame_FrameIndex(bit_ind,:,:))');
    legend(['Bitwidth=', num2str(Bitwidth(bit_ind))]);

end


subplot(Tests_Bitwidth+1,1,Tests_Bitwidth+1);
legend('Golden Reference');
plot(Ref_Corr');



%% Plot Section
% 
% figure();
% sgtitle('Framing a Signal');
% subplot(311);
% plot(Window);
% legend('Hann Window');
% subplot(312);
% plot(X(inLength+1:2*inLength,1)*scale);
% legend('Input Stream');
% subplot(313);
% plot((X(inLength+1:2*inLength,1))*scale.*Window);
% legend('Framed Input Stream');




%% Other

% %first microphone
% x_inp=X_Frame(1:inLength,1)'.*scale;
% %second microphone
% y_inp=X_Frame(1:inLength,2)'.*scale;

%problem, wenn voltage vor pwm zu klein ist dann ist der fehler hier groß
%aber wenn voltage zu groß ist wird die xcorr größer als 1

% eine möglichkeit wäre groß umwandeln von voltage zu pwm und dann
% im fpga runterskalieren

% anschauen was passiert, wenn man richtig niedrig geht mit smpling
% frequenz
% wie ist die Reflection abhängig von der Frequenz?


%% Dot can be used in stochastic loop to verify functionalyty
%dot_reference(curLag+maxLag+1) = dot(x_inp(1+shift_x:inLength+shift_x),y_inp(1+shift_y:inLength+shift_y));

