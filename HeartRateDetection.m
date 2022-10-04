%% Importing ECG signal
time=10;
fid = fread 
ecg_signal=fread(fopen('data/100.dat'),2*360*time,"float64");
Orig_Sig=ecg_signal(1:2:length(ecg_signal));
plot(time,Orig_Sig);
title("Trial 1");
plot((0.1/36):(0.1/36):0.1,Orig_Sig);
title('Acquired ECG signal');
%% Normalising the ecg signal
ecg_normalised = (ecg_signal - mean(ecg_signal))/max(ecg_signal);
Orig_Sig=ecg_normalised(1:2:length(ecg_normalised));
figure;
plot((0.1/36):(0.1/36):0.1,Orig_Sig);
%% Bandstop filter
% Powerline interference (50 or 60 Hz noise from mains supply) can be removed by using 
% a notch filter of 50 or 60 Hz cut-off frequency.
y = bandstop(ecg_normalised,[49 51],360);
Orig_Sig=y(1:2:length(y));
figure;
plot((0.1/36):(0.1/36):0.1,Orig_Sig);
%% Removing Baseline Wander
% Baseline wander is a low-frequency noise of around 0.5 to 0.6 Hz. 
% To remove it, a high-pass filter of cut-off frequency 0.5 to 0.6 Hz can be used. 
Fstp = 0.5;
Fs = 360;

dbutter = designfilt('highpassiir','FilterOrder',2, 'StopbandFrequency',Fstp,'SampleRate',Fs);

y1 = filter(dbutter,y);

Orig_Sig=y1(1:2:length(y1));
figure;
plot((0.1/36):(0.1/36):0.1,Orig_Sig);
title('ECG signal post noise removal');
%% Calculating R-peaks

% High Pass filter with cutoff = 15Hz
Fstp = 15;
Ap = 1;
Ast = 60;
Fs = 360;

dbutter_h = designfilt('highpassiir','FilterOrder',4, 'StopbandFrequency',Fstp,'SampleRate',Fs);
y2 = filter(dbutter_h,y1);

% % Low Pass filter with cutoff = 20Hz

Fp = 20;
Ap = 1;
Ast = 60;
Fs = 360;

dbutter_l = designfilt('lowpassiir','FilterOrder',4, 'PassbandFrequency',Fp,'SampleRate',Fs);
y3 = filter(dbutter_l,y2);

% Orig_Sig=y3(1:2:length(y3));
% figure;
% plot((0.1/36):(0.1/36):0.1,Orig_Sig)

% Calculating Energy Signal

y4 = y3.*y3;
%Integrator filter to smoothen out energy signal
y4 = filter(1,[1 -0.8],y4);
Orig_Sig=y4(1:2:length(y4));
figure;
plot((0.1/36):(0.1/36):0.1,Orig_Sig);
title('Energy Signal (smmothened through integrator)');

%% Autocorrelation
y_acf = xcorr(y4,y4);
[yupper,ylower] = envelope(y_acf,100,'peak');
figure;
plot(y_acf);
hold;
plot(yupper);
legend('Autocorrelation of Energy signal','Smooth envelope to detect peaks');
title('Autocorrelation of Energy signal with Peaks detected');

% Finding peaks in Autocorrelation function

[st, maximas, minimas, ed] = findextremas(yupper);
peaks_acc = maximas(:,1);
peaks_acc = peaks_acc(maximas(:,2)>=10^(-10));
peaks_acc = peaks_acc(1:floor(length(peaks_acc)/2));

% Calculating RR-intervals and Heart rate from timestamps of peaks

k = (10/7200)*peaks_acc;
[l,~] = size(k);
RR_acc = zeros(l-1,1);

for i=1:l-1
    RR_acc(i) = 60/(k(i+1)-k(i));
end

stamps_acc = k(1:end-1);

figure;
plot(stamps_acc,RR_acc);
title('Heart Rate trend');
grid;
axis on;
%% Envelope of Energy Signal

[yupper,ylower] = envelope(y4,10,'peak');
figure;
plot(y4);
hold;
plot(yupper);
legend('Energy Signal','Envelope of Energy Signal');
title('Envelope of energy signal for Peak Detection');

% Finding peaks in the envelope of energy

[st, maximas, minimas, ed] = findextremas(yupper);
peaks = maximas(:,1);
peaks = peaks(maximas(:,2)>=10^(-6));
peaks = peaks(2:end);

% Calculating RR-intervals and Heart rate from timestamps of peaks

k = (10/7200)*peaks;
[l,~] = size(k);
RR = zeros(l-1,1);

for i=1:l-1
    RR(i) = 60/(k(i+1)-k(i));
end

stamps = k(1:end-1);

% figure;
% plot(stamps,RR);
% hold;
% plot(stamps_acc,RR_acc);
% title('Heart Rate trend');
% grid;
% axis on;
% legend('Envelope Method','Autocorrelation method');

%% Threshold Method

TH = 2*mean(y4);
for i=1:length(y4)
    if y4(i)>TH
        i_start = i;
        break;
    end
end
count = 1;
peaks_th = zeros(12,1);
prev_i = i_start;
for i=i_start:length(y4)
    if i-prev_i>=470 && y4(i)>TH
        peaks_th(count) = i;
        prev_i = i;
        i=i+720;
        count = count+1;
    end
end

k = (10/7200)*peaks_th;
[l,~] = size(k);
RR_th = zeros(l-1,1);

for i=1:l-1
    RR_th(i) = 60/(k(i+1)-k(i));
end

stamps_th = k(1:end-1);

y4 = filter(1,[1 -0.8],y4);
Orig_Sig=y4(1:2:length(y4));
figure;
plot((0.1/36):(0.1/36):0.1,Orig_Sig);
hold on;
plot(k,y4(peaks_th),'*');
hold off;
title('Points obtained after threshold filtering');


figure;
plot(stamps,RR);
hold on;
plot(stamps_th,RR_th);
plot(stamps_acc,RR_acc);
hold off;
title('Heart Rate trend');
grid;
axis on;
legend('Envelope Method','Threshold Method','Autocorrelation method');