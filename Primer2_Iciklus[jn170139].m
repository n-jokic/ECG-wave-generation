clc;
close all;
clear all;

clear all
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% QRS pulse


fs =  200; %default value, Hz
bpm =  72; %default value, beats per minute
amp = 1000; %default value, micro volts
qrswidth = 0.120; %default value

QRS_wave =  QRSpulse(qrswidth,bpm,fs,amp)/1000;

n = 0: (length(QRS_wave)-1);
t = n/fs;

figure();
    plot(t, QRS_wave, 'black');
    hold on;
    stem(t,QRS_wave, 'black.','MarkerSize', 7, 'LineWidth', 0.4);
        xlim([0, t(end)]);
        title('Synthetic EKG');
        xlabel('t[s]');
        ylabel('Amplitude [mV]');
        grid('on');
        legend({['continuous'], ['samples']})
    hold off;
        
 saveas(gcf,['C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\pulse'],'epsc')  ; 
        
%d=0.070; %.07 to .135 seconds, QRS width

qrswidth = [0.07,0.1, 0.135];
figure();
N = length(qrswidth);

for k = 1:length(qrswidth)
    subplot(N,1, k);
    
    QRS_wave =  QRSpulse(qrswidth(k),bpm,fs,amp)/1000;
    n = 0: (length(QRS_wave)-1);
    t = n/fs;
    
    plot(t,QRS_wave,'black');
        xlim([0, t(end)]);
        title([ 'Synthetic EKG']);
        legend(['d = ' num2str(qrswidth(k)) ]);
        ylim([-0.2, 1]);
        xlabel('t[s]');
        ylabel('Amplitude [mV]');
        grid('on');
    
end 
saveas(gcf,['C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\pulse_d_normal'],'epsc');

qrswidth = [0.01, 0.335];

N = length(qrswidth);

for k = 1:length(qrswidth)
    figure();
    
    QRS_wave =  QRSpulse(qrswidth(k),bpm,fs,amp)/1000;
    n = 0: (length(QRS_wave)-1);
    t = n/fs;
    
    plot(t,QRS_wave,'black');
        xlim([0, t(end)]);
        title([ 'Synthetic EKG']);
        legend(['d = ' num2str(qrswidth(k)) ]);
        ylim([-0.2, 1]);
        xlabel('t[s]');
        ylabel('Amplitude [mV]');
        grid('on');
    name = ['C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\pulse_d_' num2str(k)];
    saveas(gcf,name,'epsc');
    
end    


fs =  200; %default value, Hz
bpm =  72; %default value, beats per minute
amp = 1000; %default value, micro volts
qrswidth = 0.120; %default value

QRS_wave =  QRSpulse_modif(qrswidth,bpm,fs,amp)/1000;

n = 0: (length(QRS_wave)-1);
t = n/fs;

figure();
    plot(t, QRS_wave, 'black');
        xlim([0, t(end)]);
        title('Synthetic EKG, T wave first');
        ylim([-0.2, 1]);
        xlabel('t[s]');
        ylabel('Amplitude [mV]');
        grid('on');
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\T_first','epsc');
        

%% ECGwaveGen, signal 1
t = 11; %s
bpm = 60;
amp = 1000; % uV
fs = 64; %ss

Norm = ECGwaveGen(bpm, t, fs, amp);
Abnorm = ECGwaveGen(bpm*1.5, t, fs, amp);

time = [(1:length(Norm)) / fs];

figure;
    subplot(2,1,1)
        plot(time, Norm, 'black');
            title('ECGwaveGen synthesis - normal');
            xlabel('t [s]');
            ylim([-50, 900]);
            xlim([0, time(end)]);
            grid on;
            ylabel('Amplitude [\muV]');
    subplot(2,1,2)
        time = [(1:length(Abnorm)) / fs];
        plot(time, Abnorm, 'black');
            title('ECGwaveGen synthesis - tachycardia');
            xlabel('t [s]');
            ylim([-50, 900]);
            xlim([0, time(end)]);
            grid on;
            ylabel('Amplitude [\muV]');
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\wave_tachy','epsc');

%% ecgsyn 
fs = 64;
N = 10;
%sfecg,N,Anoise,hrmean,hrstd,lfhfratio
s = ecgsyn(fs, N, 0, bpm); 
%s = ecgsyn();

t = (0:length(s)-1)/fs;
figure;
    plot(t, s, 'black')
       title('ecgsyn synthesis - 4 QRS complexes');
       xlabel('t [s]');
       xlim([0.5,4.5]);
       grid('on')
       ylabel('Amplitude [mV]'); 
       
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\syn','epsc');

    
figure;
    subplot(2,1,1)
        stem(t(33 : 225), s(33 :225), 'black.', 'MarkerSize', 2,'LineWidth',0.01)
            title('ecgsyn synthesis - 4 QRS complexes');
            xlabel('t [s]');
            xlim([t(33),3.5]);
            grid('on')
            legend(['fs = ' num2str(fs) ' Ss'])
            ylabel('Amplitude [mV]'); 
            
    subplot(2, 1, 2)
    b = t;
    s = ecgsyn(3*fs, N, 0, bpm,1 ,0.5 , 3*512);
    t = (0:length(s)-1)/fs/3;

        stem(t(97 : 673), s(97 :673), 'black.', 'MarkerSize', 2,'LineWidth',0.01)
            title('ecgsyn synthesis - 4 QRS complexes');
            xlabel('t [s]');
            xlim([t(97),3.5]);
            legend(['fs = ' num2str(3*fs) ' Ss'])
            grid('on')
            ylabel('Amplitude [mV]'); 
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\syn_fs','epsc');
       
 %% Treci zad
 
fs =  360; %default value, Hz
bpm =  72; %default value, beats per minute
amp = 1000; %default value, micro volts
qrswidth = 0.120; %default value

QRS_wave = QRSpulse(qrswidth,bpm,fs,amp)/1000;
n = 0: (length(QRS_wave)-1);
t = n/fs;


%vremenski domen

endp = round(10*fs);

figure();
    subplot(3, 1, 1)
        plot(t, QRS_wave, 'black');
            xlim([0, 10]);
            title('QRSpulse synthetic ECG');
            xlabel('t[s]');
            ylabel('Amplitude [mV]');
            ylim([-0.5, 1]);
            grid('on');
            
    subplot(3, 1, 2)
    
    t = 11; %s
    bpm = 60;
    amp = 1000; % uV

    Norm = ECGwaveGen(bpm, t, fs, amp)/1000;
    t = (0 : (length(Norm)-1))/fs;
    
        plot(t, Norm, 'black');
            xlim([0, 10]);
            title('ECGwaveGen synthetic ECG');
            xlabel('t[s]');
            ylabel('Amplitude [mV]');
            ylim([-0.5, 1]);
            grid('on');
            
    subplot(3, 1, 3)
    
    N = 20;
    bpm = 60;

    s = ecgsyn(fs, N, 0, bpm,1 ,0.5 , 3*360); 
    t = (0 : (length(s)-1))/fs;
    
        plot(t, s, 'black');
            xlim([0, 10]);
            title('ecgsyn synthetic ECG');
            xlabel('t[s]');
            ylabel('Amplitude [mV]');
            ylim([-0.5, 1]);
            grid('on');
            
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\SECG_time_comp','epsc');

%frekv domen   
Nfft =  16384;
clear sfft;
clear Normfft;
clear QRS_wavefft;

QRS_wavefft = abs(fft(QRS_wave, Nfft)/length(QRS_wave))*2;
Normfft = abs(fft(Norm(1:endp), Nfft)/length(Norm(1:endp)))*2;
sfft = abs(fft(s(1:endp), Nfft)) / length(s(1:endp))*2;

figure();
    subplot(3, 1, 1)
    
    n = 0: (length(QRS_wavefft)-1);
    t = n*fs/Nfft;
    
        plot(t, QRS_wavefft, 'black');
            xlim([0.1, fs/8]);
            title('QRSpulse synthetic ECG amplitude spectrum');
            xlabel('f[Hz]');
            ylabel('Amplitude [mV]');
            ylim([0, 0.2]);
            grid('on');
            
    subplot(3, 1, 2)
    

    t = (0 : (length(Normfft)-1))*fs/Nfft;
    
        plot(t, Normfft, 'black');
            xlim([0.1, fs/8]);
            title('ECGwaveGen synthetic ECG amplitude spectrum');
            xlabel('f[Hz]');
            ylabel('Amplitude [mV]');
            ylim([0, 0.2]);
            grid('on');
            
    subplot(3, 1, 3)
    
    t = (0 : (length(sfft)-1))*fs/Nfft;
    
        plot(t, sfft, 'black');
            xlim([0.1, fs/8]);
            title('ecgsyn synthetic ECG amplitude spectrum');
            xlabel('f [Hz]');
            ylabel('Amplitude [mV]');
            ylim([0, 0.2]);
            grid('on');    
            
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\SECG_spec_comp','epsc');

%% Merenje vremena, treci zad

total = 30;
t_wave = zeros(1,total);
t_pulse = zeros(1,total);
t_syn = zeros(1,total);

for i = 1 : total
    tic;
    QRS_wave = QRSpulse(qrswidth,bpm,fs,amp)/1000;
    t_pulse(i) = toc*1000; %ms
end

disp(['Pulse runtime: ' num2str(mean(t_pulse)) ' +- ' num2str(2*std(t_pulse)) ' ms']);

for i = 1 : total
    tic;
    Norm = ECGwaveGen(bpm, 10, fs, amp)/1000;
    t_wave(i) = toc*1000; %ms
end

disp(['Wave runtime: ' num2str(mean(t_wave)) ' +- ' num2str(2*std(t_wave)) ' ms']);
N = 60;
for i = 1 : total
    tic;
    s = ecgsyn(fs, N, 0, bpm,1 ,0.5 , 3*fs);
    t_syn(i) = toc*1000; %ms
end

disp(['Syn runtime: ' num2str(mean(t_syn)) ' +- ' num2str(2*std(t_syn)) ' ms']);
%% 4. zadatak 
s = ecgsyn(); 

fs = 256;

t = (0 : length(s)-1)/fs;
noise1 = sin(2*pi*0.5*t);
noise2 = 0.5*sin(2*pi*5*t);
noise3 = 0.1*sin(2*pi*50*t);
n = length(s);
noise4 = randn(1,n)*0.1/3;

%a)
n = 0 : length(s)-1;
[pk, lk] = findpeaks(s,n,'MinPeakDistance', 120, 'MinPeakHeight',1);

figure();
    plot(n/256, s,'black',lk/256,pk,'ro');
        grid on;
        title('ecgsyn synthetic ECG R peak detection');
        xlim([1, 11]);
        xlabel('t[s]');
        ylim([-2, 2]);
        ylabel('Amplitude [mV]');

s = s' + noise1 + noise2 + noise3 + noise4;
[pk, lk] = findpeaks(s,n,'MinPeakDistance', 120, 'MinPeakHeight',1);

saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\synR_det','epsc');

figure();
    plot(n/256, s,'black',lk/256,pk,'ro');
        grid on;
        title('ecgsyn synthetic ECG R peak detection, with noise');
        xlim([1, 11]);
        xlabel('t[s]');
        ylabel('Amplitude [mV]');
        ylim([-2, 2]);
saveas(gcf,'C:\Users\milos\OneDrive\VII semestar\MAS\vezbe\vezba 2\izvestaj\synR_det_noise','epsc');
