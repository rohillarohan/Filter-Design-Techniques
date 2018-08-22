%% Filter Design Using Different Techniques
clc
close all
clearvars

%% Specifications

sampFreq = 20000; % Sampling Frequency
pbFreq = 4000; % Pass Band Frequency
sbFreq = 4500; % Stop Band Frequency
pbRipl = 0.1; % Pass Band Ripple
sbRipl = 0.05; % Stop Band Ripple
delta = min(pbRipl,sbRipl); % Finds Minimum Ripple
fCuts = [pbFreq sbFreq]; % Band edges
normFreq = fCuts/(sampFreq/2); % Normalising the frequency

%% Kaiser Window Method

mags = [1 0]; % For low pass filter
devs = [pbRipl sbRipl]; % Pass-band and Stop-band ripple

[M_k, Wn, beta, ftype] = kaiserord(normFreq, mags, devs); 
% M_k = Filter order(for Kaiser window), Wn = Normalised band edge, beta, ftype = Filter Type

if mod(M_k,2) == 1 % Checks the parity of the order and makes it even
    M_k = M_k+1;
end
hK = fir1(M_k,Wn,ftype,kaiser(M_k+1,beta),'noscale');
[HK,wK] = freqz(hK,1);
freqK = 0:length(HK)-1;
fhzK = (freqK*(sampFreq/2))/length(HK); % Converting the frequency in Hertz
figure;
plot(fhzK,20*log10(abs(HK)));
grid on;
title('Frequency Response(dB) for Kaiser Window Method');
ylabel('Magnitude (dB)');
xlabel('Frrequency (Hz)');
fvtool(hK, 'polezero'); title('Pole-Zero Plot for Kaiser Window Method')% Pole-Zero Plot for Kaiser Window Method

%% Parks-McClellan Method

[M_pm,fo,ao,W] = firpmord(normFreq, mags, devs);
%Computes filter order M_pm, Normalized Frequency Band Edges fo, Amplitude Response in ao, and Band Weights W

if mod(M_pm,2) == 1
    M_pm = M_pm+1;
end

[hPm,deltaNew]=firpm(M_pm,fo,ao,W);
%Determines the Impulse Response hPm of the Equiripple Filter

% Checks the Stop Band Ripple Height and Changes the Order of the Filter Accordingly
while min(delta, deltaNew) >= delta 
    M_pm = M_pm + 2;
    [hPm,deltaNew]=firpm(M_pm,fo,ao,W);
end

[HPm,wPm] = freqz(hPm,1);
freqPm = 0:length(HPm)-1;
fhzPm = (freqPm*(sampFreq/2))/length(HPm); % Converting the frequency in Hertz
ampPm = abs(HPm); % Finding the Amplitude 
figure;
subplot(2,1,1)
plot(fhzPm,20*log10(abs(HPm)));
grid on;
title('Frequency Response(dB) for Parks-McClellan Algorithm');
ylabel('Magnitude (dB)');
xlabel('Frrequency (Hz)');
subplot(2,1,2)
plot(fhzPm,ampPm);
title('Amplitude Response from Parks-Mcclellen Algorithm');
ylabel('Magnitude');
xlabel('Frrequency (Hz)');
grid on;
fvtool(hPm,'polezero');title('Pole-Zero Plot for Parks Mccllelen') % Pole-Zero Plot for Parks Mccllelen

%% Analysis (Finding the Minimum Stop Band Ripple)

for n = 1:length(HPm) 
    if ampPm(n) < 10^(-3)
        delmin = max(ampPm(n:length(HPm)));
        break
    end
end

%% Minimum Phase Filter Design

hPm(M_pm/2+1) = hPm(M_pm/2+1) + delmin; 
% Adding the Stop Band Ripple of the Amplitude Response in the Middle Sample 

[HPm2, wpm2] = freqz(hPm,1);
freqpm2 = 0:length(HPm2)-1;
fhzPm2 = (freqpm2*(sampFreq/2))/length(HPm2); % Converts the frequency in Hertz
amp2 = abs(HPm2); % Finds the Amplitude 
zero = roots(hPm);
insUc = find(abs(zero) <= 1.000); % Finds the Zeros On and Inside the Unit Circle

hMinIns = poly(zero(insUc)); 
% Construts the Impulse Response for Minimum Phase FIlter from the (New) Zeros

hMin = sqrt(hPm(1)/prod(abs(roots(hMinIns(1:length(hMinIns))))))*hMinIns; 
% Scales the Impulse Response

[Hmin,W] = freqz(hMin,1);
freqmin = 0:length(Hmin)-1;
fhzmin = (freqmin*(sampFreq/2))/length(Hmin); % Converts the frequency in Hertz
ampmin = abs(Hmin); % Finds the Amplitude
figure;
plot(fhzmin,ampmin); 
hold on;
plot(fhzPm2,amp2);
title('Amplitude Response Comparision from Parks-Mcclellen Algorithm with Minimum Phase Filter');
ylabel('Magnitude');
xlabel('Frrequency (Hz)');
grid on;
legend('Minimum Phase Filter', 'Linear Phase Filter');
fvtool(hMin, 'polezero'); title('Pole-Zero Plot for Minimum Phase Filter') % Pole-Zero Plot for Minimum Phase Filter