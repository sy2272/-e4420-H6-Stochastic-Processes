clc; clear all; close all;
load hrs.mat;
N = length(hrs);
P = 10; % set AR model order

%% Autocorrelation
% compute AR model using autocorrelation (solve for a1)
autocorrx = xcorr(hrs);
[m n] = size(autocorrx);
r0 = (m+1)/2;
rx = autocorrx(r0+1:r0+P);
Rxx = toeplitz(autocorrx(r0:r0+P-1));
a1 = Rxx\(-1*rx);

% compute signal power spectrum
Px = (1/N) * (abs(fft(hrs,N))).^2;

% compute frequency response
pSpecL = (N+1)/2;
[h,w] = freqz(1,[1 a1'],pSpecL);
Ph = abs(h).^2;

% plot signal power spectrum and AR model power spectrum
figure(),
plot(w,10*log10(Px(1:pSpecL))); hold on;
plot(abs(w),10*log10(Ph),'r--');
legend('signal power spectrum', 'AR model power spectrum');
title(['Power spectra of signal and AR model, P = ' num2str(P)]);
axis([0 pi -50 50]); xlabel('Normalized frequency'); ylabel('10*log10(P)');
saveas(gcf,'autocorr_Ph.jpg');

% compute innovation process with inverse filter (solve for Pv)
v = filter(-1*[1 a1'],1,hrs);
Pv = (1/N) * (abs(fft(v,N))).^2;

% plot innovation process and its power spectrum
figure;
subplot(2,1,1);plot(v);title('Innovation process');axis([1 N -50 50]);
subplot(2,1,2);plot(w,10*log10(Pv(1:pSpecL)));title('Power spectrum of innovation process');
axis([0 pi -50 50]); ylabel('10*log10(P)'); xlabel('Normalized frequency');
saveas(gcf,'autocorr_Pv.jpg');

%% Linear Prediction
% compute AR model with lpc
a = lpc(hrs,P);

% compute linear prediction with AR model
est_x = filter(-1*[0 a(2:end)],1,hrs);

figure;
plot(hrs);hold on;plot(est_x,'r-');
title(['Signal and linear prediction estimation result, P = ' num2str(P)]);
legend('Signal', 'Linear prediction');axis([1 N 50 100]);
saveas(gcf,'lpc_estx.jpg');

%% LPC Estimation Error
% show linear prediction estimation error
err = est_x - hrs;
figure;
subplot(2,1,1);plot(err);
title('Linear prediction estimation error');
axis([1 N -50 50]);
saveas(gcf,'lpc_err.jpg');

% estimated error power spectrum
Perror = (1/N) * (abs(fft(err,N))).^2;
subplot(2,1,2);plot(w,10*log10(Perror(1:pSpecL)));
axis([0 pi -50 50]);
title(['Linear prediction estimation error power spectrum , P = ' num2str(P)]);
ylabel('10*log10(P)');xlabel('Normalized frequency');
saveas(gcf,'lpc_Perr.jpg');

% Compute Signal to Noise Ratio of LPC model
SNR=sum(est_x.^2)/sum(err.^2);
% SNR for P=10 is 348.69

% SNR in dB as function of model order P
P=100;
SNR=[];
dBSNR=[];

for i=1:P
    % compute AR 
    a=lpc(hrs,i);
    
    % compute linear prediction with AR model
    est_x = filter(-1*[0 a(2:end)],1,hrs);
   
    % solve for e
    err=hrs-est_x;
    
    %solve for SNR in dB
    SNR(i)=sum(est_x.^2)/sum(err.^2);
    dBSNR(i)=db(SNR(i));
    
end

figure(),plot(1:100,dBSNR)
title('SNR in dB of LPC');
ylabel('SNR in dB'); xlabel(['model order, P=' num2str(P)]);
saveas(gcf,'lpc_dbSNR.jpg');

