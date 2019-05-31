function [filtered_trace, event_times] = wiener_filter_PPR( trace, params, template, nfft)

%
% Reference: Richb's paper
% Created: Tue May 4 16:24:06 CDT 1999, Huipin Zhang
% Altered by CRW: 8 March 2019

dt = params.dt;
ar_noise_params = params.init_method.ar_noise_params;
threshold = params.init_method.threshold;
min_window = params.init_method.min_interval;

Fs = 1/dt;

N = length(trace);

%Fast fourier transforms
   %frequency domain representation = fft(input, transform length)
trace_f = fft(trace,nfft); % replace with multitaper
template_f = fft(template,nfft);

%[PSD estimatea, cyclical freq vector] = pmtm(input signal, time-halfbandwidth product,  # of DFT points, sample rate)
[trace_P, freq_mt] = pmtm(trace,9,nfft*2,Fs); 

trace_P = trace_P(1:end-1)';
trace_P = trace_P/sum(trace_P);

% get parameterized ar noise psd
%frequency response of digital filter
[noise_P, freqs_noise_P] = freqz(sqrt(ar_noise_params.sigma_sq),ar_noise_params.phi,nfft,Fs);

% min(freqs_noise_P)
% max(freqs_noise_P)

% direct implementation of the regularized inverse filter, 
% when alpha = 1, it is the Wiener filter
% Gf = conj(template_f).*Pxf./(abs(template_f.^2).*Pxf+alpha*noise_P);
%
% Since we don't know Pxf, the following 
% handle singular case (zero case)
% template_f_clean = template_f.*(abs(template_f)>0)+1/gamma*(abs(template_f)==0);
template_f_clean = template_f;
inverse_template_f = 1./template_f_clean;
% inverse_template_f = inverse_template_f.*(abs(template_f)*gamma>1)+gamma*abs(template_f_clean).*inverse_template_f.*(abs(template_f_clean)*gamma<=1);



noise_P = noise_P'/sum(noise_P);
% size(trace_P)
% size(noise_P)
first_i = find(trace_P<=noise_P,1,'first');
if isempty(first_i)
    first_i = nfft;
end
% trace_P = trace_P.*(trace_P>noise_P)+noise_P.*(trace_P<=noise_P);
trace_P = trace_P.*((1:nfft)<first_i)+noise_P.*((1:nfft)>=first_i);
wien_filter = inverse_template_f.*(trace_P-noise_P)./(trace_P); %in denom: -(1-alpha)*noise_P


% loglog(freqs_noise_P,abs(trace_P),'b'); hold on;
% % loglog(freqs_noise_P,abs(trace_P_period),'g'); hold on;
% loglog(freqs_noise_P,abs(noise_P),'r'); hold on;

% max(max(abs(Gf).^2)) % should be equal to gamma^2
% Restorated image without denoising
filtered_trace = wien_filter.*trace_f;
filtered_trace = real(ifft2(filtered_trace));
filtered_trace = filtered_trace(1:length(trace));


[~, event_times] = findpeaks(zscore(filtered_trace),'MinPeakHeight',threshold,'MinPeakDistance',min_window);


end

