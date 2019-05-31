function [filtered_trace, event_times, event_sizes] = wiener_filter_post_bessel(trace, unfilt, params, template, nfft)
%
% ex = wienerFilter(y,h,sigma,gamma,alpha);
%
% Generalized Wiener filter using parameter alpha. When
% alpha = 1, it is the Wiener filter. It is also called
% Regularized inverse filter.
%
% Reference: Richb's paper
% Created: Tue May 4 16:24:06 CDT 1999, Huipin Zhang

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

event_sizes = zeros(size(event_times));

if params.event_sign == -1
    for j = 1:length(event_times)
        baseline = mean(unfilt(max(1,event_times(j)-40):event_times(j))); %changed this to make amplitude more realistic
        event_sizes(j) = (min(unfilt(event_times(j):min(event_times(j)+ min_window,length(unfilt)))) - baseline) * params.event_sign; %edited this to make the amplitude more realistic
        if min(unfilt(event_times(j):min(event_times(j)+min_window, length(unfilt)))) > baseline
            event_sizes(j) = 0;
        end
    end
else
    for k = 1:length(event_times)
        baseline = mean(unfilt(max(1,event_times(k)-40):event_times(k)));
        event_sizes(k) = (max(unfilt(event_times(k):min(event_times(k)+min_window, length(unfilt)))) - baseline);
        if max(unfilt(event_times(k):min(event_times(k)+min_window, length(unfilt)))) < baseline
            event_sizes(k) = 0;
        end
    end
end

event_times(event_sizes < 2) = [];
event_sizes(event_sizes < 2) = [];

end

