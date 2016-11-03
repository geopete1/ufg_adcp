function [x_f] = jfpa_fourierfilter(x,indices,plotting)
x = detrend(x);
N = length(x);

% % % [X,f]=jfpa_fourier(t,x); % Fourier transform of x

% Fast Fourier transform
X = fft(x);

% User input
% SELECT INDICES TO EXCLUDE FROM FOURIER TRANSFORM
if nargin < 2
    % Ask user for the indices
    prompt_1 = 'Select initial index: ';
    n_1 = input(prompt_1);

    prompt_2 = 'Select final index: ';
    n_2 = input(prompt_2);
else
    % Indices are provided when the function is called
    n_1 = min(indices);
    n_2 = max(indices);
end

% Indices for the first half of the transform
index01 = n_1:n_2;

% Indices for the second half of the transform
delta_index01 = index01-1;
delta_index02 = delta_index01-1;
index02 = N-delta_index02;

% Indices to exclude (both sides)
index0 = [index01 index02(end:-1:1)]';

if index0(end)>N
    index0 = index0(1:1:end-1);
end

% Kind of a Heaveside function (one everywhere except where we say it is
% equal to zero, in this case between frequencies 0.1 and 0.25 Hz)
% We need to get rid of both sides of spectrum
% Look for the indices you want to substract from the plot of the absolute
% values of X
H = ones(length(X),1);
H(index0,1) = 0;


% Convolution of that H function with X
X_p = X.*H;

% Comparison of FFT before and after filtering
if plotting == 1
    figure
    title('Filtering of x using Fourier transform')
    plot(abs(X),'-b'),hold on
    plot(abs(X_p),'-r')
    xlabel('Index of frequency')
    legend('Initial time series','After filtering')
    grid on
end

% Inverse Fourier transform to X and Xp (first gives us x back, second
% generates a time series of x without the frequencies we took out)
x_f = ifft(X_p);