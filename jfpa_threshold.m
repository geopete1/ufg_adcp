function [x1_int,N_ind,screening] = jfpa_threshold(x0,time0,plotting)
%% PHASE-SPACE THRESHOLDING METHOD
% Goring & Nikora 2002, J. Hydraul. Eng. 128
%
% Input:
% x0: time series of variable to be cleaned.
% time0: time vector associated to x0.
% plotting: variable that indicates whether the plots are made or not
%   If plotting == 1, then the function makes phase plots.
%
% Output:
% x1_int: time series of corrected values of x0.
%
% Coded by J.F. Paniagua-Arroyave
% Geomorph Lab, University of Florida
% Summer 2016
%
% "The method combines three concepts: (1) that differentiation enhances 
%   the high frequency portion of a signal, (2) that the expected maximum 
%   of a random series is given by the Universal threshold, and (3) that 
%   good data cluster in a dense cloud in phase space or Poincare maps. 
%   These concepts are used to construct an ellipsoid in three-dimensional 
%   phase space, then points lying outside the ellipsoid are designated as 
%   spikes."

% Do this whilst 'index_exclude' is not empty.
index_exclude = NaN;

% Iteration number
iter = 0;

% Do screening while 'index_exclude', but max 50% of the vector length
while isempty(index_exclude) == 0 && iter < 50
    % Erase previous values
    clear index_exclude index_exclude0 index_exclude2
    
    % Iteration
    iter = iter + 1;
    
    %% 0. Demeaning, and high-pass filtering
    % Look for 10% of the Nyquist frequency for the frequency cutoff
    % If f_Nyq = 1 Hz, f_cutoff = 0.1 Hz.
    x_m = mean(x0);
    x0 = detrend(x0);
    N_x0 = length(x0);
    dt = (time0(2)-time0(1))*86400;
    f = ((1/(N_x0*dt))*(1:N_x0))';
    f_Nyq = max(f)/2;
    f_cutoff = 0.1*f_Nyq;

    % Indices for high-pass filtering (erase frequencies less than 10% of
    %   Nyquist frequency)
    indices_highpass = ...
        [find(f(1:end/2)<f_cutoff,1,'first') find(f(1:end/2)<f_cutoff,1,'last')];
% % %     indices_lowpass = ...
% % %         [find(f(1:end/2)>f_cutoff,1,'first') find(f(1:end/2)>f_cutoff,1,'last')];

    % Filtering using Fourier transform
    [x0_hp] = jfpa_fourierfilter(x0,indices_highpass,0);
% % %     [x0_lp] = jfpa_fourierfilter(x0,indices_lowpass,0);
    
    % % % figure
    % % % plot(x0,'-r'),hold on
    % % % plot(x0_f,'-b')
    % % % 
    % % % figure
    % % % plot(x0,'-r'),hold on
    % % % plot(x0,'.r')

    %% 1. Detecting the spikes
    % Surrogate for the first derivative (Dx01)
    % "Slide" vector to get the difference between u_(i+1) - u_(i-1) by
    % using the 'end' index
    Dx = (x0_hp(3:end) - x0_hp(1:end-2))/2;

    % Surrogate for the second derivative (D2x01)
    D2x = (Dx(3:end) - Dx(1:end-2))/2;
    
    % Discard first 2, and last 2 points within 'x0_hp' vector, and first
    % and last points in 'Dx' vector
    % This is due to the surrogates given by Dx(i) = ( x(i+1) - x(i-1) )/2
    x = x0_hp(3:end-2);
    Dx = Dx(2:end-1);
    
    % Use new vectors to quantify the 'expected absolute maximum'
    % eq. 1 (also, see Theorem 2 in Donoho & Johnstone 1994)
    if isequal(length(x), length(Dx), length(D2x))
        N = length(x);
    else
        fprintf('Error: lengths of vectors not equal.');
        return
    end

    % N is the same for all vectors
    % Maxima for each time series
    % Expected absolute maximum, lambda_U = sqrt(2*log(N))
    % Donoho and Johnstone 1994, Biometrika 081 after eq. 15
    max_x = sqrt(2*log(N))*std(x);
    max_Dx = sqrt(2*log(N))*std(Dx);
    max_D2x = sqrt(2*log(N))*std(D2x);

    % Rotation angle of the principal axis of D2x versus x
    % (Ellipse is rotated only for D2x versus x, other ellipses has theta = 0)
    theta = atan2(sum(x.*D2x),sum(x.^2));
    stheta = sin(theta); ctheta = cos(theta);
    
    % Axis of ellipses from maxima and minima
    % Vector of angles 'phi' covering 2*pi
    phi = -2*pi:0.01:2*pi;
    sphi = sin(phi); cphi = cos(phi);
    
        % D2x and x (solution to eq. 9)
        a_x_D2x = sqrt((max_D2x^2*(stheta)^2 - max_x^2*(ctheta)^2)/...
            ((stheta)^2-(ctheta)^2));
        b_x_D2x = sqrt((max_D2x^2*(ctheta)^2 - max_x^2*(stheta)^2)/...
            ((ctheta)^2-(stheta)^2));

            % Coordinates of ellipse
            x_x_D2x = a_x_D2x.*ctheta.*cphi ...
                - b_x_D2x.*stheta.*sphi;
            y_x_D2x = a_x_D2x.*stheta.*cphi...
                + b_x_D2x.*ctheta.*sphi;

        % Dx and x
        a_x_Dx = max_x;
        b_x_Dx = max_Dx;

            % Coordinates of ellipse (theta = 0)
            x_x_Dx = a_x_Dx.*cphi;
            y_x_Dx = b_x_Dx.*sphi;

        % D2x and Dx
        a_Dx_D2x = max_Dx;
        b_Dx_D2x = max_D2x;

            % Coordinates of ellipse (theta = 0)
            x_Dx_D2x = a_Dx_D2x.*cphi;
            y_Dx_D2x = b_Dx_D2x.*sphi;

    % Indentification of points outside the ellipses
    % Quantify y values as a function of x instead of phi
        % Ellipses major axis coinciding with x axis
        % x vs. Dx, and Dx vs. D2x
        % Add 2 to account for the lost indices at the beginning
        index_out_x_Dx = ...
            2 + find((x.^2+Dx.^2*a_x_Dx^2/b_x_Dx^2) > a_x_Dx^2);
        index_out_Dx_D2x = ...
            2 + find((Dx.^2+D2x.^2*a_Dx_D2x^2/b_Dx_D2x^2) > a_Dx_D2x^2);

        % Ellipse's major axis rotated by an angle 'theta'
        % x vs. D2x
        % Rotate x and D2x values by an angle '-theta'
        % Remember: cos(-theta)=cos(theta); sin(-theta)=-sin(theta)
        x_prime = x.*ctheta - D2x.*(-stheta);
        D2x_prime = x.*(-stheta) + D2x.*ctheta;

        % Compare rotated values with maxima expected values in the ellipse
        % Add 2 to account for the lost indices at the beginning
        index_out_x_D2x = ...
            2 + find((x_prime.^2+D2x_prime.^2*a_x_D2x^2/b_x_D2x^2) > a_x_D2x^2);

        % Complete list of indices to be potentially excluded
        index_exclude0 = ...
            sort([index_out_x_Dx; index_out_Dx_D2x; index_out_x_D2x],'ascend');

        % Find indices that are repeated within the complete list...
        unique_index = unique(index_exclude0);
        if length(unique_index) ~= 1
            % Create a histogram for values in 'index_exclude0'. In other 
            %   words, count how many times each value appears in that index.
            count_index = hist(index_exclude0,unique_index);
        else
            % If there is only one value in 'unique_index', the number of times
            %   that value appears is equal to the length of 'index_exclude0'
            count_index = length(index_exclude0);
        end

        % three times
        index_exclude3 = unique_index(count_index == 3);
        % twice
        index_exclude20 = unique_index(count_index == 2);

        % Exclude values only if:
        %   1. They are outside of 3 ellipses
        %   2. They are outside of x vs. D2x and Dx vs. D2x (a spike that is
        %   not very large in value x)
        %   3. They are outside of x vs. Dx and x vs. D2x (a spike with
        %   large curvature and large value)
        jj = 1;
        for ii = 1:length(index_exclude20)
            % Check condition 1
            if any(index_exclude20(ii) == index_out_x_D2x)==1 && ...
                    any(index_exclude20(ii) == index_out_Dx_D2x)==1

                index_exclude2(jj,1) = index_exclude20(ii);
                jj = jj+1;
                
            % Check condition 2
            elseif any(index_exclude20(ii) == index_out_x_Dx)==1 && ...
                    any(index_exclude20(ii) == index_out_x_D2x)==1

                index_exclude2(jj,1) = index_exclude20(ii);
                jj = jj+1;
            end
        end
        
        
        % If there are not indices repeated twice
        if exist('index_exclude2','var') == 0
            index_exclude2 = [];
        end

        % Arrange indices from 2 and 3 repetitions
        index_exclude = sort([index_exclude2; index_exclude3],'ascend');
        if iter == 1
            N_ind = length(index_exclude);
        end
        
    % Replacement of values
    % All indices
    index = (1:length(x0))';
    % Put 'NaNs' in indices that need to be excluded
    index(index_exclude) = NaN;
    % Store '1' in variable 'screening' if there are indices to be excluded
    if iter == 1
        if isempty(index_exclude) ~= 1
            screening = 1;
        else
            screening = 0;
        end
    end
    
    % Convert to a matrix of integers
    screening = int16(screening);
    
    % Indices and values that will be included (or not excluded)
    index_include = index(isnan(index) == 0);
    x1 = x0(index_include);
    time1 = time0(index_include);

    % Interpolate in time in blank spaces that were subtracted using the second
    %   derivative criterion
    % Interpolation using a third order polynomial (cubic, or 'pchip') with at
    %   least 12 points on either side of the spike.
    % Add the mean back of the time series
    x1_int = interp1(time1,x1,time0,'pchip','extrap') + x_m;

    % Comparison of time series (before and after discarding spikes)
    if plotting == 1 && iter < 5
        figure
        plot(x0+x_m,'-r'),hold on
        title('Screening of spikes by the phase-space thresholding method',...
            'FontName','Arial','FontSize',20)
        plot(index_exclude,x0(index_exclude)+x_m,'*r','Markersize',20)
        plot(x1_int,'-b')
        plot(x1_int,'.b')
        plot(index_exclude0,x1_int(index_exclude0),'.m','Markersize',20)
        xlabel('Indices','FontName','Arial','FontSize',20)
        ylabel('Value of x','FontName','Arial','FontSize',20)
        set(gca,'FontName','Arial','FontSize',20)
        grid on
    end

    % Assign the new values 'x1_int' to 'x0'
    x0 = x1_int;
    
    % % % figure
    % % % plot(time0,x0,'-r'),hold on
    % % % plot(time0,x0,'.r')
    % % % plot(time0,x11_int,'-b')
    % % % plot(time0(index_exclude),x11_int(index_exclude),'*r','Markersize',20)

    %% 3. Plotting of ellipses and data (discarded data in red)
    if plotting == 1 && iter < 5

        fontsi = 18;

        num_bins = sqrt(N); % Devore & Berk 2012 Math Stats (pp 16)
        % Check: Birge & Rozenholc 2006 ESAIM PS 10, eq. 1.1 and 1.2
        % Check: Wand 1997 Amer. Statist. 51

        Lim_x = [min(x) max(x)]*1.5;
        Lim_Dx = [min(Dx) max(Dx)]*1.5;
        Lim_D2x = [min(D2x) max(D2x)]*1.5;

        % x and Dx, x and D2x
        figure
        subplot(3,4,1),hold on
        title('Phase-Space Thresholding','FontSize',fontsi,'FontName','Arial')
        hist(x,num_bins)
        set(gca,'XLim',Lim_x,'FontSize',fontsi,'FontName','Arial')
        xlabel('x','FontSize',fontsi,'FontName','Arial')
        ylabel('Counts','FontSize',fontsi,'FontName','Arial')
        grid on

        subplot(3,4,5)
        plot(x,Dx,'.b','Markersize',15),hold on
        plot(x(index_out_x_Dx-2),Dx(index_out_x_Dx-2),'.r','Markersize',15)
        plot(x_x_Dx,y_x_Dx,'-r')
        set(gca,'XLim',Lim_x,'YLim',Lim_Dx,'FontSize',fontsi,'FontName','Arial')
        xlabel('x','FontSize',fontsi,'FontName','Arial')
        ylabel('\Deltax','FontSize',fontsi,'FontName','Arial')
        text(min(Lim_x)-0.1*min(Lim_x),min(Lim_Dx)-0.1*min(Lim_Dx),num2str(length(index_out_x_Dx)),'FontSize',20)
        grid on

        subplot(3,4,6)
        hist(Dx,num_bins)
        set(gca,'view',[90 -90])
        set(gca,'XLim',Lim_Dx,'FontSize',fontsi,'FontName','Arial')
        xlabel('\Deltax','FontSize',fontsi,'FontName','Arial')
        ylabel('Counts','FontSize',fontsi,'FontName','Arial')
        grid on

        subplot(3,4,9)
        plot(x,D2x,'.b','Markersize',15),hold on
        plot(x(index_out_x_D2x-2),D2x(index_out_x_D2x-2),'.r','Markersize',15)
        plot(x_x_D2x,y_x_D2x,'-r')
        xlabel('x','FontSize',fontsi,'FontName','Arial')
        ylabel('\Delta^2x','FontSize',fontsi,'FontName','Arial')
        set(gca,'XLim',Lim_x,'YLim',Lim_D2x,'FontSize',fontsi,'FontName','Arial')
        text(min(Lim_x)-0.1*min(Lim_x),min(Lim_D2x)-0.1*min(Lim_D2x),num2str(length(index_out_x_D2x)),'FontSize',20)
        grid on

        subplot(3,4,10)
        hist(D2x,num_bins)
        set(gca,'view',[90 -90])
        set(gca,'XLim',Lim_D2x,'FontSize',fontsi,'FontName','Arial')
        xlabel('\Delta^2x','FontSize',fontsi,'FontName','Arial')
        ylabel('Counts','FontSize',fontsi,'FontName','Arial')
        grid on

        % Dx and D2x
        subplot(3,4,7)
        hist(Dx,num_bins)
        set(gca,'XLim',Lim_Dx,'FontSize',fontsi,'FontName','Arial')
        xlabel('\Deltax','FontSize',fontsi,'FontName','Arial')
        ylabel('Counts','FontSize',fontsi,'FontName','Arial')
        grid on

        subplot(3,4,11)
        plot(Dx,D2x,'.b','Markersize',15),hold on
        plot(Dx(index_out_Dx_D2x-2),D2x(index_out_Dx_D2x-2),'.r','Markersize',15)
        plot(x_Dx_D2x,y_Dx_D2x,'-r')
        set(gca,'XLim',Lim_Dx,'YLim',Lim_D2x,'FontSize',fontsi,'FontName','Arial')
        xlabel('\Deltax','FontSize',fontsi,'FontName','Arial')
        ylabel('\Delta^2x','FontSize',fontsi,'FontName','Arial')
        text(min(Lim_Dx)-0.1*min(Lim_Dx),min(Lim_D2x)-0.1*min(Lim_D2x),num2str(length(index_out_Dx_D2x)),'FontSize',20)
        grid on

        subplot(3,4,12)
        hist(D2x,num_bins)
        set(gca,'view',[90 -90])
        set(gca,'XLim',Lim_D2x,'FontSize',fontsi,'FontName','Arial')
        xlabel('\Delta^2x','FontSize',fontsi,'FontName','Arial')
        ylabel('Counts','FontSize',fontsi,'FontName','Arial')
        grid on
    end
end
end