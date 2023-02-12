function [mu,kappa,p,S,Diff,iter] = nl_lse_mldrw(curve, time, init_guess, maxiter, tolerance, display)


if isequal(display,'on')
    figure;
end

% Initilization variables
mu = zeros(1,maxiter+1);
kappa = zeros(1,maxiter+1);
p = zeros(length(time),maxiter); % probablity distribution
S = zeros(1,maxiter+1); % Sum of squared
Diff = inf*ones(1,maxiter+1); % Difference in sum of squared between two iterations

% Assign kappanown values to variables
iter = 1; % iteration
mu(1,iter) = init_guess(1);
kappa(1,iter) = init_guess(2);
AUC = init_guess(3);
alpha = 0.1; % learning rate of Gauss-Newton method
t0 =0;
% We estimate the probability distribution after injection time
start_t = min(find(time-t0 > 0))+1;
Time = time(start_t:end);
while iter <= maxiter
    
    % Probability distribution of the trasmit time given in Eq 2
    p(start_t:end,iter) = sqrt(kappa(1,iter)./(2*pi*(Time))).*exp(-kappa(1,iter)*(Time-mu(1,iter))...
        .^2./(2*(Time)));
    
    % First derivates in terms of mu and kappa
    dpdmu = p(start_t:end,iter).*kappa(1,iter).*(1-mu(1,iter)./(Time));
    
    dpdkappa = p(start_t:end,iter).*0.5.*(1./kappa(1,iter) - (Time-mu(1,iter)).^2./(Time));
    
    % Prepare the Jacobian matrix
    J = [dpdmu dpdkappa];
    
    % Calculate the sum of square of residuals for the initial guess
    S(1,iter) = sum((curve - p(:,iter)).^2);
    
    % Calculate the difference of the sum of square between current and
    % periouvs iterations
    if iter > 1
        Diff(:,iter) = abs(S(1,iter-1) - S(1,iter));
    elseif iter == 1
        Diff(:,iter) = +inf;
    end
    
    % Stop iteration when the difference is smaller than the tolerance
    if Diff(:,iter) < tolerance
        return;
    end
    
    % least square estimation
    Theta = [mu(1,iter);kappa(1,iter)] + alpha*inv(transpose(J)*J)*transpose(J)*(curve(start_t:end,1) - p(start_t:end,iter));
    mu(1,iter+1) = Theta(1,1);
    kappa(1,iter+1) = Theta(2,1);
    
    % Update estimation curve
    if isequal(display,'on')
        plot(time,AUC*curve,'b*')
        hold on
        plot(time,AUC*p(:,iter),'r')
        hold off
        xlabel('Time (s)')
        ylabel('Video Intensity')
        legend('Real Time-Intensity curve', 'Estimated Time-Intensity curve')
        title(['Iteration ', num2str(iter), '. With mu = ', num2str(mu(1,iter)), ' and kappa = ',...
            num2str(kappa(1,iter))])
        drawnow
        pause(0.01)
    end
    iter = iter + 1;
end

end