function paper(s_0, s_1, p_0, L)
    
    p_1 = 1 - p_0;
    L_star = (s_0 .* p_0 - s_1 .* p_1) ./ (p_1 - p_0);

    print('p_1', p_1)
    print('p_0', p_0)
    print('L_star', L_star)

    lambda_0 = (s_0 + L).^2;
    lambda_1 = (s_1 + L).^2;

    lambda_avg = lambda_0 .* p_0 + lambda_1 .* p_1;
    lambda_wave = sqrt(lambda_0) .* p_0 + sqrt(lambda_1) .* p_1;

    I = - lambda_avg .* log(lambda_avg) ...
        + p_0 .* lambda_0 .* log(lambda_0) ...
        + p_1 .* lambda_1 .* log(lambda_1);

    I(abs(I) < 1e-15) = 0;

    % dI_old = p_0 .* sqrt(lambda_0) .* log2(lambda_0/lambda_avg) + ...
    %          p_1 .* sqrt(lambda_1) .* log2(lambda_1/lambda_avg);

    dI_new = p_0.* sqrt(lambda_0).*log(lambda_0) ...
             + p_1.* sqrt(lambda_1).*log(lambda_1) ...
             - log(lambda_avg).*lambda_wave;
    
    [maxVal, maxIndex] = max(I);
    print('L(I Max)', L(maxIndex))

    % plot(L, dI_old);
    % title('dI_old');
    % grid on;
    
    ddI = p_0 .* log(lambda_0) ...
          + p_1 .* log(lambda_1) ... 
          - log(lambda_avg) ... 
          - 2 .* lambda_wave.^2./lambda_avg ...
          + 2;

    figure();
    plot(L, dI_new);
    title('dI_new');
    grid on;
    
    % figure();
    % plot(L, ddI);
    % title('ddI');
    % grid on;
    
    figure();
    plot(L, I);
    title('I');
    grid on;

    figure();
    plot(L, dI_new,'r');
    hold on
    plot(L, I,'g');
    title('All');
    grid on;

    % figure();
    % plot(L, lambda_0,'r');
    % hold on
    % plot(L, lambda_1,'g');
    % hold on
    % plot(L, lambda_avg,'b');
    % title('lambdas');
    % grid on;

end
     

function print(str, val)
    disp([str ':  ' num2str(val)]);
end