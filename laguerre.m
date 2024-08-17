%% Question 1

function F = laguerre(p, n, d)
    % Returns:
    %   F - Cell array of transfer functions representing Laguerre basis filters
    if p <= 0 || p >= 1
        error('Pole must be in the range (0, 1)');
    end
    % Define symbolic variable z for transfer functions
    z = tf('z', 1);

    % Initialize the Laguerre basis filters
    F = cell(1, n);

    % Calculate the Laguerre basis filters
    for i = 1:n
        F{i} = z^(-1-d) * (sqrt(1 - p^2) / (1 - p * z^(-1))) * ((z^(-1) - p) / (1 - p * z^(-1)))^(i - 1);
    end
end
