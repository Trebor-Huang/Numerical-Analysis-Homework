%% Problem 2
hs = 2 .^ (-1:-1:-60);
ds = (exp(1+hs) - exp(1-hs)) ./ (2 * hs);
d = exp(1);

figure(1);
loglog(hs, abs(ds - d));

%% Problem 4
pts = (-1:0.01:1);  % Evaluation points

f = @(x) 1/(1+25*x^2);  % Hat
% f = @(x) exp(-x^2);  % Gauss

figure(2);
hold on
for N = [10 25 50 80 100]
    nodes = (0:N) * 2 / N - 1;  % Equidistant
    ylim([-2 2]);  % Prevent extreme errors to drown everything else
    % nodes = cos(pi * ((0:N) * 2 + 1) / (N + 1) / 2);  % Chebyshev
    
    plot(pts, arrayfun(@(x) interpolate(f, nodes, x) - f(x), pts), ...
        "DisplayName", "N = " + string(N))
end
hold off

function val = lagrange(Xs, Ys, x)
    n = length(Xs);
    val = 0;
    for i = 1:n
        factors = (x - Xs) ./ (Xs(i) - Xs);
        factors(i) = [];
        val = val + Ys(i) * prod(factors);
    end
end

function val = interpolate(f, Xs, x)
    val = lagrange(Xs, arrayfun(f, Xs), x);
end
