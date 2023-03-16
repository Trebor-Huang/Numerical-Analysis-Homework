%% Problem 1
format long
NUM_MAX = 20;
Jpi_accurate = besselj(0:NUM_MAX, pi);

% (i)
Jpi = Jpi_accurate(1:2);  % Initial values
for n = 2:NUM_MAX
    % Fills in the rest via recurrence
    Jpi(n+1) = 2*(n-1)*Jpi(n) / pi - Jpi(n-1);
end

figure(1);
semilogy(abs(Jpi_accurate - Jpi));
% The accuracy deteriorates quickly, going from zero up to 1e-3.


% (ii)
Jpi = [zeros(1,NUM_MAX), 1, 0];
% Since we know a total, and the recurrence is linear,
% we can be off by a constant proportion, and correct
% for it at the very end. Now we can do the recurrence
% backwards, starting with [..., c, 0], where c is an
% arbitrary constant. (This results in very large numbers
% near the beginning.)
for n = NUM_MAX:-1:1
    % Backwards recurrence
    Jpi(n) = 2*n*Jpi(n+1) / pi - Jpi(n+2);
end
Jpi = Jpi / (Jpi(1) + 2*sum(Jpi(3:2:NUM_MAX+1)));
Jpi = Jpi(1:NUM_MAX+1);

figure(2);
plot(Jpi_accurate - Jpi);
% The error is now stably around 1e-16.

%% Problem 2

I1_accurate = sqrt(pi);
I2_accurate = sqrt(pi) * erf(2);

Ns = (2:50); fun = @(x) exp(-x.^2);
% Here we cut off the infinite integral at sqrt(N)
I1_values = arrayfun(@(N) trapz(fun, -sqrt(N), sqrt(N), N), Ns);
I2_values = arrayfun(@(N) trapz(fun, -2, 2, N), Ns);

I1_errors = abs(I1_values - I1_accurate);
I2_errors = abs(I2_values - I2_accurate);
% We see the effect of cutoff error
I1_cutoff = sqrt(pi) * (1 - erf(sqrt(Ns)));

figure(3);
semilogy(Ns, I1_errors, Ns, I2_errors, Ns, I1_cutoff);
legend("I_1-errors", "I_2-errors", "I_1-cutoff");

function int = trapz(fun, left, right, N)
    h = (right-left) / N;
    pts = h * (0:N) + left;
    fpts = fun(pts);
    int = h/2 * sum(fpts(1:N) + fpts(2:N+1));
end
