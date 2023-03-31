%% 3.
% (ii)
figure(1)
% We evaluate 1/(1+x^2).
I_acc = pi / 2;
Ns = 2 .^ (2:10);
errors = arrayfun(@(N) CCquadrature(@(x) 1./(1+x.^2), N), Ns) - I_acc;
loglog(Ns, abs(errors), 'DisplayName', 'f_1')

hold on
% We evaluate exp(x)/cos(x).
I_acc = 2.9558446348045468559570113589764801300271657138917130259488856278;
errors = arrayfun(@(N) CCquadrature(@(x) exp(x)./cos(x), N), Ns) - I_acc;
loglog(Ns, abs(errors), 'DisplayName', 'f_2')
hold off

% (iii)
figure(2)
xs = (-1:0.02:1);
integral = arrayfun(@(b) CCquadratureAlt(@(x)1./(1+x.^2), 20, -1, b), xs);
true_integral = atan(xs) - atan(-1);
plot(xs, integral, 'DisplayName', 'numerical');
hold on
plot(xs, true_integral, 'DisplayName', 'accurate');
hold off

figure(3)
xs = (0:0.01:1);
integral = arrayfun(@(b) CCquadratureAlt(@(x)sin(20*x.^2), 30, 0, b), xs);
accurate = fresnels(xs * sqrt(40/pi)) * sqrt(pi/40);
plot(xs, abs(integral - accurate));


%% 4.
% (i)

figure(4)
Ns = 2.^(2:20);
hs = 1 ./ Ns;
integral_midpoint = arrayfun(@(N) midpoint(N), Ns);
integral_trapez = arrayfun(@(N) trapez(N), Ns);
integral_simpson = arrayfun(@(N) simpson(N), Ns);
loglog(hs, abs(integral_midpoint - pi), 'DisplayName', 'midpoint');
hold on
loglog(hs, abs(integral_trapez - pi), 'DisplayName', 'trapez');
loglog(hs, abs(integral_simpson - pi), 'DisplayName', 'simpson');
hold off

% (ii)
N = 10;
points = integral_midpoint(1:N);
figure(5)
semilogy(abs(points - pi), 'DisplayName', 'R_0');
hold on
for k = (1:N-1)
  % points have size (N-k+1) ~> (N-k)
  points = points(2:N-k+1) + (points(2:N-k+1) - points(1:N-k)) / (4^k - 1);
  semilogy(abs(points - pi), 'DisplayName', 'R_' + string(k));
end
hold off


%% Functions

function val = CCquadrature(f, N)
  samples = arrayfun(f, cos((0:2*N-1) / N * pi));
  coef = fft(samples);
  coef(1) = coef(1)/2;
  val = dot(coef(1:2:N), 2 ./ (1 - (0:2:N-1).^2)) / N;
end

function val = CCquadratureAlt(f, N, a, b)
  samples = arrayfun(f, cos((0:2*N-1) / N * pi));
  coef = fft(samples);
  coef(1) = coef(1)/2;
  val = dot(coef(1:2:N), iChebyshev((0:2:N-1), b) - iChebyshev((0:2:N-1), a)) / N;
end

function val = iChebyshev(k, x)
  val = 0.5 .* (cos((k+1) .* acos(x)) ./ (k+1) ...
            - cos((k-1) .* acos(x)) ./ (k-1));
end

function val = midpoint(N)
  pts = ((1:N) - 0.5)./N;
  val = sum(4./(1+pts.^2)) / N;
end

function val = trapez(N)
  pts = (1:N-1) ./ N;
  val = (sum(4./(1+pts.^2)) + 3) / N;
    % 3 = (f(0) + f(1)) / 2
end

function val = simpson(N)
  pts_maj = ((1:N) - 0.5)./N;
  pts_min = (1:N-1)./N;
  val = (sum(4./(1+pts_maj.^2)) * 2/3 + sum(4./(1+pts_min.^2)) / 3 + 1) / N;
end
