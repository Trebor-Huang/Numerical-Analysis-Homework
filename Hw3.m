global Xs n xs;
% Points of evaluation
Xs = (-1 : 0.001 : 1);

Ns = (10:10:100);

figure(1)
semilogy(Ns, differences(@(x) abs(x)), 'DisplayName', 'f_1');

figure(2)
semilogy(Ns, differences(@(x) x*abs(x)), 'DisplayName', 'f_2');
hold on
semilogy(Ns, differences(@(x) exp(-3*x^2)), 'DisplayName', 'f_3');
hold off

figure(3)
semilogy(Ns, differences(@(x) tanh(50*pi*x)), 'DisplayName', 'Actual');

M = tanh(50*pi);
rho = sqrt(1.0001) + 0.01;
invrm1 = (sqrt(10001) + 99)/2;  % inverse rho - 1 for accuracy
predicted = 4 * M * invrm1 * rho .^ (-Ns);
hold on
semilogy(Ns, predicted, 'DisplayName', 'Predicted');
hold off

function setup(new_n)
  % Setting the number of nodes to precompute weights
  global n xs;
  n = new_n;
  xs = cos((0:n-1) * pi / (n-1));
end

function y = barycentric(ys, x)
  global n xs;
  dx = (-1) .^ (0:n-1) ./ (x - xs);
  dx([1;n]) = dx([1;n]) ./ 2;
  y = sum(ys .* dx) / sum(dx);
end

function interpolate(f)
  % Plots the difference
  global xs Xs;
  ys = arrayfun(f, xs);
  plot(Xs, arrayfun(@(x) barycentric(ys, x) - f(x), Xs));
end

function d = difference(f)
  % Only computes the maximal difference
  global xs Xs;
  ys = arrayfun(f, xs);
  d = max(abs(arrayfun(@(x) barycentric(ys, x) - f(x), Xs)));
end

function ds = differences(f)
  % Collects data for different setups
  ds = zeros(1,10);
  for k=(1:10)
    setup(k * 10);
    ds(k) = difference(f);
  end
end
