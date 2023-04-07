%% Problem 2

% I can't stand MatLab as a programming language!!!
% It's so awful!!! I'm using python instead!!!

%% Problem 3

h = 10.^(-1:-1:-10);
x = 0.7;
d = 1/x;  % Accurate derivative
dI = (log(x + h) - log(x))./h - d;
dIb = (log(x + h) - log(x))./(x+h-x) - d;
dII = (log(x+h) - log(x-h))./(2*h) - d;
dIIb = (log(x+h) - log(x-h))./((x+h)-(x-h)) - d;
dIII = (log(x-2*h) - 8*log(x-h) + 8*log(x+h) - log(x+2*h))./(12*h) - d;

figure(2)
loglog(h, abs(dI), "DisplayName", "f'_i", "Color", "red");
hold on
loglog(h, abs(dIb), ":", "DisplayName", "f'_i-bonus", "Color", "red");
loglog(h, abs(dII), "DisplayName", "f'_{ii}", "Color", "blue");
loglog(h, abs(dIIb), ":", "DisplayName", "f'_{ii}-bonus", "Color", "blue");
loglog(h, abs(dIII), "DisplayName", "f'_{iii}");
hold off
