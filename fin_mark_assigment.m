% TPPE29 - Financial Markets and Instruments, Linkoping University.
% Submission task: Option pricing with binomial method


% Task 1 - Value an European call option without dividend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
r = 0.0278;                     % risk free interest rate
T1 = datetime(2022,11,28);      % Start date
T2 = datetime(2023,7,21);       % End date
T = days252bus(T1,T2);          % Business days
S0 = 62;                         % Underlying Value
K = 55;                         % Strike price
sigma = 0.2;                    % Spread
fiscal_periods = 8;             % Number of periods.

fprintf("\nTask 1:")

CalculateCallOption(T, sigma, fiscal_periods, S0, K, r, 0, 0, 0, "C-EU")


% Task 2 - Value an European call option on OMSX30 using Black-Scholes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
r = 0.0278;                     % risk free interest rate
T1 = datetime(2022,11,28);      % Start date
T2 = datetime(2023,3,17);       % End date
T = days252bus(T1,T2);          % Business days 
S0 = 2096.462;                  % Spot price
K = 2100;                       % Strike price
C = 95.50;                      % Ask
fiscal_periods = 10;            % Number of periods.

% Task 2a
disp("Task 2a:")

% This could be done in a more efficient way
for sigmaOpt = 0:0.00001:1
    c = BlackScholes(r, sigmaOpt, T, S0, K);
    if (C == round(c,2))
        disp("Implied Volatility: " + sigmaOpt)
        break;
    end
end

bsCallPrice = BlackScholes(r, sigmaOpt, T, S0, K);
disp("Call price (BS): " + bsCallPrice)
disp(" ") % new line

% Task 2b
disp("Task 2b:")
CalculateCallOption(T, sigmaOpt, fiscal_periods, S0, K, r, 0, 0, 0, "C-EU")

% Task 2c
disp("Task 2c:")
optionPrices = zeros(1,195); % Instantiate empty vector
convergingValue = 0;
for fiscal_periods = 5:200
    optionTree = CalculateCallOption(T, sigmaOpt, fiscal_periods, S0, K,...
        r, 0, 0, 0, "C-EU");
    optionPrice = optionTree(1,1);
    optionPrices(1,fiscal_periods) = optionPrice;
    if bsCallPrice * 0.995 <= optionPrice && optionPrice...
            <= bsCallPrice * 1.005 && convergingValue == 0
        convergingValue = fiscal_periods;
    end
end
plot(optionPrices);
disp("Converging Value: " + convergingValue)

% Task 3 - Value an European or American call option with dividend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
r = 0.0278;                     % risk free interest rate
T1 = datetime(2022,11,28);      % Start date
T2 = datetime(2023,7,21);       % End date
T2div = datetime(2023,3,21);    % Dividend end date
T = days252bus(T1,T2);          % Business days
Tdiv = days252bus(T1, T2div);   % Divident time
S0 = 62;                        % Spot price
K = 55;                         % Strike price
sigma = 0.2;                    % Spread
DIV = 8;                        % Divident
fiscal_periods = 8;             % Number of periods.

tdiv = Tdiv / 252;        
Sstar = S0 - DIV * exp(-r*tdiv);
fprintf("\nTask 3: \n")
% Prints bintree for European Call option
CalculateCallOption(T, sigma, fiscal_periods, Sstar, K, r, DIV,...
    tdiv, 0, "C-EU")
% Prints bintree for American Call option
CalculateCallOption(T, sigma, fiscal_periods, Sstar, K, r, DIV,...
    tdiv, 0, "C-AM")


% Task 4 - Value an American call option for SAAB AB.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
r = 0.0278;                     % risk free interest rate
T1 = datetime(2022,11,28);      % Start date
T2 = datetime(2023,9,15);       % End date
T2div = datetime(2023,4,13);    % Dividend end date
T = days252bus(T1,T2);          % Business days
Tdiv = days252bus(T1, T2div);   % Divident time
S0 = 385.90;                    % Spot price
K = 430;                        % Strike price
sigma = 0.3593;                 % Spread
DIV = 4.90;                     % Divident
fiscal_periods = 200;           % Number of periods.

tdiv = Tdiv / 252;        
Sstar = S0 - DIV * exp(-r*tdiv);
disp("Task 4:")
optTree = CalculateCallOption(T, sigma, fiscal_periods, Sstar, K,...
    r, DIV, tdiv, 0, "C-AM");
disp("Option Price: " + optTree(1,1));
disp(" ") % newline


% Task 5 - Value an Up-and-In barrier Option.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
r = 0.0278;                     % risk free interest rate
T1 = datetime(2022,11,28);      % Start date
T2 = datetime(2023,3,17);       % End date
T = days252bus(T1,T2);          % Business days 
S0 = 2096.462;                  % Spot price
K = 2100;                       % Strike price
fiscal_periods = 10;            % Number of periods.
sigma = sigmaOpt;               % Spread (from task 2)

disp("Task 5:")

H = S0 * 1.05; % Barrier

% Using black & scholes
upAndOut = UpAndOut(S0, K, r, T, sigma, H); % Calc up and out price with BS
plainVanilla = BlackScholes(r, sigma, T, S0, K); % price of plain option


upAndIn = plainVanilla - upAndOut;
disp("Up and in model using Black Scholes")
disp("Plain vanilla price " + plainVanilla);
disp("Up and in option price " + upAndIn);
disp(" ") % new line

% Using binomial option tree
upAndOut = CalculateCallOption(T, sigma, fiscal_periods, S0, K, r, 0, 0,...
    H, "C-EU");
plainVanilla = CalculateCallOption(T, sigma, fiscal_periods, S0, K, r,...
    0, 0, 0, "C-EU");
disp("Up and in model using binomial Tree")
upAndIn = plainVanilla(1,1) - upAndOut(1,1);
disp("Plain vanilla price " + plainVanilla(1,1));
disp("Up and in option price " + upAndIn)


% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Returns Binomial option tree for Amarican or European Call options
function optionBinTree = CalculateCallOption(T, sigma, ...
    fiscal_periods, S, K, r, DIV, tdiv, H, type)
    
    % Boolean to handle Up and In in task 5
    upOut = false;
    if H ~= 0
        upOut = true;
    end
    
    deltat = T / (252 * fiscal_periods);
    u = exp(sigma * sqrt(deltat));              % Up value underlying
    d = exp(-sigma * sqrt(deltat));             % Down value underlying
    q = (exp(r * deltat) - d) / (u - d);
    
    % Binomial stock tree for underlying stock.
    underlyingTree = zeros(fiscal_periods + 1); % instanciate empty matrix
    for i = 0:fiscal_periods
        for k = 0:i
            underlyingTree(k+1, i+1) = S * u^(i-k) * d^(k);
        end
    end
    
    optionTree = zeros(fiscal_periods + 1);     % instanciate empty matrix
    
    if type == "C-EU"
    
        % Binomial tree for European Call option (with or without dividend)
        for i = fiscal_periods:-1:0
            for k = 0:i
                if i == fiscal_periods
                    if underlyingTree(k+1, i+1) > H && upOut
                        optionTree(k+1, i+1) = 0;
                    else
                        optionTree(k+1, i+1) = ...
                            max(underlyingTree(k + 1,i + 1) - K, 0);
                    end
                else
                    if underlyingTree(k+1, i+1) > H && upOut
                        optionTree(k+1, i+1) = 0;
                    else
                        optionTree(k+1, i+1) = exp(-r*deltat) *...
                            (q * optionTree(k+1, i+2)...
                            + (1-q)*optionTree(k+2,i+2));
                    end
                end
            end
        end
        optionBinTree = optionTree;
    end

    if type == "C-AM"
        % Binomial tree for American Call option with dividend
        for i = fiscal_periods:-1:0
            for k = 0:i
                if i == fiscal_periods
                    optionTree(k+1, i+1) = ...
                        max(underlyingTree(k + 1,i + 1) - K, 0);
                else
                    if floor(tdiv/deltat) == i
                        optionTree(k+1, i+1) = ...
                            max(underlyingTree(k + 1,i + 1) + DIV - K,...
                            exp(-r*deltat) * (q * optionTree(k+1, i+2)...
                            + (1-q)*optionTree(k+2,i+2)));
                    else
                    optionTree(k+1, i+1) = exp(-r*deltat)...
                        * (q * optionTree(k+1, i+2)...
                        + (1-q)*optionTree(k+2,i+2));
                    end
                end
            end
        end
        optionBinTree = optionTree;
    end


end

% Returns option price based on Black & Scholes
function optionPrice = BlackScholes(r, sigma, T, S, K)
    T = T / 252; % Time in year
    d1 = (log(S/K) + (r + (sigma^2)/2) * T) / (sigma * sqrt(T));        
    d2 = d1 - sigma * sqrt(T);
    optionPrice = S * normcdf(d1) - K * exp(-r*T) * normcdf(d2);
end

% Returns upAndOut option price using Black & Scholes.
function price = UpAndOut(S, K, r, T, sigma, H)
    
    %Calculating call price using Black & Scholes.
    callPrice = BlackScholes(r, sigma, T, S, K);
    % Prob that the option will be knocked out
    prob = normcdf((log(H/S) + (r + sigma^2/2)*T)/(sigma*sqrt(T)));
    % Returns option price
    price = callPrice * (1 - prob);
end
