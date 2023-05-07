%
% https://en.wikipedia.org/wiki/Probability_distribution#Absolutely_continuous_probability_distribution
% https://en.wikipedia.org/wiki/List_of_probability_distributions
%
%==============================================================
% Distribuições probabilísticas
% Autor: Thomas W. Rauber trauber@gmail.com , 2008 - 2022
%%==============================================================

clear all;
close all;
raiz = '/home/thomas/Nextcloud/algnum/soft/';
figdir = './';


addpath([raiz 'util'], [raiz 'integral'], [raiz 'raiz']);

is_octave = (exist('OCTAVE_VERSION','builtin')>1); % Octave ou Matlab
if is_octave
	available_graphics_toolkits();
	%graphics_toolkit gnuplot;
	graphics_toolkit fltk;
	%graphics_toolkit qt;

	pkg load symbolic;
	pkg load statistics;
	%
	% Se um pacote, por exemplo, o pacote 'symbolic' do Octave ainda
	% não estiver instalado, pegue aqui:
	% https://sourceforge.net/projects/octave/
	% instale com o comando dentro do Octave: pkg install PACOTE.tar.gz
	% e depois carregue com: pkg load PACOTE
end;


salvar_grafico = @(nome, formato) ...
    {fprintf('Salvando gráfico em %s no formato %s ...\n', nome, formato(3:end)); print(nome, formato);};


salvar_eps = @(figdir, nome) ...
    {...
        fnam = [figdir nome]; fnamsvg = [fnam '.svg'];...
        fprintf('Salvando gráfico em %s ...\n', fnamsvg);...
        print([fnam '.svg'], '-dsvg');...
        bashcatch = 'echo ''Não consigo converter de SVG para EPS. Inkscape instalado?''';
        cmd = ['inkscape --without-gui --export-eps=' fnam '.eps ' fnamsvg ' || ' bashcatch];...
        fprintf('Executando ''%s'' ...\n', cmd);...
        salvar_eps = system(cmd);...
    };
  
fontsize_title = 16;
fontsize_legend = 12;
fontsize_ax_label = 12;
linewidth_plot = 1.25;

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distribuição Normal (Gaussiana)                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://pt.wikipedia.org/wiki/Distribui%C3%A7%C3%A3o_normal

syms f(x) finv(y) F(x) a b offset mu sigma;

mu = sym(-2);
sigma = sqrt(sym(1)/sym(2));

f(x) = (1/(sigma*sqrt(sym(2)*pi)))*exp(-(sym(1)/sym(2))*((x-mu)/sigma)^2)
fstr = 'sym(1)/(sqrt(sym(1)/sym(2))*sqrt(sym(2)*pi))*exp(-(sym(1)/sym(2))*((x-sym(-sym(2)))/sqrt(sym(1)/sym(2)))^sym(2))';
F(x) = (sym(1) + erf((x-mu)/(sqrt(sym(2)*sigma)))) / sym(2)
Finv(y) = mu + sqrt(sym(2))*sigma * erfinv(sym(2)*y-sym(1))


percentil50 = Finv(sym(50)/sym(100))
fprintf('Percentil 50 = %20.10e\n', double(percentil50));
percentil99 = Finv(sym(99)/sym(100))
fprintf('Percentil 99 = %20.10e\n', double(percentil99));

%return

offset = sym(10) * sigma;
a = mu-offset;
b = percentil99;


I = F(b)-F(a);
II = double(I);
%trapeziosPlot( fstr, a, b, [1, 2, 3, 5, 10], 'lixo.eps' );
%simpsonPlot( fstr, a, b, [1, 2, 3, 5, 10], 'lixo.eps' );
%return

FDP_Normal = @(x, mu, sigma) 1./(sigma.*sqrt(2*pi)) .*exp(-0.5*((x-mu)./sigma).^2);
FDA_Normal = @(x, mu, sigma) 0.5*(1+erf((x-mu)./(sigma*sqrt(2))));
FDA_Normal_inv = @(y, mu, sigma) mu+sqrt(2)*sigma*erfinv(2*y-1);


mu = double(mu);
sigma = double(sigma);
a = double(a);
b = double(b);
y = 0.99;
x = FDA_Normal_inv(y, mu, sigma);
fprintf('Gaussiana: FDA(%.2f; %.2f, %.2f)=%.10e\n', x, mu, sigma, y);

fprintf('Integração numérica ...\n');

func = @(x) FDP_Normal(x, mu, sigma);	% Wrapper
%func = @(x) normpdf(x, mu, sigma);	% Wrapper
nmax = 10;	% Subdivisões
verbose = true;
verbose = false;
C = coefGaussLegendre( nmax+1 );
[T, A] = tabelaAbcissasPesosGaussLegendre( C );
% showGaussLegendre( C, T, A );

if 1
for n=1:nmax
    ITR = integralTrapeziosRepetidaFunc( func, a, b, n, verbose );
    errTR = II - ITR;
    errSR = NaN;
    ISR = NaN;
    if ~mod(n, 2)
        ISR = integralSimpsonRepetidaFunc(func, a, b, n, verbose );
        errSR = II - ISR;
    end
    % IGQ = quad(func, a, b);
    IGQ = integralGaussLegendreFunc( func, a, b, n, T, A, verbose );
    errGQ = II - IGQ;
    fprintf('n=%2d ITR=%15.6e ISR=%15.6e IGQ=%15.6e --- Erro T=%15.6e  Erro S=%15.6e  Erro QG=%15.6e\n',...
            n, ITR, ISR, IGQ, errTR, errSR, errGQ)
end
end

% Determinar x por bissecção para que FDA = 99%
func = @(x) -0.99 + FDA_Normal(x, mu, sigma);
a = -1;
b = 0;
maxiter = 10;
eps = 0.0;
posfalsa = false;
r = raizBisecPosFalsa( posfalsa, func, a, b, eps, maxiter );

%return



% Reproduzir desenho do Wikipedia

close all;
x = -5.0:0.01:5.0;
clf;
hold on;
leg = {};
% Pares (mu, sigma)
param = {
    {0.0, sqrt(0.2)},
    {0.0, 1},
    {0.0, sqrt(5)},
    {-2.0, sqrt(0.5)},
};
for i=1:length(param)
    mu = param{i}{1};
    sigma = param{i}{2};
    y = FDP_Normal(x, mu, sigma);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\mu=%3.1f, \\sigma=%3.1f', mu, sigma);
end
ylim([0 1.0]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('f(x)', 'fontsize', fontsize_ax_label);
title('Distribuição Normal (Gaussiana): FDP', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'northeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDP_Normal');
%input('...');
hold off;

fig = figure;
hold on;
leg = {};
for i=1:length(param)
    mu = param{i}{1};
    sigma = param{i}{2};
    y = log(FDP_Normal(x, mu, sigma));
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\mu=%3.1f, \\sigma=%3.1f', mu, sigma);
end
ylim([-20 1.0]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('F(x)', 'fontsize', fontsize_ax_label);
title('Logaritmo da Distribuição de Normal (Gaussiana): FDA', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'southeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDP_logNormal');
hold off;

fig = figure;
hold on;
leg = {};
for i=1:length(param)
    mu = param{i}{1};
    sigma = param{i}{2};
    y = FDA_Normal(x, mu, sigma);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\mu=%3.1f, \\sigma=%3.1f', mu, sigma);
end
ylim([0 1.0]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('F(x)', 'fontsize', fontsize_ax_label);
title('Distribuição de Normal (Gaussiana): FDA', 'fontsize', fontsize_title);
h = legend(leg);
set(h, 'interpreter', 'tex', 'Location', 'southeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDA_Normal');
hold off;

% Geração de números aleatórios
semente = 2021;
randn('seed', semente);
mu = -2;
sigma = 0.7;
n = 10000;
x = randn(1, n) * sigma + mu; % -2.6242  -1.2363  -2.2256 ... 
x(1:3)

muhat = sum(x)/n;
sigmahat = sqrt((1/(n-1))*sum((x-muhat).^2));

fprintf('Parâmetros:\nVerdadeiros=%.20e, %.20e\nEstimados=  %.20e, %.20e\n',...
    mu, sigma, muhat, sigmahat)
end
%%% END Distribuição Normal (Gaussiana)                    %%%%
%return

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distribuição exponencial                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://en.wikipedia.org/wiki/Weibull_distribution

FDP_exponencial = @(x, lambda) lambda.*exp(-lambda.*x);
FDA_exponencial = @(x, lambda) 1 - exp(-lambda.*x);
% Reproduzir desenho do Wikipedia

close all;
x = 0:0.01:5;
clf;
hold on;
leg = {};
% Parâmtro único: lambda
param = {
    {0.5},
    {1.0},
    {1.5},
};
for i=1:length(param)
    lambda = param{i}{1};
    y = FDP_exponencial(x, lambda);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\lambda=%3.1f', lambda);
end
ylim([0 2.5]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('f(x)', 'fontsize', fontsize_ax_label);
title('Distribuição exponencial: FDP', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'northeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDP_exponencial');
hold off;

fig = figure;
hold on;
leg = {};
for i=1:length(param)
    lambda = param{i}{1};
    y = FDA_exponencial(x, lambda);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\lambda=%3.1f', lambda);
end
ylim([0 1.0]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('F(x)', 'fontsize', fontsize_ax_label);
title('Distribuição exponencial: FDA', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'southeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDA_exponencial');
hold off;
end
%%% END Distribuição de exponencial                    %%%%


if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distribuição de Weibull                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://en.wikipedia.org/wiki/Weibull_distribution

FDP_Weibull = @(x, lambda, k) (k/lambda)*(x./lambda).^(k-1).*exp(-(x./lambda).^k);
FDA_Weibull = @(x, lambda, k) 1 - exp(-(x./lambda).^k);
% Reproduzir desenho do Wikipedia

close all;
x = 0:0.01:2.5;
clf;
hold on;
leg = {};
% Pares (lambda, k)
param = {
    {1.0, 0.5},
    {1.0, 1.0},
    {1.0, 1.5},
    {1.0, 5.0},
    {2.0, 1.0},
    {2.0, 5.0},
};
for i=1:length(param)
    lambda = param{i}{1};
    k = param{i}{2};
    y = FDP_Weibull(x, lambda, k);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\lambda=%3.1f, k=%3.1f', lambda, k);
end
ylim([0 2.5]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('f(x)', 'fontsize', fontsize_ax_label);
title('Distribuição de Weibull: FDP', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'northeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDP_Weibull');
hold off;

fig = figure;
hold on;
leg = {};
for i=1:length(param)
    lambda = param{i}{1};
    k = param{i}{2};
    y = FDA_Weibull(x, lambda, k);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\lambda=%3.1f, k=%3.1f', lambda, k);
end
ylim([0 1.0]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('F(x)', 'fontsize', fontsize_ax_label);
title('Distribuição de Weibull: FDA', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'southeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'FDA_Weibull');
hold off;
end
%%% END Distribuição de Weibull                    %%%%


if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distribuição de Pareto Generalizada                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===== Distribuição de Pareto Generalizada =========\n');

% Amostragem

% LOCATION = mu
mu = 0.0;
% SCALE
s = 1;
% SHAPE
ksi = 1;


% FDP
gpd = @(x, mu, s, ksi) (1.0/s) * (1 + ksi*(x-mu)/s).^(-1-1/ksi)

% FDA
gpdfda = @(x, mu, s, ksi)   ifelse(ksi ~= 0,...
                            1-(1+ksi*(x-mu)/s).^(-1/ksi),...
                            1-exp(-(x-mu)/s));


% fig = gcf;	% handle to current figure
% clf(fig);
% close(fig);
% close all;
fig = figure;

x = mu:0.01:5;
clf;
hold on;
leg = {};
% Pares (sigma, ksi)
param = {
    {1.0, 1.0},
    {1.0, 5.0},
    {1.0, 20.0},
    {2.0, 1.0},
    {2.0, 5.0},
    {2.0, 20.0},
};
for i=1:length(param)
    s = param{i}{1};
    ksi = param{i}{2};
    y = gpd(x, mu, s, ksi);
    plot(x, y, 'linewidth', linewidth_plot);;
    leg{end+1} = sprintf('\\mu=%3.1f, \\sigma=%3.1f, \\xi=%3.1f', mu, s, ksi);
end
ylim([0 1.0]);

xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('f(x)', 'fontsize', fontsize_ax_label);
title('Distribuição de Pareto Generalizada: FDP', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'northeast', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'GPDFDP');
hold off;

fig = figure;
hold on;
leg = {};
% Pares (sigma, ksi)
for i=1:length(param)
    s = param{i}{1};
    ksi = param{i}{2};
    y = gpdfda(x, mu, s, ksi);
    plot(x, y, 'linewidth', linewidth_plot);
    leg{end+1} = sprintf('\\mu=%3.1f, \\sigma=%3.1f, \\xi=%3.1f', mu, s, ksi);
end
ylim([0 1.0]);
xlabel('x', 'fontsize', fontsize_ax_label);
ylabel('F(x)', 'fontsize', fontsize_ax_label);
title('Distribuição de Pareto Generalizada: FDA', 'fontsize', fontsize_title);
h = legend(leg);
set (h, 'interpreter', 'tex', 'Location', 'northwest', 'fontsize', fontsize_legend);
salvar_eps(figdir, 'GPDFDA');
hold off;
end
%%% END Distribuição de Pareto Generalizada                    %%%%


