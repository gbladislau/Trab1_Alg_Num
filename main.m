pkg load symbolic
clc
addpath('./edo')
addpath('./util')
warning('off', 'all');
close all;
clear all;

function SolveLetraA()
    clear;
    fprintf('---------Solucao Letra A---------\n')
    % definicao das variavies
    syms y(x) x0 y0;

    % definindo o lado direito da EDO
    RHS = '(1 - (2* exp(x)* y ))/exp(x)'

    % condicoes iniciais
    x0 = 0;
    y0 = 1;

    [f, sol, PVIstr, yx, yxstr] = solveEDO( RHS, x0, y0 );

    %    f   - O lado direito f(x, y(x)) da ODE como funcção numérica
    %    sol - Solução simbólica do problema do valor inicial
    %    PVIstr - Descrição textual do PVI
    %    yx - Solução numérica do PVI
    %    yxstr - Descrição textual da função que é solução do PVI

    fprintf('\n')
    disp(PVIstr)
    fprintf('\n1 - Solucao Analítica do PVI: ')
    disp(yxstr)
    fprintf('\n2 - Solucao Numérica: ')
    disp(yx)
    fprintf('\n')
    n = 5.0;
    passo = 0.5;

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0,'A')

end

function SolveLetraB()
    clear;
    fprintf('---------Solucao Letra B---------\n')
    % definicao das variavies
    syms y(x) x0 y0;

    % definindo o lado direito da EDO

    RHS = "((sin(x)/x^2) - (3*y))/x"

    % condicoes iniciais
    x0 = pi;
    y0 = 1;

    [f, sol, PVIstr, yx, yxstr] = solveEDO( RHS, x0, y0 );

    %    f   - O lado direito f(x, y(x)) da ODE como funcção numérica
    %    sol - Solução simbólica do problema do valor inicial
    %    PVIstr - Descrição textual do PVI
    %    yx - Solução numérica do PVI
    %    yxstr - Descrição textual da função que é solução do PVI

    fprintf('\n')
    disp(PVIstr)
    fprintf('\n1 - Solucao Analítica do PVI: ')
    disp(yxstr)
    fprintf('\n2 - Solucao Numérica: ')
    disp(yx)

    fprintf('\n')
    n = 5.0;
    passo = 1.0;

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0,'B')
end

function SolveLetraC()
    clear;
    fprintf('---------Solucao Letra C---------\n')
    % definicao das variavies
    syms y(x) x0 y0;

    % definindo o lado direito da EDO


    RHS = "(power(cos(x),2) - (tan(x)*y))"


    % condicoes iniciais
    x0 = (pi/8);
    y0 = 1;

    [f, sol, PVIstr, yx, yxstr] = solveEDO( RHS, x0, y0 );

    %    f   - O lado direito f(x, y(x)) da ODE como funcção numérica
    %    sol - Solução simbólica do problema do valor inicial
    %    PVIstr - Descrição textual do PVI
    %    yx - Solução numérica do PVI
    %    yxstr - Descrição textual da função que é solução do PVI

    fprintf('\n')
    disp(PVIstr)
    fprintf('\n1 - Solucao Analítica do PVI: ')
    disp(yxstr)
    fprintf('\n2 - Solucao Numérica: ')
    disp(yx)

    fprintf('\n')
    n = 5.0;
    passo = (pi/16);

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0,'C')
end

function SolveLetraD()
    clear;
    fprintf('---------Solucao Letra D---------\n')
    % definicao das variavies
    syms y(x) x0 y0;

    % definindo o lado direito da EDO

    RHS = "(1 - (1/x) - (2*y))/x"

    % condicoes iniciais
    x0 = 1; %*
    y0 = 1; %*

    [f, sol, PVIstr, yx, yxstr] = solveEDO( RHS, x0, y0 );

    %    f   - O lado direito f(x, y(x)) da ODE como funcção numérica
    %    sol - Solução simbólica do problema do valor inicial
    %    PVIstr - Descrição textual do PVI
    %    yx - Solução numérica do PVI
    %    yxstr - Descrição textual da função que é solução do PVI

    fprintf('\n')
    disp(PVIstr)
    fprintf('\n1 - Solucao Analítica do PVI: ')
    disp(yxstr)
    fprintf('\n2 - Solucao Numérica: ')
    disp(yx)

    fprintf('\n')
    n = 5.0;
    passo = 0.1;

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0,'D')
end

function null = PlotaGraficoComSolucoes(f,yx, passo, n,x0,y0,nome)
    Y_Solucoes = [];
    x = x0;
    output = [];
    fprintf("3 - Discretizando Solucao")
    for i = 0 : n
        fprintf("Valor de x = %f // Valor de y = %f\n",x, yx(x))
        output = [output ; x, yx(x)];
        x = x+passo;
    end

    figure
    hold on;
    title1 = strcat("3 - Discretizacao da Solucao, Letra ",nome);
    title(title1)

    scatter(output(:,1),output(:,2), 20, 'k', 'filled')

    ax = [output(1,1) output(end,1)];
    x = ax(1) : passo/10 : ax(2);
    plot(x,yx(x),'r')

    legend("(x,y)","y(x)",'location','northeastoutside');
    epsfilename = strcat('Letra_',nome,'_fig_1');
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n\n', epsfilename );
    print( epsfilename ,'-depsc2');

    hold off;
    figure
    hold on;
    #MOSTRANDO ONDE estao os passos x e y
    title1 = strcat("5 - Solucao por metodo, Letra ", nome);
    title(title1)

    scatter(output(:,1),output(:,2), 20, 'k', 'filled')

    leg = {};
    leg{end+1} = sprintf('(x,y)');
    leg{end+1} = sprintf('y(x)');
    # mostando a funcao simbolica
    ax = [output(1,1) output(end,1)];
    x = ax(1) : passo/10 : ax(2);
    plot(x,yx(x),'b')

    Y_Solucoes = [Y_Solucoes output];

    fprintf("4 - Calculando valores por métodos e salvando para plotar\n\n")
    leg{end+1} = sprintf('Euler');
    [Euler_x Euler_y] = Euler(f, x0, y0, passo, n);
    plot(Euler_x,Euler_y,'r+-');
    Y_Solucoes = [Y_Solucoes Euler_y(:)];

    leg{end+1} = sprintf('Euler Melhorado');
    [EulerMe_x EulerMe_y] = EulerMelhorado(f, x0, y0, passo, n);
    plot(EulerMe_x,EulerMe_y,'g');
    Y_Solucoes = [Y_Solucoes EulerMe_y(:)];

    leg{end+1} = sprintf('Euler Modificado');
    [EulerMod_x EulerMod_y] = EulerModificado(f, x0, y0, passo, n);
    plot(EulerMod_x,EulerMod_y,'m--');
    Y_Solucoes = [Y_Solucoes EulerMod_y(:)];

    leg{end+1} = sprintf('Van der Houwen’s/Wray third-order');
    butcher.a = zeros(3,3);
    butcher.a(2,1) = 8/15;
    butcher.a(3,1) = 1/4; butcher.a(3,2) = 5/12;
    butcher.c = [0 8/15 2/3];
    butcher.b = [1/4 0 3/4];
    butcher.isEmbedded = false;
    [Van_T Van_Y Y_low] = RungeKutta(f, x0, y0, passo, n, butcher, 3);
    plot(Van_T,Van_Y,'cd-');
    Y_Solucoes = [Y_Solucoes Van_T(:)];

    leg{end+1} = sprintf('Ralston’s fourth-order method');
    butcher.a = zeros(4,4);
    butcher.a(2,1) = 0.4;
    butcher.a(3,1) = .29697761; butcher.a(3,2) = .15875964;
    butcher.a(4,1) =  .21810040; butcher.a(4,2) = -3.05096516; butcher.a(4,3) = 3.83286476;
    butcher.c = [0 .4 .45573725 1];
    butcher.b = [.17476028 -0.55148066 1.20553560 .17118478];
    butcher.isEmbedded = false;
    [Rals_T Rals_Y Y_low] = RungeKutta(f, x0, y0, passo, n, butcher, 4);
    plot(Rals_T ,Rals_Y,'rs-');
    Y_Solucoes = [Y_Solucoes Rals_T(:)];

    leg{end+1} = sprintf('Dormand_Prince RungeKutta');
    [DP_x DP_y Y_low] = RungeKutta_Dormand_Prince45(f, x0, y0, passo, n);
    plot(DP_x, DP_y,'bo-');
    Y_Solucoes = [Y_Solucoes DP_y(:)];

    leg{end+1} = sprintf('Dormand_Prince PassoFixo');
    [DPpf_x DPpf_y] = RungeKutta_Dormand_Prince_ode45(f, x0, y0, passo, n,true);
    plot(DPpf_x, DPpf_y,'md-');
    Y_Solucoes = [Y_Solucoes DPpf_y(:)];

    leg{end+1} = sprintf('Dormand_Prince PassoAdaptativo');
    [DPpa_x DPpa_y] = RungeKutta_Dormand_Prince_ode45(f, x0, y0, passo, n,false);
    plot(DPpa_x, DPpa_y,'rd-');

    legend(leg,'location','northeastoutside');
    epsfilename = strcat('Letra_',nome,'_fig_2');
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n\n', epsfilename );
    print( epsfilename ,'-depsc2');


    Erros = [];
    Erros = [Erros output(:,1)];
    for i = 2 : size(Y_Solucoes,2)
        Erros = [Erros abs(Y_Solucoes(:,2) - Y_Solucoes(:,i))];
    endfor


    %QUESTÃO 6:
    figure
    hold on;
    title1 = strcat("6 - Erro por metodo, Letra ", nome);
    title(title1)
    markers = {'r+-','g', 'm--', 'cd-', 'rs-', 'bo-', 'md-','rd-'};
    mk = 1;
    for i = 3:9
        semilogy(Euler_x(2:end),Erros(2:end,i),markers(mk))
        mk++;
    endfor

    ErroPasso_Adaptativo = abs(yx(DPpa_x) - DPpa_y);
    semilogy(DPpa_x,ErroPasso_Adaptativo,markers(mk))


    legend("Euler", "Euler Melhorado", "Euler Modificado",
           "Van der Houwen’s/Wray third-order", "Ralston’s fourth-order method",
           "Dormand_Prince RungeKutta", "Dormand_Prince PassoFixo","Dormand_Prince PassoAdaptativo",
           'location','northeastoutside');

    epsfilename = strcat('Letra_',nome,'_fig_3_erro');
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n\n', epsfilename );
    print( epsfilename ,'-depsc2');

    fprintf("7 e 8 - Tabela com valores das soluções e Erros dos metodos.\n")
    fprintf('%12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s | %12s\n',
    'x', 'Valor Exato',
    'Euler', 'Euler Mel.', 'Euler Mod.',
    'V d Houven/Wray', 'Ralston' ,
    'Dorm.-Pr45-Bu', 'ODE45 fixo', 'ODE45 adap');
    fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------\n');

    for i=1:length(Euler_x)
        fprintf('%12.7f | %12.7f | %12.7f | %12.7f | %12.7f | %12.7f    | %12.7f | %12.7f  | %12.7f | ----\n',
        Euler_x(i) , yx(Euler_x(i)) , Euler_y(i), EulerMe_y(i), EulerMod_y(i),Van_Y(i), Rals_Y(i), DP_y(i),
        DPpf_y(i));
    end

    nome_arq = cstrcat("letra_",nome,"_vetores.txt");
    save("-text",nome_arq, 'Y_Solucoes' );

    fprintf("ERROS\n")
    for i=1:length(Euler_x)
        fprintf('%12.7f | %12.7f | %12.7f | %12.7f | %12.7f | %12.7f    | %12.7f | %12.7f  | %12.7f | ----\n',
        Euler_x(i) , yx(Euler_x(i)) - yx(Euler_x(i)) , Erros(i,3) , Erros(i,4), Erros(i,5) ,Erros(i,6),Erros(i,7), Erros(i,8),
        Erros(i,9));
    end
    nome_arq = cstrcat("letra_",nome,'_erros.txt');
    vet_erros = [Euler_x(:) (yx(Euler_x(:))-yx(Euler_x(:)))(:) Erros(:,3:9)];
    save("-text", nome_arq, 'vet_erros');

end

function null = SolveQuestao31()

    % C = Qin - Qout, sendo uma constante
    %SOLUÇÃO NUMÉRICA DO PVI

    syms y(t) C t0 v0 Qin Qout;

    DE = diff(y, t) == C;
    cond = y(t0) == v0;
    fprintf("1. Solução Analítica V(t)\n")
    V = dsolve(DE,cond)

    %SOLUÇÃO SIMBÓLICA DO PVI
    fprintf("2. Função numérica\n")
    yt = matlabFunction(V)
    %C, t, t0, v0

    %3 CASOS:

    t0 = 0;
    v0 = 2000;
    Vmax = 5000;
    Vmin = 0;

    %Cenário de esvaziamento: Qin = 45, Qout = 50
    X = [0 400];
    fprintf("3. Tempo para esvaziamento: 400\n")
    Qin = 45;
    Qout = 50;
    t = 0:0.1:400;

    sol = yt((Qin-Qout), t, t0, v0);

    figure
    hold on;
    axis([0,400,0,5000]);
    line(X,[v0 v0], "linestyle", "--", "color", "g");
    line(X,[Vmax Vmax], "linestyle", "--", "color", "r");
    plot(t, sol);

    title("Evolucao temporal do volume no tanque");
    xlabel("t [min]");
    ylabel("V(t) [L]");

    a = 400;
    plot([a a], [0 Vmax])

    legend("v0 = 2000.00 L", "Vmax = 5000.00 L","V (t)", "Vazamento completo do tanque",'location','northeastoutside');
    epsfilename = '3.1_esvaziamento.eps';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    print( epsfilename ,'-depsc2');



    hold off;

    %Cenário de transbordamento: Qin = 50, Qout = 45
    Qin = 50;
    Qout = 45;
    fprintf("3. Tempo para transbordamento: 600\n")

    figure
    hold on;
    t = 0:0.1:600;
    sol = yt(Qin-Qout, t, t0, v0);

    axis([0,600, 0, 5000]);
    X = [0 600];
    line(X,[v0 v0], "linestyle", "--", "color", "g");
    line(X,[Vmax Vmax], "linestyle", "--", "color", "r");
    plot(t, sol);

    title("Evolucao temporal do volume no tanque")
    xlabel("t [min")
    ylabel("V(t) [L]")
    a = 600;
    plot([a a], [0 Vmax])
    legend("v0 = 2000.00 L", "Vmax = 5000.00 L","V (t)", "Transbordagem do tanque",'location','northeastoutside');
    epsfilename = '3.1_transbordamento.eps';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    print( epsfilename ,'-depsc2');
    hold off;

    %Cenário de constância no volume: Qin = 50, Qout = 50


    figure
    hold on;
    Qin = 50;
    Qout = 50;

    t = 0:0.1:500;

    sol = yt(Qin-Qout, t, t0, v0);
    axis([0,500, 0, 5000]);
    X = [0 500];
    line(X,[v0 v0], "linestyle", "--", "color", "g");
    line(X,[Vmax Vmax], "linestyle", "--", "color", "r");

    plot(t, sol);

    title("Evolucao temporal do volume no tanque");
    xlabel("t [min]");
    ylabel("V(t) [L]");
    legend("v0 = 2000.00 L", "Vmax = 5000.00 L","V (t)",'location','northeastoutside');
    epsfilename = '3.1_volumecte.eps';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    print( epsfilename ,'-depsc2');
    hold off;
end




function null = Solve32_Qin_diff_Qout()


    syms cin c(t) c0 t0 v0 qin qout;

    v(t) = (qin-qout)*(t-t0) + v0;


    DE = diff(c, t) == (qin * (cin - c(t)))/v(t);
    cond = c(0) == c0;
    fprintf("1. Calcule analiticamente a solução c(t) do PVI da concentração\n")
    c(t) = dsolve(DE, cond)

    syms m(t);
    fprintf("2. Defina a função do aditivo m(t)\n")
    m(t)= v(t) * c(t)

    fprintf("3. Converta as funções c(t) e m(t) em funções numéricas pela função auxiliar matlabFunction\n")
    cNum = matlabFunction(c(t))
    mNum = matlabFunction(m(t))
    vNum = matlabFunction(v(t));

    v0 = 2000;
    c0 = 0.05;
    cin = 2;
    vMax = 5000;
    t0 = 0;

    fprintf("4. Caso de esvaziamento: Qin = 40, Qout = 45\n")

    %Caso de esvaziamento

    figure;
    subplot(2,1, 1)
    hold on;

    title("Caso de esvaziamento - Qin = 40, Qout = 45\n\nEvolucao temporal da concentracao para o caso vazamento")
    xlabel("t [min]")
    ylabel("m(t) [kg] V(t) [L]")
    qin = 40;
    qout = 45;

    t = 0:0.1:400;
    x = [0 400];
    cVet = cNum(c0, cin, qin, qout, t, t0, v0);


    plot(t, cVet)
    line(x, [c0 c0], "linestyle", "--", "color", "g")
    line(x, [cin cin], "linestyle", "--", "color", "r")



    legend("c(t)", "cin = 2.00 kg/L", "c0 = 0.05 kg/L",'location','northeastoutside')
    hold off;

    subplot(2,1,2)
    hold on;
    title("Evolucao temporal do material e volume do tanque para o caso vazamento")
    ylabel("c(t) [Kg/L]")
    xlabel("t [min]")

    mVet = mNum(c0, cin, qin, qout, t, t0, v0);
    axis([0,400, 0, 5000])
    plot(t, mVet)

    line(x, [v0 v0], "linestyle", "--", "color", "g")
    line(x, [vMax vMax], "linestyle", "--", "color", "r")

    a = 400;
    plot([a a], [0 vMax])

    vVet = vNum(qin, qout, t, t0, v0);
    plot(t, vVet)

    legend("m(t)", "v0 = 2000.00 L","Vmax = 5000.00, L", "Vazamento Completo do tanque","V(t)",'location','northeastoutside')

    epsfilename = '3.2_esvaziamento.eps';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    print( epsfilename ,'-depsc2');
    hold off;

    fprintf("4. Caso de transbordamento: Qin = 45, Qout = 40\n")


    %Caso de transbordamento

    figure;
    subplot(2,1,1)
    hold on;

    qin = 45;
    qout = 40;

    t = 0:0.1:600;
    x = [0 600];

    cVet = cNum(c0, cin, qin, qout, t, t0, v0);

    title("Caso de transbordamento - Qin = 45, Qout = 40\n\nEvoluo temporal da concentracao para o caso transbordamento\n\n")
    ylabel("c(t) [Kg/L]")
    xlabel("t [min]")
    plot(t, cVet)
    line(x, [c0 c0], "linestyle", "--", "color", "g")
    line(x, [cin cin], "linestyle", "--", "color", "r")

    legend("c(t)", "c0 = 0.05 kg/L", "cin = 2.00 kg/L",'location','northeastoutside')
    hold off;

    subplot(2,1,2)

    hold on;
    title("Evolucao temporal do material e volume do tanque para o caso transbordamento\n\n")
    ylabel("c(t) [Kg/L]")
    xlabel("t [min]")
    mVet = mNum(c0, cin, qin, qout, t, t0, v0);
    axis([0,600])
    plot(t, mVet)
    line(x, [v0 v0], "linestyle", "--", "color", "g")
    line(x, [vMax vMax], "linestyle", "--", "color", "r")

    vVet = vNum(qin, qout, t, t0, v0);
    plot(t, vVet)

    a = 600;
    plot([a a], [0 vMax])

    legend("m(t)", "v0 = 2000.00 L","Vmax = 5000.00 L",  "V(t)","Transbordagem do tanque",'location','northeastoutside')
    epsfilename = '3.2_transbordamento.eps';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    print( epsfilename ,'-depsc2');

end

function null = Solve32_Qin_equals_Qout()

    syms cin c(t) c0 t0 v0 qin qout;

    v(t) = v0;


    DE = diff(c, t) == (qin * (cin - c(t)))/v(t);
    cond = c(0) == c0;
    c(t) = dsolve(DE, cond);

    syms m(t);
    m(t)= v(t) * c(t);

    cNum = matlabFunction(c(t));
    mNum = matlabFunction(m(t));
    vNum = matlabFunction(v(t));

    v0 = 2000;
    c0 = 0.05;
    cin = 2;
    vMax = 5000;

    fprintf("4. Caso de constância: Qin = 45, Qout = 45\n")


    %CASO DE CONSTÂNCIA

    figure;
    subplot(2,1,1)
    hold on;



    qin = 45;
    qout = 45;

    t = 0:0.1:500;
    x = [0 500];

    cVet = cNum(c0, cin, qin, t, v0);

    title("Caso de constancia - Qin = 45, Qout = 45\n\nEvolucao temporal da concentracao para o caso constante\n\n")
    ylabel("c(t) [Kg/L]")
    xlabel("t [min]")
    plot(t, cVet)
    line(x, [c0 c0], "linestyle", "--", "color", "g")
    line(x, [cin cin], "linestyle", "--", "color", "r")

    legend("c(t)", "c0 = 0.05 kg/L", "cin = 2.00 kg/L",'location','northeastoutside')
    hold off;

    subplot(2,1,2)

    hold on;
    title("Evolucao temporal do material e volume do tanque para o caso constante\n\n")
    ylabel("c(t) [Kg/L]")
    xlabel("t [min]")

    mVet = mNum(c0, cin, qin, t, v0);

    axis([0,500])
    plot(t, mVet)
    line(x, [v0 v0], "linestyle", "--", "color", "g")
    line(x, [vMax vMax], "linestyle", "--", "color", "r")

    vVet = vNum(v0);
    %Como vVet é um vetor constante, temos que:
    line(x, [v0, v0], "linestyle", "-", "color", "k")


    legend("m(t)", "v0 = 2000.00 L","Vmax = 5000.00 L",  "V(t)",'location','northeastoutside')
    epsfilename = '3.2_volumecte.eps';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    print( epsfilename ,'-depsc2');
end

