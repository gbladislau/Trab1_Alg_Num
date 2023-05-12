pkg load symbolic
clc
addpath('./edo')
addpath('./util')

function null = SolveLetraA()
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
    fprintf('\nSolucao do PVI: ')
    disp(yxstr)
    fprintf('\n')
    n = 5.0;
    passo = 0.5;

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0)

end

function null = SolveLetraB()
    clear;
    fprintf('---------Solucao Letra B---------\n')
    % definicao das variavies
    syms y(x) x0 y0;

    % definindo o lado direito da EDO

    %%%%% -CONSERTAR- %%%%%
    RHS = "((sin(x)/x^2) - (3*y))/x" 
    %%%%% ----------- %%%%%

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
    fprintf('\nSolucao do PVI: ')
    disp(yxstr)
    disp(sol)
    fprintf('\n')
    n = 5.0; 
    passo = 1.0; 

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0)
end

function null = SolveLetraC()
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
    fprintf('\nSolucao do PVI: ')
    disp(yxstr)
    disp(sol)
    fprintf('\n')
    n = 5.0;
    passo = (pi/16); 

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0)
end

function null = SolveLetraD()
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
    fprintf('\nSolucao do PVI: ')
    disp(yxstr)
    disp(sol)
    fprintf('\n')
    n = 5.0; 
    passo = 0.1; 

    PlotaGraficoComSolucoes(f,yx,passo,n,x0,y0)
end

function null = PlotaGraficoComSolucoes(f,yx, passo, n,x0,y0)
    Y_Solucoes = []
    x = passo;
    output = [];
    output = [output ; x0, y0];
    fprintf("Valor de x = %f // Valor de y = %f\n",x0,y0)
    for i = 0 : n-1
        fprintf("Valor de x = %f // Valor de y = %f\n",x, yx(x))
        output = [output ; x, yx(x)];
        x = x+passo;
    end

    figure(1);
    hold on;
    #MOSTRANDO ONDE estao os passos x e y

    scatter(output(:,1),output(:,2), 20, 'k', 'filled')

    leg = {};
    leg{end+1} = sprintf('(x,y)');
    leg{end+1} = sprintf('y(x)');
    # mostando a funcao simbolica
    ax = [output(1,1) output(end,1)];
    x = ax(1) : 0.01 : ax(2);
    plot(x,yx(x),'b')

    Y_Solucoes = [Y_Solucoes output];

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
    butcher.a(3,1) = 1/4; butcher.a(3,2) = 5/12
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

    set(legend(leg),'fontsize',18);
    epsfilename = 'Letra A';
    fprintf('Gerando grafico vetorial em arquivo EPS ''%s''...\n', epsfilename );
    %print(epsfilename, '-depsc2');


    Erros = [];
    Erros = [Erros output(:,1)];
    for i = 2 : size(Y_Solucoes,2)
        Erros = [Erros abs(Y_Solucoes(:,2) - Y_Solucoes(:,i))];
    endfor


    %QUESTÃO 6:

    figure (2);
    hold on;
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
           "Dormand_Prince RungeKutta", "Dormand_Prince PassoFixo","Dormand_Prince PassoAdaptativo")

    fprintf('%s | %s | %s | %s | %s | %s | %s | %s | %s | %s\n', 'x', 'Valor Exato', "Euler", "Euler Melhorado", "Euler Modificado",
           "Van der Houwen’s/Wray", "Ralston",
           "DP RungeKutta", "DP PassoFixo","DP PassoAdaptativo");
    fprintf('-------------------------------------------------------------------------------------------\n');

    for i=1:length(Euler_x)
        fprintf('%10.2f | %10.2f | %10.2f\n', x(i), Euler_x(i) , yx(Euler_x(i)) , Euler_y(i));
    end

end
SolveLetraD()