addpath("./edo")
function null = SolveLetraB()
    fprintf('---------Solucao Letra B---------\n')
    % definicao das variavies
    syms y(x) x0 y0;

    % definindo o lado direito da EDO
    RHS = "(sin(x)/x^3) - (3*y/x)"

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
    disp(sol)
    fprintf('\n')
    n = 5.0;
    passo = 1.0;

    x = passo;
    output =[];
    output = [output ; x0, y0];
    fprintf("Valor de x = %f // Valor de y = %f\n",x0,y0)
    for i = 0 : n-1
        fprintf("Valor de x = %f // Valor de y = %f\n",x, yx(x))
        output = [output ; x, yx(x)];
        x = x+passo;
    end

    plot(output(:,1),output(:,2))

end



SolveLetraB()









