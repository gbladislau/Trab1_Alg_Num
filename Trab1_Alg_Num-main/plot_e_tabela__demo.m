

% Definição de duas funções: seno e cosseno
f1 = @(x) sin(x);
f2 = @(x) cos(x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       T A B E L A      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Eixo x como vetor de valores discretos entre -pi e pi com intervalo de 1/pi
% A divisão é mais grossa, comparado com o gráfico. Então temos menos pontos
x = -pi:1/pi:pi;


% Gerando uma tabela

fprintf('%10s | %10s | %10s\n', 'x', 'SEN', 'COS' );
fprintf('---------------------------------------------------\n');
for i=1:length(x)
	fprintf('%10.2f | %10.2f | %10.2f\n', x(i), f1(x(i)), f2(x(i)) );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       G R Á F I C O    %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Eixo x como vetor de valores discretos entre -pi e pi com intervalo de 0.01
% A divisão é mais fina para evitar curvas não suaves
xx = -pi:0.01:pi;

clf;  % Limpa gráfico
legenda = {};   % Inicializa célula que contém os rótulos dos objetos gráficos
hold on

plot(xx, f1(xx), '-r', 'linewidth', 3);
legenda{end+1} = 'sen(x)';    % Insere rótulo no final das células dos rótulos

plot(xx, f2(xx), '-g', 'linewidth', 3);
legenda{end+1} = 'cos(x)';    % Insere rótulo no final das células dos rótulos

xlabel('x');
ylabel('y');
title('Funcoes trigonometricas')
legend(legenda, 'location', 'northwest');
hold off
shg;  % Mostre gráfico
