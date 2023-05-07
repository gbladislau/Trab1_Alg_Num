%%
%% (ESTA FUNÇÃO ESTAVA IMPLEMENTADA NA VERSÃO ANTERIOR DE OCTAVE)
%%
%% Troca de valores de escalares, vetores ou matrizes
%%
%% Input: Elementos a, b
%% Output: Elementos a, b com os valores trocados
%%
function [a,b] = swap(a,b)
    [lin, col] = size(a);
    if sum( [lin,col] ~= size(b) )
        error('Erro de dimensao de matriz e/ou vetor');
    end
    for i=1:lin
        for j=1:col
            s = a(i,j); a(i,j) = b(i,j); b(i,j) = s;
        end
    end
end


