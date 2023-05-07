%%
%% Convers�o de um n�mero decimal para fra��o
%%
%% Input: N�mero decimal 
%% Output: vetor de solu��o, e indica��o de sucesso
%%
% Usage: [n, d, success] = dec2frac( 4.0/3 )
function [numerator, denominator, success] = dec2frac( decimal )
    if ~isscalar(decimal)
        error('Argumento tem que ser escalar (um valor somente).');
    end
    accuracy = 1E-5;    % precis�o de convers�o
    maxiter = 1000;     % n�mero m�ximo de itera��es
    numerator = 0;
    denominator = 1;
    negative = sign( decimal );
    fraction = 0.0;
    success = true;

    iter = 0;
    while abs( fraction - decimal ) > accuracy && iter < maxiter
        if abs(fraction) > abs(decimal)
            denominator = denominator + 1;
        else
            numerator = numerator + negative;
        end
        fraction = numerator / denominator;
        iter = iter + 1;
    end
    if iter == maxiter
        success = false;
        %fprintf('dec2frac> Warning: Could not convert %f in less than %d iterations\n',...
        %	 decimal, maxiter+1 );
    end
    %fprintf('dec2frac> num=%d denom=%d success=%d\n', numerator, denominator, success );
end

