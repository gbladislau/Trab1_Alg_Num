%%
%% Geração de string de formatação
%% Tentativa de converter float em inteiro, mas sem gerar fração
%%
function fstr = f2i(f)
        [n, d, success] = dec2frac( f );
        fstr = '***';
        if success
            if d == 1
                fstr = sprintf('%s', '%d');
            else
                fstr = sprintf('%s', '%.3f');
            end
        end
end
