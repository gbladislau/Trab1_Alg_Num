%%
%% Mostrar vetor/matriz em forma de frações
%%
%%
%% Input: Vetor/Matriz de números decimais, Comentário a ser impresso
%% Índice a partir de qual os multiplicadores serão coloridos
%% Output: Impressão

%	maxcase        		número máximo campos após a vírgula (ponto)
%	totalspace        	número total do campos
%   
%%
function showMatLU( A, comment, colunascoloridas )
    % Tenta ler variável de ambiente 'OCTAVE_SEM_COR'. Se achar imprima sem cor
    if length(getenv('OCTAVE_SEM_COR')) > 0
        showMatDecAndFrac( A, comment );
        return;
    end
    if ispc % Windows ?
        % Para ter cor no terminal, esta chave no registro do Windows
        % tem que existir, e ter valor igual a um
        % https://superuser.com/questions/413073/windows-console-with-ansi-colors-handling
        enableANSIcolor = 'reg add HKEY_CURRENT_USER\\Console /v VirtualTerminalLevel /t REG_DWORD /d 1 /f\n';
        enableANSIcolorCmd = 'reg add HKEY_CURRENT_USER\Console /v VirtualTerminalLevel /t REG_DWORD /d 1 /f';
        info1 = '\nFazendo registro no Windows para permitir cores:\n';
        info2 = 'Na proxima vez apos carregar o Octave, ja deve funcionar!\n';
        try
            if ~winqueryreg('HKEY_CURRENT_USER', 'Console', 'VirtualTerminalLevel')
                printf(info1);
                printf(enableANSIcolor);
                printf(info2);
                dos('reg add HKEY_CURRENT_USER\Console /v VirtualTerminalLevel /t REG_DWORD /d 1 /f');
                showMatDecAndFrac( A, comment );
                return;
            end
        catch
              printf(info1);
              printf(enableANSIcolor);
              printf(info2);
              dos('reg add HKEY_CURRENT_USER\Console /v VirtualTerminalLevel /t REG_DWORD /d 1 /f');
              showMatDecAndFrac( A, comment );
        end
    end
    maxcase = 3;
    [row, col] = size(A);
    nums = zeros(row,col); denoms = zeros(row,col); % Guardar os numeradores e denominadores em matrizes
    maxPlacesNum = 0; maxPlacesDenom = 0;
    for i = 1:row
            for j = 1:col
                    decimal = A(i,j);
                    [numerator, denominator, success] = dec2frac( decimal );
                    nums(i,j) = numerator;
                    denoms(i,j) = denominator;
                    % detectar número de casas decimais
                    if numerator > 0
                            numPlacesNum = 1 + floor(log10(numerator));		
                    else
                            numPlacesNum = 1;
                    end
                    if denominator > 0
                            numPlacesDenom = 1 + floor(log10(denominator));
                    else
                            numPlacesDenom = 1;
                    end
                    if numPlacesNum > maxPlacesNum
                            maxPlacesNum = numPlacesNum;
                    end
                    if numPlacesDenom > maxPlacesDenom
                            maxPlacesDenom = numPlacesDenom;
                    end
            end
    end
    maxdenoms = max(max(denoms)); % se for maior que um, tem que imprimir a forma decimal

    fprintf('%s\n\n', comment );
    
    casas = 6 + maxPlacesNum + maxPlacesDenom;
    totalspace = max(10, casas+2);
    % https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences
    % https://pt.wikipedia.org/wiki/C%C3%B3digo_escape_ANSI
    cor = '\033[102m';
    cor = '\033[43m';
    correset = '\033[0m';
    %fprintf('Casas max: Nominador=%d, Denomindador=%d, casas=%d\n', maxPlacesNum, maxPlacesDenom, casas)
    for i = 1:row
            for j = 1:col
                        pintar = j < i && j <= colunascoloridas;
                        if( denoms(i,j) ~= 1 && nums(i,j) ~=0 )
                                numstr = sprintf('%d/%d', nums(i,j), denoms(i,j) );
                        else
                                numstr = sprintf('%d',  nums(i,j) );
                        end
                        if success
                                blankformstr = sprintf('%%%ds', casas - length(numstr));
                                blanks = sprintf(blankformstr, ' ');
                                %fprintf('Num casas=%d blankformstr=%s\n', casas, blankformstr)
                                if ~pintar
                                    fprintf(blanks); fprintf(numstr);
                                else
                                    numstrcolorido = [cor blanks numstr correset];
                                    fprintf(numstrcolorido);
                                end
                        else
                                if ~pintar
                                    fprintf('  ************');
                                else
                                    numstrcolorido = ['  ' cor '************' correset];
                                    fprintf(numstrcolorido);
                                end
                        end
            end
            if maxdenoms > 1 
                fprintf('\t\t==\t');
                for j = 1:col
                            pintar = j < i && j <= colunascoloridas;
                            s = trimdec( A(i,j), maxcase, totalspace );
                            if ~pintar
                                fprintf('%s', s );
                            else
                                numstrcolorido = [cor s correset];
                                fprintf(numstrcolorido);
                            end
                end
            end
            fprintf('\n');
    end
    fprintf('\n');
end
