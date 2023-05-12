%%
%% Mostrar Matriz aumentada na Eliminação de Gauss em forma de frações e decimal
%%
%%
%% Input: Matriz e vetor, Comentário a ser impresso
%% Output: Impressão

%	maxcase        		número máximo campos após a vírgula (ponto)
%	totalspace        	número total do campos
%   
%%
function showMatAug( A, b, comment )
    maxcase = 3;
    totalspace = 10;
    [row, col] = size(A);
    [rowb, colb] = size(b);
    if row ~= rowb || colb ~= 1
        error('Dimensao de A e b nao ok: A(%d,%d) b(%d,%d). Saindo ...', row, col, rowb, colb);
    end
    A = [A b];  % Matriz aumentada
    [row, col] = size(A);

    nums = zeros(row,col); denoms = zeros(row,col); % Guardar os numeradores e denominadores em matrizes
    maxPlacesNum = 0; maxPlacesDenom = 0;
    for i = 1:row
            for j = 1:col
                    decimal = A(i,j);
                    [numerator, denominator, success] = dec2frac( decimal );
                    nums(i,j) = numerator;
                    denoms(i,j) = denominator;
                    if numerator > 0
                            numPlacesNum = 1 + floor(log10(numerator));		% detectar número de casas decimais
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

    elemformstr = sprintf( '%%%ds', 4+maxPlacesNum+maxPlacesDenom  ); % formato para um elemento
    fprintf('%s\n\n', comment );

    for i = 1:row
            for j = 1:col
                        if( denoms(i,j) ~= 1 && nums(i,j) ~=0 )
                                elemStr = sprintf('%d/%d', nums(i,j), denoms(i,j) );
                        else
                                elemStr = sprintf('   %d', nums(i,j) );
                        end
                        %fprintf('elemformstr=>>>%s<<< elemStr=>>>%s<<<\n', elemformstr, elemStr );
                        if j == col
                            fprintf('%s', ' | ');
                        end
                        if success
                            fprintf(elemformstr, elemStr );
                        else
                            fprintf('************  ');
                        end
            end
            %fprintf('\n++++ denoms='), denoms'
            if maxdenoms > 1 
                fprintf('\t\t==\t');
                for j = 1:col
                            s = trimdec( A(i,j), maxcase, totalspace );
                            if j == col
                                fprintf('%s', ' | ');
                            end
                            fprintf('%s', s );
                end
            end
            fprintf('\n');
    end
    fprintf('\n');
end
