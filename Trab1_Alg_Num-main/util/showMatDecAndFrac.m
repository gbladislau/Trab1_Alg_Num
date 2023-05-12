%%
%% Mostrar vetor/matriz em forma de fra��es
%%
%%
%% Input: Vetor/Matriz de n�meros decimais, Coment�rio a ser impresso
%% Output: Impress�o

%	maxcase        		n�mero m�ximo campos ap�s a v�rgula (ponto)
%	totalspace        	n�mero total do campos
%   
%%
function showMatDecAndFrac( A, comment )
    maxcase = 3;
    totalspace = 10;
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
                            numPlacesNum = 1 + floor(log10(numerator));		% detectar n�mero de casas decimais
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
                        if success
                                fprintf(elemformstr, elemStr );
                        else
                                fprintf('************  ');
                        end
            end
            % fprintf('\n++++ maxdenoms=%d denoms=', maxdenoms), denoms'
            if maxdenoms > 1 || ~success
                fprintf('\t\t==\t');
                for j = 1:col
                            s = trimdec( A(i,j), maxcase, totalspace );
                            fprintf('%s', s );
                end
            end
            fprintf('\n');
    end
    fprintf('\n');
end
