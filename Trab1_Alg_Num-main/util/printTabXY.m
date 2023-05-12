%
% Imprimir tabela X-Y horizontal
% Input:
%	Variáveis X, Y, os seus nomes, string de formatação e título da tabela
% 	Exemplo: printTabXY( X, 'X', Y, 'Y', '%5.2f', 'Tabela X-Y' )
%
function printTabXY( X, Xname, Y, Yname, formstr, titulo )
	fprintf('%s\n', titulo );
	maxlen = max(length(Xname), length(Yname));
	namestr = sprintf('%%%ds%%2s', maxlen+1 );
	fprintf( namestr, Xname, ' |' );
        
        maxlenXY = 0;
        for i=1:length(X)
            maxlenXY = max(maxlenXY, length(sprintf(formstr, X(i))));
            maxlenXY = max(maxlenXY, length(sprintf(formstr, Y(i))));
        end
        
        for i=1:length(X)
             buf = sprintf(formstr, X(i));
             fmt = sprintf('%%%ds', maxlenXY+2 );
             fprintf( fmt, buf );
        end
        fprintf('\n');
        fprintf( '%s', repmat(['-'], 1, maxlen+3) );
        sep = repmat(['-'], 1, maxlenXY+2);
        for i=1:length(X)
            fprintf( '%s', sep );
        end
        fprintf('\n');
        fprintf( namestr, Yname, ' |' );
        for i=1:length(Y)
             buf = sprintf(formstr, Y(i));
             fmt = sprintf('%%%ds', maxlenXY+2 );
             fprintf( fmt, buf );
        end
        fprintf('\n\n');
end


