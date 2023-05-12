%%
%% Conversão de um número decimal em número expresso como fração
%% ocupando 'numcasas' casas decimais
%% algoritmo de conversão com precisão e número máximo de iterções
%%
function s = dec2fracstr( x, numcasas )
	[num,denom,success] = dec2frac( x );

	if success
		if denom == 1
			formstr = sprintf('%%d'); 
			fracstr = sprintf(formstr,num);
		elseif num == 0
			fracstr= '0';
		else
			formstr = sprintf('%%d/%%d');
			fracstr = sprintf(formstr,num,denom);
		end
		casasstr = sprintf('%%%ds',numcasas);
%		x
%		num
%		denom
%		formstr
%		fracstr
%		casasstr
		s = sprintf(casasstr, fracstr);
	else		
		for i=1:numcasas s(i) = '*'; end
	end

end		

