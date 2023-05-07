%%
%% String com somente as partes fracionárias diferente de zero de um número decimal
%%
%% Input: Número decimal, número máximo de casas
%% Output: string
%%
function s = trimdec( x, maxcase, totalspace )
	frac = x - floor(x);
	if frac == 0
		formstr = sprintf('%%%dd',totalspace);
		s = sprintf(formstr, x );
		return;
	end


	formstr = sprintf('%%.%df',maxcase);
	s = sprintf(formstr, x );
	[beforeperiod, fromperiod] = strtok(s,'.');

	fs = fromperiod(2:end);
	c = 1;

	while c < maxcase && c < totalspace-length(beforeperiod)-1 && c < length(fs) && fs(c+1) ~= '0'
		c = c+1;
	end	

	s = [beforeperiod '.' fs(1:c)];		% length(s)

	formstr = sprintf('%%%ds',totalspace);
	s = sprintf(formstr,s);				% length(s)
end

