function sstr = signString( x )
	if x < 0.0
		sstr = '-';
	elseif x > 0.0
		sstr = '+';
	else
		sstr = '';
	end
end


