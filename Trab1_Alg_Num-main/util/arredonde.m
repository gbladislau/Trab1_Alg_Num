function arredondado = arredonde( a, casas_decimais )

is_octave = (exist('OCTAVE_VERSION','builtin')>1); % Octave ou Matlab

if is_octave
    a = a.*(10^(casas_decimais));
    a = floor(a);
    a = a.*(10^(-casas_decimais));
else
    a = round(a, casas_decimais);
end
    arredondado = a;
end
