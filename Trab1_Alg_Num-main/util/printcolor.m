
% Imprimir x colorido no formato especificado em formstr 
function printcolor(x, formstr)
  %  colstr = "\033[31;1;4m"; % red, bold, underline
    colstr = "\033[31;1m";    % red bold
    valstr = sprintf(formstr,x);
    resetstr = "\033[0m";
    termstr = sprintf("%s%s%s%s", resetstr, colstr, valstr, resetstr );
    fprintf(1,termstr)
end
%  fprintf(1,"\033[31;1mHello Pi=%10.5e\n\033[0m",pi);



% http://ascii-table.com/ansi-escape-sequences.php

% https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences



% The ANSI escape sequences are the Select Graphic Rendition subset.
% All of these have the form

% \033[XXXm

% where XXX is a series of semicolon-separated parameters.

% To say, make text red, bold, and underlined in Octave/Matlab:

% fprintf(1,"\033[31;1;4mHello\033[0m");
