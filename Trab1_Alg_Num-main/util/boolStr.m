function bs = boolStr(b)
    if b
        bs='VERDADEIRO';
    else
        bs='FALSO';
    end
end


% Anonymous implementation

% iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
% boolstr = @(bool) iif( bool==true, @() 'TRUE', bool==false, @() 'FALSE');

