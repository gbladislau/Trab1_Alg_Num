function c = signChar(x)
    if x < 0.0
        c='-';
    else
        c='+';
    end
end


% Anonymous implementation

% iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
% signChar = @(x) iif( x < 0.0, @() '-', x >= 0.0, @() '+');

