function printdecandfrac( decimal, newline )
    [num, denom, success] = dec2frac( decimal );

    %fprintf('printdecandfrac> num=%d denom=%d success=%d\n', num, denom, success );

    dec = sprintf('%.4f', decimal );
    if success
        if denom == 1
            frac = sprintf('%d', num );
        else
            %fprintf('%.4f', decimal );
            frac = sprintf('%d/%d', num, denom );
        end;
    %else
    %    frac = sprintf('=*****');
    end
    
    if success
        if denom == 1
            fprintf('%s', frac)
        else
            fprintf('%s=%s', dec, frac)
        end
    else
        fprintf('%s', dec)
    end    
    
    if newline fprintf('\n'); end
end


