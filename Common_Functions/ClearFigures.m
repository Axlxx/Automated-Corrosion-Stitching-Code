function ClearFigures()
    % Clear all Figures
    figs = findall(groot,'Type','figure');
    for f = figs
        clf(f);
    end
end

