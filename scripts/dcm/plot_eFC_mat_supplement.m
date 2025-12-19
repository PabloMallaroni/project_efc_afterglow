function plot_eFC_mat_supplement(eFC, Pp, regionNames, cmap, clim)
% Plot eFC with low-Pp cells set to 0 (neutral). Show numbers for Pp>=0.90 (★ for >=0.99).

    % --- thresholds
    p_show = 0.5;  p_star = 0.99;  tol = 1e-8;

    % --- coerce & normalise
    eFC = double(eFC);
    Pp  = double(Pp);
    if max(Pp(~isnan(Pp)),[],'all') > 1+1e-6, Pp = Pp/100; end

    % --- colour limits
    if nargin < 5 || isempty(clim)
        v = eFC(isfinite(eFC));
        L = max(abs(v)); if isempty(L) || L==0, L = 1e-3; end
        clim = [-L L];
    elseif isscalar(clim)
        clim = [-abs(clim) abs(clim)];
    else
        clim = clim(:).';
    end

    % --- build mask FIRST, zero low-Pp cells
    mask = isfinite(Pp) & (Pp + tol >= p_show);
    eFC_plot = eFC;
    eFC_plot(~mask) = 0;

    % --- axes & image (use eFC_plot!)
    ax = gca; set(ax,'Color','w');
    imagesc(ax, eFC_plot, clim); axis(ax,'square'); hold(ax,'on');

    % --- colormap
    if (ischar(cmap) || (isstring(cmap) && isscalar(cmap)))
        cname = char(string(cmap));
        if strcmpi(cname,'coolwarm') && exist('coolwarm','file')
            colormap(ax, coolwarm(256));
        else
            colormap(ax, cname);
        end
    elseif isnumeric(cmap) && size(cmap,2)==3
        colormap(ax, cmap);
    else
        colormap(ax, parula(256));
    end

    % --- lock CLim & colourbar ticks (endpoints + 0) in black
    caxis(ax, clim);
    c = colorbar(ax); c.Limits = clim;
    if clim(1) <= 0 && 0 <= clim(2)
        c.Ticks = [clim(1) 0 clim(2)];
    else
        c.Ticks = clim;
    end
    c.Color = 'k'; c.LineWidth = 1.2;

    % --- axis tick labels (always black)
    N = size(eFC,1);
    if nargin < 3 || isempty(regionNames) || numel(regionNames)~=N
        regionNames = arrayfun(@(k) sprintf('R%d',k), 1:N, 'uni', 0);
    else
        if isnumeric(regionNames), regionNames = cellstr(string(regionNames)); end
        if isstring(regionNames),  regionNames = cellstr(regionNames); end
        if ischar(regionNames),    regionNames = cellstr(regionNames); end
    end
    set(ax,'XTick',1:N,'YTick',1:N, ...
           'XTickLabel',regionNames,'YTickLabel',regionNames, ...
           'TickDir','out','LineWidth',1.2,'FontSize',14,'FontWeight','bold', ...
           'Layer','top','XColor','k','YColor','k');

    % --- overlay numbers & stars where mask==true
    fmt = '%.2f';
    for i = 1:N
        for j = 1:N
            if mask(i,j) && isfinite(eFC(i,j))
                text(ax, j, i, sprintf(fmt, eFC(i,j)), ...
                     'HorizontalAlignment','center','VerticalAlignment','middle', ...
                     'FontSize',14,'FontWeight','bold','Color','k','Clipping','on');
                if Pp(i,j) + tol >= p_star
                    text(ax, j+0.32, i-0.32, '*', ...
                         'HorizontalAlignment','center','VerticalAlignment','middle', ...
                         'FontSize',14,'FontWeight','bold','Color','k','Clipping','on');
                end
            end
        end
    end
    hold(ax,'off');
end
