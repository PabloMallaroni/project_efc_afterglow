function plot_eFC_delta(eFC_A, Pp_A, eFC_B, Pp_B, regionNames, cmap, clim, fmt, maskMode, p_show)
% Plot ΔeFC = eFC_B - eFC_A using the *raw* matrices (no pre-masking).
% Optional display mask: 'none' (default), 'conj' (both pass), 'union' (either passes).
% Numbers shown in every cell unless you choose to hide non-masked cells.
%
% eFC_A/B : NxN raw effective connectivity (no zeros injected by thresholding)
% Pp_A/B  : NxN posterior probabilities (0..1 or 0..100); can be [] if not used
% regionNames, cmap, clim, fmt : as in your original helper
% maskMode : 'none' | 'conj' | 'union'  (default 'none')
% p_show   : probability threshold (default 0.50)

    if nargin < 9 || isempty(maskMode), maskMode = 'none'; end
    if nargin < 10 || isempty(p_show),   p_show   = 0.50;  end
    if nargin < 8 || isempty(fmt),       fmt      = '%.2f'; end

    eFC_A = double(eFC_A);  eFC_B = double(eFC_B);
    if ~isequal(size(eFC_A), size(eFC_B)), error('A and B must be same size'); end
    N = size(eFC_A,1);

    % --- Δ from RAW matrices (this preserves sign properly)
    D = eFC_B - eFC_A;

    % --- normalise Pp if provided
    havePp = ~(isempty(Pp_A) || isempty(Pp_B));
    if havePp
        Pp_A = double(Pp_A); Pp_B = double(Pp_B);
        if max(Pp_A(:),[],'omitnan') > 1+1e-6, Pp_A = Pp_A/100; end
        if max(Pp_B(:),[],'omitnan') > 1+1e-6, Pp_B = Pp_B/100; end
        switch lower(maskMode)
            case 'conj' % both pass
                mask = (Pp_A >= p_show) & (Pp_B >= p_show);
            case 'union' % either passes
                mask = (Pp_A >= p_show) | (Pp_B >= p_show);
            otherwise
                mask = true(size(D)); % show all
        end
    else
        mask = true(size(D));
    end

    % --- colour limits (symmetric around 0)
    if nargin < 7 || isempty(clim)
        v = D(isfinite(D));
        L = max(abs(v)); if isempty(L) || L==0, L = 1e-3; end
        clim = [-L L];
    elseif isscalar(clim)
        clim = [-abs(clim) abs(clim)];
    else
        clim = clim(:).';
    end

    % --- axes & image
    ax = gca; set(ax,'Color','w');
    imagesc(ax, D, clim); axis(ax,'square'); hold(ax,'on');

    % --- colormap
    if nargin >= 6 && ~isempty(cmap)
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
    else
        colormap(ax, parula(256));
    end

    % --- colourbar
    caxis(ax, clim);
    c = colorbar(ax); c.Limits = clim;
    if clim(1) <= 0 && 0 <= clim(2), c.Ticks = [clim(1) 0 clim(2)]; else, c.Ticks = clim; end
    c.Color = 'k'; c.LineWidth = 1.2; ylabel(c, '\Delta eFC (B - A)');

    % --- labels
    if nargin < 5 || isempty(regionNames) || numel(regionNames) ~= N
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

    % --- numbers in every cell; dim those outside mask (grey text)
    for i = 1:N
        for j = 1:N
            if isfinite(D(i,j))
                col = 'k'; if ~mask(i,j), col = [0.4 0.4 0.4]; end
                text(ax, j, i, sprintf(fmt, D(i,j)), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',14,'FontWeight','bold','Color',col,'Clipping','on');
            end
        end
    end

    % Optional: hatch/outline non-masked cells (commented)
    % if any(~mask(:))
    %     [r,c] = find(~mask);
    %     plot(ax, c+0.5*[ -1 -1 1 1 -1 ], r+0.5*[ -1 1 1 -1 -1 ], 'Color',[.6 .6 .6], 'LineWidth',1);
    % end

    hold(ax,'off');
    title(ax, '\Delta eFC (Condition B - A)', 'FontWeight','bold');
end
