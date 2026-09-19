function progress_bar(label, k, total)
    % PROGRESS_BAR  In-place terminal progress bar.
    %
    %   progress_bar(label, k, total) draws a fixed-width bar showing the
    %   fraction k/total, overwriting its own line via a carriage return.
    %   Call once per iteration inside a loop; the caller is responsible
    %   for printing a trailing newline (fprintf('\n')) once the loop
    %   finishes (including on early break) so the cursor moves on.
    %
    %   Example:
    %       for i = 1:K
    %           ...
    %           progress_bar('OMP', i, K);
    %       end
    %       fprintf('\n');

    width = 30;
    frac  = max(0, min(1, k / total));
    nfill = round(frac * width);
    bar   = [repmat('#', 1, nfill), repmat('.', 1, width - nfill)];
    fprintf('\r%-6s [%s] %3d%%', label, bar, round(frac * 100));
end
