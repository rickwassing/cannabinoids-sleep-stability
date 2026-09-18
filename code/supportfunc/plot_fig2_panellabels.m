function plot_fig2_panellabels(Fig, roi)
abc = {'A', 'B', 'C', 'D', 'E', 'F'};
i = 0;

% A
if strcmpi(roi, 'fz')
    i = i+1;
    Ax = axes(Fig, 'Position', [0, 1, 0.01, 0.01]);
    Ax.Visible = 'off';
    tx = text(Ax, 0, 0, abc{i}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
    tx.Units = 'normalized';
    tx.Position = [0 1 0];
end
% B
if strcmpi(roi, 'fz')
    i = i+1;
    Ax = axes(Fig, 'Position', [0.495, 1, 0.01, 0.01]);
    Ax.Visible = 'off';
    tx = text(Ax, 0, 0, abc{i}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
    tx.Units = 'normalized';
    tx.Position = [0 1 0];
end
% C
i = i+1;
Ax = axes(Fig, 'Position', [0, 0.6, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, abc{i}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% D
i = i+1;
Ax = axes(Fig, 'Position', [0.11, 0.6, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, abc{i}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% E
i = i+1;
Ax = axes(Fig, 'Position', [0.5, 0.6, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, abc{i}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% F
i = i+1;
Ax = axes(Fig, 'Position', [0.89, 0.6, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, abc{i}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
end