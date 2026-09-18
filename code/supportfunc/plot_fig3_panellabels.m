function plot_fig3_panellabels(Fig)
% A
Ax = axes(Fig, 'Position', [0, 0.95, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, 'A', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% B
Ax = axes(Fig, 'Position', [0.48, 0.95, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, 'B', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% C
Ax = axes(Fig, 'Position', [0.735, 0.95, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, 'C', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% D
Ax = axes(Fig, 'Position', [0, 0.52, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, 'D', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% E
Ax = axes(Fig, 'Position', [0.48, 0.52, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, 'E', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];
% F
Ax = axes(Fig, 'Position', [0.735, 0.52, 0.01, 0.01]);
Ax.Visible = 'off';
tx = text(Ax, 0, 0, 'F', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment','left', 'VerticalAlignment','top');
tx.Units = 'normalized';
tx.Position = [0 1 0];

end