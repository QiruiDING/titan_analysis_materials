% ============================================================================
% Titan Submersible Comprehensive Failure Analysis
% Word-ready (5.0 x 3.5 in @ 300 dpi), clean style, with preview windows on
% Figures 1/2/3 have no "box-analysis" panels
% ============================================================================

clear; close all; clc;

% ------------------------ Minimalist global styling -------------------------
set(groot, ...
    'DefaultFigureColor','w', ...
    'DefaultAxesFontName','Helvetica', ...
    'DefaultAxesFontSize',8, ...
    'DefaultTextFontName','Helvetica', ...
    'DefaultLineLineWidth',1.4, ...
    'DefaultAxesLineWidth',0.9, ...
    'DefaultAxesTickDir','out', ...
    'DefaultAxesBox','off');

% Word single-column size
WORD_W_IN = 5.0;    % width in inches
WORD_H_IN = 3.5;    % height in inches
DPI       = 300;

% ----------------------------- Data / material ------------------------------
mesh_sizes          = [1024, 5838, 112271, 357132, 500000];  % 添加第5个数据点
total_deformation   = [0.03748, 0.035229, 0.039663, 0.039664, 0.039665];
ti_stress           = [3.9854e8, 4.1649e8, 4.2447e8, 4.2476e8, 4.2478e8];
cfrp_hoop           = [-1.5792e8, -1.5764e8, -1.5689e8, -1.5677e8, -1.5676e8];
cfrp_axial          = [-9.2523e7, -8.4993e7, -7.9049e7, -7.8166e7, -7.8100e7];
equivalent_stress   = [3.9854e8, 4.1649e8, 4.2447e8, 4.2476e8, 4.2478e8];  % 添加等效应力数据

CFRP.sigma_L = 1560e6;
CFRP.sigma_T = 50e6;
CFRP.tau_LT  = 90e6;

depth = 4000;
pressure_operational = depth * 9.81 * 1000;

% ============================================================================
% FIGURE 1: Mesh convergence (2x3布局, 显示5个数据集)
% ============================================================================
f1 = newfig('Figure 1: Mesh Convergence (Five Data Sets)', WORD_W_IN, WORD_H_IN);
t1 = tiledlayout(f1,2,3,'Padding','compact','TileSpacing','compact');

% (a) deformation convergence with subtle band
nexttile(t1,1); hold on;
unc = total_deformation*0.01;
fill_between(mesh_sizes, total_deformation-unc, total_deformation+unc, [0.85 0.90 1.00], 1.0);
loglog(mesh_sizes, total_deformation, 'o-','Color',hex2rgb('#1f77b4'),'MarkerSize',6);
xlabel('Mesh elements'); ylabel('Total deformation (mm)');
grid on; set(gca,'GridAlpha',0.12);

% (b) relative error decay
nexttile(t1,2);
err_def = [abs((total_deformation(1)-total_deformation(2))/total_deformation(2))*100, ...
           abs((total_deformation(2)-total_deformation(3))/total_deformation(3))*100, ...
           abs((total_deformation(3)-total_deformation(4))/total_deformation(4))*100, ...
           abs((total_deformation(4)-total_deformation(5))/total_deformation(5))*100];
semilogy(1:4, err_def, 'o-','Color',hex2rgb('#d62728'),'MarkerSize',6);
xlim([0.8 4.2]); xlabel('Refinement step'); ylabel('Relative error (%)');
grid on; set(gca,'GridAlpha',0.12);

% (c) stress convergence
nexttile(t1,3); hold on;
loglog(mesh_sizes, ti_stress/1e8, 'o-','Color',hex2rgb('#ff7f0e'),'MarkerSize',6,'DisplayName','Ti-6Al-4V');
loglog(mesh_sizes, abs(cfrp_hoop)/1e8, 's-','Color',hex2rgb('#9467bd'),'MarkerSize',6,'DisplayName','CFRP hoop');
loglog(mesh_sizes, abs(cfrp_axial)/1e8,'^-','Color',hex2rgb('#2ca02c'),'MarkerSize',6,'DisplayName','CFRP axial');
xlabel('Mesh elements'); ylabel('Stress (×10^8 Pa)');
legend('Location','best','Box','off');
grid on; set(gca,'GridAlpha',0.12);

% (d) CFRP hoop stress (单独显示)
nexttile(t1,4);
plot(1:5, abs(cfrp_hoop)/1e8, 'o-','Color',hex2rgb('#9467bd'),'MarkerSize',6,'LineWidth',1.4);
xlabel('Index'); ylabel('Value');
title('CFRP hoop stress');
grid on; set(gca,'GridAlpha',0.12);
xlim([0.8 5.2]);

% (e) CFRP axial stress (单独显示)  
nexttile(t1,5);
plot(1:5, abs(cfrp_axial)/1e7, 'o-','Color',hex2rgb('#2ca02c'),'MarkerSize',6,'LineWidth',1.4);
xlabel('Index'); ylabel('Value');
title('CFRP axial stress');
grid on; set(gca,'GridAlpha',0.12);
xlim([0.8 5.2]);

% (f) Equivalent stress (单独显示)
nexttile(t1,6);
plot(1:5, equivalent_stress/1e8, 'o-','Color',hex2rgb('#ff7f0e'),'MarkerSize',6,'LineWidth',1.4);
xlabel('Index'); ylabel('Value');
title('Equivalent stress');
grid on; set(gca,'GridAlpha',0.12);
xlim([0.8 5.2]);

exportgraphics(f1, 'Figure1_MeshConvergence.png','Resolution',DPI);
drawnow;

% ============================================================================
% FIGURE 2: Stress field visualization (clean 1x3)
% ============================================================================
f2 = newfig('Figure 2: Stress Field', WORD_W_IN, WORD_H_IN);
t2 = tiledlayout(f2,1,3,'Padding','compact','TileSpacing','compact');

% Lamé solution through thickness
r_inner = 1.0; r_outer = 1.2;
theta_grid = linspace(0,2*pi,160);
r_grid     = linspace(r_inner,r_outer,80);
[Theta, R] = meshgrid(theta_grid, r_grid);
A = pressure_operational * r_outer^2 / (r_outer^2 - r_inner^2);
B = -pressure_operational * r_outer^2 * r_inner^2 / (r_outer^2 - r_inner^2);
sigma_r     = A + B ./ (R.^2);
sigma_theta = A - B ./ (R.^2);
VM = sqrt(sigma_r.^2 + sigma_theta.^2 - sigma_r.*sigma_theta);

% (a) hoop stress heatmap
axA = nexttile(t2,1);
H = sigma_theta/1e8;
p = pcolor(axA, Theta, R, H); set(p,'EdgeColor','none');
colormap(axA, parula); cb = colorbar(axA); cb.Label.String = 'Hoop stress (×10^8 Pa)';
xlabel('\theta (rad)'); ylabel('r (m)');

% (b) through-thickness lines
axB = nexttile(t2,2); hold(axB,'on');
r_line = linspace(r_inner,r_outer,200);
sig_r_line = A + B./(r_line.^2);
sig_t_line = A - B./(r_line.^2);
plot(axB, r_line, sig_r_line/1e8, 'o-','MarkerSize',5,'Color',hex2rgb('#1f77b4'),'DisplayName','\sigma_r');
plot(axB, r_line, sig_t_line/1e8, 's-','MarkerSize',5,'Color',hex2rgb('#9467bd'),'DisplayName','\sigma_\theta');
xlabel('r (m)'); ylabel('Stress (×10^8 Pa)');
legend(axB,'Location','best','Box','off'); grid(axB,'on'); set(axB,'GridAlpha',0.12);

% (c) von Mises vs r
axC = nexttile(t2,3); hold(axC,'on');
VM_line = sqrt(sig_r_line.^2 + sig_t_line.^2 - sig_r_line.*sig_t_line);
plot(axC, r_line, VM_line/1e8, '-','Color',hex2rgb('#d62728'),'LineWidth',1.6);
xlabel('r (m)'); ylabel('Von Mises (×10^8 Pa)');
grid(axC,'on'); set(axC,'GridAlpha',0.12);

exportgraphics(f2, 'Figure2_StressField.png','Resolution',DPI);
drawnow;

% ============================================================================
% FIGURE 3: Composite failure criteria (clean 2x2; no box panels)
% ============================================================================
f3 = newfig('Figure 3: Failure Criteria', WORD_W_IN, WORD_H_IN);
t3 = tiledlayout(f3,2,2,'Padding','compact','TileSpacing','compact');

% normalized stresses at finest mesh
sigma_L = abs(cfrp_axial(end)) / CFRP.sigma_L;
sigma_T = abs(cfrp_hoop(end))  / CFRP.sigma_T;
F_LT    = -0.5;
tsai_wu_fail       = sigma_L^2 + sigma_T^2 + F_LT*sigma_L*sigma_T;
hashin_fiber_fail  = sigma_L^2;
hashin_matrix_fail = (abs(cfrp_hoop(end))/CFRP.tau_LT)^2 + (abs(cfrp_axial(end))/CFRP.sigma_T)^2;

% (a) Tsai-Wu envelope
nexttile(t3,1);
sigL = linspace(0,1.5,220); sigT = linspace(0,3.0,220);
[SigL, SigT] = meshgrid(sigL, sigT);
TW = SigL.^2 + SigT.^2 + F_LT.*SigL.*SigT;
contourf(sigL, sigT, TW, 18, 'LineStyle','none'); hold on;
contour(sigL, sigT, TW, [1 1], 'k-','LineWidth',1.2);
plot(sigma_L, sigma_T, 'k*','MarkerSize',9,'LineWidth',1.1);
colormap(parula); cb = colorbar; cb.Label.String = 'Tsai-Wu index';
xlabel('\sigma_L / \sigma_{L,allow}'); ylabel('\sigma_T / \sigma_{T,allow}');

% (b) Hashin matrix envelope
nexttile(t3,2);
HM = (SigT*CFRP.sigma_T/CFRP.tau_LT).^2 + (SigL*CFRP.sigma_L/CFRP.sigma_T).^2;
contourf(sigL, sigT, HM, 18, 'LineStyle','none'); hold on;
contour(sigL, sigT, HM, [1 1], 'k-','LineWidth',1.2);
plot(sigma_L, sigma_T, 'k*','MarkerSize',9,'LineWidth',1.1);
colormap(parula); cb = colorbar; cb.Label.String = 'Hashin matrix index';
xlabel('\sigma_L / \sigma_{L,allow}'); ylabel('\sigma_T / \sigma_{T,allow}');

% (c) Failure indices vs mesh
nexttile(t3,3); hold on;
tsai_wu_hist = zeros(1,numel(mesh_sizes));
hashin_fiber_hist = zeros(1,numel(mesh_sizes));
hashin_matrix_hist = zeros(1,numel(mesh_sizes));
for i = 1:numel(mesh_sizes)
    sL = abs(cfrp_axial(i))/CFRP.sigma_L;
    sT = abs(cfrp_hoop(i))/CFRP.sigma_T;
    tsai_wu_hist(i) = sL^2 + sT^2 + F_LT*sL*sT;
    hashin_fiber_hist(i) = sL^2;
    hashin_matrix_hist(i) = (abs(cfrp_hoop(i))/CFRP.tau_LT)^2 + (abs(cfrp_axial(i))/CFRP.sigma_T)^2;
end
plot(1:5, tsai_wu_hist, 'o-','Color',hex2rgb('#1f77b4'),'MarkerSize',6,'DisplayName','Tsai-Wu');
plot(1:5, hashin_fiber_hist, 's-','Color',hex2rgb('#2ca02c'),'MarkerSize',6,'DisplayName','Hashin fiber');
plot(1:5, hashin_matrix_hist, '^-','Color',hex2rgb('#d62728'),'MarkerSize',6,'DisplayName','Hashin matrix');
yline(1,'k:');
xlim([0.8 5.2]); ylim([0 max([tsai_wu_hist,hashin_fiber_hist,hashin_matrix_hist])*1.15]);
set(gca,'XTick',1:5,'XTickLabel',{'1k','6k','112k','357k','500k'});
xlabel('Mesh refinement'); ylabel('Failure index');
legend('Location','best','Box','off'); grid on; set(gca,'GridAlpha',0.12);

% (d) Polar comparison (use PolarAxes in the same tile position)
tmpAx = nexttile(t3,4); pos = get(tmpAx,'Position'); delete(tmpAx);
axP = polaraxes('Position',pos); hold(axP,'on');
theta_vals = [0, 2*pi/3, 4*pi/3, 2*pi];
r_vals = [tsai_wu_fail, hashin_fiber_fail, hashin_matrix_fail, tsai_wu_fail];
polarplot(axP, theta_vals, r_vals, 'o-','MarkerSize',6,'Color',hex2rgb('#1f77b4'));
polarplot(axP, theta_vals, ones(size(theta_vals)), 'k--','LineWidth',1.0);
thetaticks(axP,[0 120 240]); thetaticklabels(axP, {'Tsai-Wu','Hashin-F','Hashin-M'});
rlim(axP,[0 max(2, max(r_vals)*1.1)]);

exportgraphics(f3, 'Figure3_FailureCriteria.png','Resolution',DPI);
drawnow;

% ============================================================================
% FIGURE 4: Heatmap • Fit • 3D (same compact size; optional advanced views)
% ============================================================================
f4 = newfig('Figure 4: Heatmap • Fit • 3D', WORD_W_IN, WORD_H_IN);
t4 = tiledlayout(f4,1,3,'Padding','compact','TileSpacing','compact');

% (a) failure-index heatmap across meshes
nexttile(t4,1);
Hdata = [tsai_wu_hist; hashin_fiber_hist; hashin_matrix_hist]'; % 5x3
h = heatmap({'Tsai-Wu','Hashin-F','Hashin-M'}, {'1k','6k','112k','357k','500k'}, Hdata);
h.Colormap = parula; h.ColorbarVisible = 'on';
h.Title = ''; h.XLabel = 'Criterion'; h.YLabel = 'Mesh';

% (b) power-law fit (log-log)
nexttile(t4,2); hold on;
X = mesh_sizes(:); Y = total_deformation(:);
pf = polyfit(log10(X), log10(Y), 1);
yhat = 10.^(polyval(pf, log10(X)));
loglog(X, Y, 'o','MarkerSize',6,'Color',hex2rgb('#d62728'));
loglog(X, yhat, '-','LineWidth',1.4,'Color',hex2rgb('#1f77b4'));
xlabel('Mesh elements'); ylabel('Total deformation (mm)');
legend('FEA','Power-law fit','Location','best','Box','off');
grid on; set(gca,'GridAlpha',0.12);

% (c) compact 3D von Mises surface
nexttile(t4,3);
[Xc, Yc] = pol2cart(Theta, R);
s = surf(Xc, Yc, VM/1e8, 'EdgeColor','none'); 
shading interp; axis equal tight; view(38,28);
colormap(parula); cb = colorbar; cb.Label.String = 'Von Mises (×10^8 Pa)';
xlabel('x (m)'); ylabel('y (m)'); zlabel('VM (×10^8 Pa)');

exportgraphics(f4, 'Figure4_Heatmap_Fit_Complex.png','Resolution',DPI);
drawnow;

disp('Exported PNGs (300 dpi): Figure1_MeshConvergence.png, Figure2_StressField.png, Figure3_FailureCriteria.png, Figure4_Heatmap_Fit_Complex.png');

% ============================================================================
% Helper functions
% ============================================================================
function f = newfig(name, w_in, h_in)
    f = figure('Name',name,'Visible','on','Units','inches','PaperPositionMode','auto');
    f.Position = [1 1 w_in h_in];  % exact physical size for Word + preview window sized accordingly
end

function rgb = hex2rgb(hex)
    if hex(1) == '#', hex = hex(2:end); end
    rgb = [hex2dec(hex(1:2))/255, hex2dec(hex(3:4))/255, hex2dec(hex(5:6))/255];
end

function fill_between(x, lower, upper, face_rgb, face_alpha)
    x = x(:)'; lower = lower(:)'; upper = upper(:)';
    X = [x, fliplr(x)];
    Y = [lower, fliplr(upper)];
    p = fill(X, Y, face_rgb, 'EdgeColor','none');
    p.FaceAlpha = face_alpha;
end
