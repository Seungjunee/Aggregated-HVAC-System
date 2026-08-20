function files = plot_paper_figures(series, outputDir)
%PLOT_PAPER_FIGURES Generate Figures 5-13. Figure 4 is intentionally excluded.

if ~isfolder(outputDir)
    mkdir(outputDir);
end

files = strings(9, 1);
files(1) = save_figure(figure5(series), outputDir, 'figure_05_baseline_power');
files(2) = save_figure(figure6(series), outputDir, 'figure_06_vess_power');
files(3) = save_figure(figure7(series), outputDir, 'figure_07_soc_regulation');
files(4) = save_figure(figure8(series), outputDir, 'figure_08_vb_power');
files(5) = save_figure(figure9(series), outputDir, 'figure_09_temperature');
files(6) = save_figure(figure10(series), outputDir, 'figure_10_soz');
files(7) = save_figure(figure11(series), outputDir, 'figure_11_soc_setpoint');
files(8) = save_figure(figure12(series), outputDir, 'figure_12_soc_mass');
files(9) = save_figure(figure13(series), outputDir, 'figure_13_soc_vess');
end

function fig = new_figure(width, height)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 width height]);
end

function style_axes(ax)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
xlim(ax, [10 18]);
xlabel(ax, 'Time [h]');
end

function fig = figure5(s)
fig = new_figure(760, 360);
ax = axes(fig);
plot(ax, s.time, s.vb.b1.baselineActual, '-', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.10]); hold(ax, 'on');
plot(ax, s.time, s.vb.b1.baseline, ':', 'LineWidth', 2.5, 'Color', [0.85 0.33 0.10]);
plot(ax, s.time, s.vb.b11.baselineActual, '-', 'LineWidth', 2.5, 'Color', [0.00 0.45 0.75]);
plot(ax, s.time, s.vb.b11.baseline, ':', 'LineWidth', 2.5, 'Color', [0.00 0.45 0.75]);
style_axes(ax); ylabel(ax, 'Baseline Power [kW]');
legend(ax, {'P_1^{base} (Eplus)', 'P_1^{base} (predict)', ...
    'P_{11}^{base} (Eplus)', 'P_{11}^{base} (predict)'}, ...
    'Location', 'south', 'NumColumns', 2);
end

function fig = figure6(s)
fig = new_figure(780, 450);
ax = axes(fig);
plot(ax, s.time, s.powerSignal, '-', 'LineWidth', 3, 'Color', [0.58 0.26 0.96]); hold(ax, 'on');
plot(ax, s.time, s.baseline, 'k-', 'LineWidth', 2.5);
plot(ax, s.time, s.referencePower, 'r--', 'LineWidth', 2);
plot(ax, s.time, s.actualPowerSetpoint, '-', 'LineWidth', 2.5, 'Color', [0.00 0.45 0.74]);
plot(ax, s.time, s.actualPowerMass, '-', 'LineWidth', 2.5, 'Color', [0.93 0.69 0.13]);
plot(ax, s.time, s.baseline - s.dischargeLimit, '-', 'LineWidth', 0.7, 'Color', [0.65 0.65 0.65]);
plot(ax, s.time, s.baseline + s.chargeLimit, '-', 'LineWidth', 0.7, 'Color', [0.65 0.65 0.65]);
style_axes(ax); ylabel(ax, 'Power [kW]'); ylim(ax, [-100 200]);
legend(ax, {'P_{VESS}^{ch/dch}', 'P_{VESS}^{baseline}', 'P_{VESS}^{ref}', ...
    'P_{VESS}^{setpoint}', 'P_{VESS}^{mass}', 'Limits'}, ...
    'Location', 'southoutside', 'NumColumns', 3);
end

function fig = figure7(s)
fig = new_figure(700, 340);
ax = axes(fig);
plot(ax, s.lambdaTime, s.lambda, 'k--', 'LineWidth', 2.2);
style_axes(ax); ylabel(ax, 'SOC_{VESS} [-]'); ylim(ax, [-2 2]);
legend(ax, {'\lambda'}, 'Location', 'northeast');
end

function fig = figure8(s)
fig = new_figure(1050, 380);
for panel = 1:2
    ax = subplot(1, 2, panel, 'Parent', fig);
    if panel == 1, v = s.vb.b1; else, v = s.vb.b11; end
    plot(ax, s.time, v.baseline, 'k-', 'LineWidth', 2.5); hold(ax, 'on');
    plot(ax, s.time, v.reference, 'r--', 'LineWidth', 2);
    plot(ax, s.time, v.actualSetpoint, '-', 'LineWidth', 2.5, 'Color', [0.00 0.45 0.74]);
    plot(ax, s.time, v.actualMass, '-', 'LineWidth', 2.5, 'Color', [0.93 0.69 0.13]);
    style_axes(ax); ylabel(ax, 'Power [kW]'); ylim(ax, [5 9]);
    title(ax, sprintf('(%c)', char('a' + panel - 1)));
    legend(ax, {'P^{baseline}', 'P^{ref}', 'P^{setpoint}', 'P^{mass}'}, ...
        'Location', 'best', 'NumColumns', 2);
end
end

function fig = figure9(s)
fig = new_figure(1050, 380);
plot_temperature_panel(subplot(1, 2, 1, 'Parent', fig), s.time, ...
    s.temperature.setpoint.b1, s.temperature.setpoint.b11, '(a) Setpoint regulation');
plot_temperature_panel(subplot(1, 2, 2, 'Parent', fig), s.time, ...
    s.temperature.mass.b1, s.temperature.mass.b11, '(b) Mass flow rate');
end

function plot_temperature_panel(ax, time, threeZone, fiveZone, panelTitle)
colors = lines(8);
values = [threeZone, fiveZone];
for k = 1:size(values, 2)
    plot(ax, time, values(:, k), 'LineWidth', 2.2, 'Color', colors(k, :)); hold(ax, 'on');
end
style_axes(ax); ylabel(ax, 'Temperature [degC]'); ylim(ax, [20 30]); title(ax, panelTitle);
end

function fig = figure10(s)
fig = new_figure(1050, 380);
plot_soz_panel(subplot(1, 2, 1, 'Parent', fig), s.time, ...
    s.soz.setpoint.b1, s.soz.setpoint.b11, '(a) Setpoint regulation');
plot_soz_panel(subplot(1, 2, 2, 'Parent', fig), s.time, ...
    s.soz.mass.b1, s.soz.mass.b11, '(b) Mass flow rate');
end

function plot_soz_panel(ax, time, threeZone, fiveZone, panelTitle)
colors = lines(8);
values = [threeZone, fiveZone];
for k = 1:size(values, 2)
    plot(ax, time, values(:, k), 'LineWidth', 2, 'Color', colors(k, :)); hold(ax, 'on');
end
plot(ax, time, ones(size(time)), 'r--', 'LineWidth', 1.5);
plot(ax, time, -ones(size(time)), 'r--', 'LineWidth', 1.5);
style_axes(ax); ylabel(ax, 'SOZ [-]'); ylim(ax, [-2 2]); title(ax, panelTitle);
end

function fig = figure11(s)
fig = soc_profiles(s.time, s.socByBuildingSetpoint);
end

function fig = figure12(s)
fig = soc_profiles(s.time, s.socByBuildingMass);
end

function fig = soc_profiles(time, profiles)
fig = new_figure(760, 360);
ax = axes(fig);
plot(ax, time, profiles', 'LineWidth', 1.8); hold(ax, 'on');
hLimit = plot(ax, time, ones(size(time)), 'r--', 'LineWidth', 1.5);
plot(ax, time, -ones(size(time)), 'r--', 'LineWidth', 1.5);
style_axes(ax); ylabel(ax, 'SOC [-]'); ylim(ax, [-2 2]);
legend(ax, hLimit, {'Limits'}, 'Location', 'northeast');
end

function fig = figure13(s)
fig = new_figure(760, 360);
ax = axes(fig);
plot(ax, s.lambdaTime, s.lambda, 'k--', 'LineWidth', 2); hold(ax, 'on');
plot(ax, s.time, s.socVessSetpoint, '-', 'LineWidth', 2.3, 'Color', [0.00 0.45 0.74]);
plot(ax, s.time, s.socVessMass, '-', 'LineWidth', 2.3, 'Color', [0.93 0.69 0.13]);
plot(ax, s.time, ones(size(s.time)), 'r--', 'LineWidth', 1.3);
plot(ax, s.time, -ones(size(s.time)), 'r--', 'LineWidth', 1.3);
style_axes(ax); ylabel(ax, 'SOC_{VESS} [-]'); ylim(ax, [-2 2]);
legend(ax, {'\lambda', 'soc_{VESS}^{setpoint}', 'soc_{VESS}^{mass}', 'Limits'}, ...
    'Location', 'northeast', 'NumColumns', 2);
end

function pngPath = save_figure(fig, outputDir, baseName)
set(fig, 'PaperPositionMode', 'auto');
pngPath = fullfile(outputDir, [baseName '.png']);
figPath = fullfile(outputDir, [baseName '.fig']);
print(fig, pngPath, '-dpng', '-r200');
savefig(fig, figPath);
close(fig);
end
