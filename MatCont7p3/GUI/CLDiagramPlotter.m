function CLDiagramPlotter(diagram, dim, solutionname)
names = diagram.getSolutionNames();
if isempty(names); return; end

if nargin < 3 || ~any(strcmp(solutionname, names))
    solutionname = names{end};
end

%fprintf(2, 'default curve for diagram: %s\n', solutionname);

solution = CompSolution.load(fullfile(diagram.getPath(), solutionname));

oi = solution.compbranch.getOutputInterpreter();
[~, ~, plotselections] = oi.interpret(solution, solution.settings, solution.compbranch);

wm = GUIWindowManager([]);
previousplots = solution.settings.previousplots;
plotconf = previousplots.recover(dim, wm);
if isempty(plotconf)
    plotconf = GUIPlotConf(dim, wm);
end


plotconf.configurePlotSelection(plotselections);
plotconf.setSolutionHandler(CLSolutionHandler.fromSolution(diagram, solution)); %FIXME, path might not be set correctly.

plotconf.generateFigure();

layoutpanel = plotconf.showLayoutWindow();

layoutpanel.waitForClose();
plotconf.isValid()

plotconf.plotDiagram()

end
