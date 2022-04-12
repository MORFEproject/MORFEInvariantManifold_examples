function stop = GUIODEOutput(tpart, xpart, type, varargin)


global sOutput
global Matcont_ODEDATAx;
global Matcont_ODEDATAt;
global Matcont_interval;
global Matcont_eventdata;
global Matcont_options;
global Matcont_eventf;
stop = 0;


%if isempty(sOutput) || ~sOutput.enabled; return; end


if isempty(type)
    Matcont_ODEDATAx = [Matcont_ODEDATAx, xpart];
    Matcont_ODEDATAt = [Matcont_ODEDATAt, tpart];
    
    if length(Matcont_ODEDATAt) >= Matcont_interval
        sOutput.outputODE(Matcont_ODEDATAt, Matcont_ODEDATAx');
        processEvents(Matcont_ODEDATAt, Matcont_ODEDATAx, varargin);
        drawnow;
        Matcont_ODEDATAx = [];
        Matcont_ODEDATAt = [];
    end

    stop = sOutput.checkPauseResumeStop([] , [] , 0);
    
    
elseif strcmp(type, 'done')
    if ~isempty(Matcont_ODEDATAt)
        sOutput.outputODE(Matcont_ODEDATAt, Matcont_ODEDATAx');
        processEvents(Matcont_ODEDATAt, Matcont_ODEDATAx, varargin);
        
    end
    clear Matcont_ODEDATAx Matcont_ODEDATAt Matcont_ODEDATAt Matcont_options Matcont_eventf Matcont_interval;
    sOutput.endRun();
elseif strcmp(type, 'init')
    Matcont_ODEDATAx = [];
    Matcont_ODEDATAt = [];
    Matcont_interval = sOutput.getPlotPointsInterval();
    Matcont_eventdata = [];
    
    if isempty(Matcont_options.Events)
        Matcont_eventf = [];
    else
        Matcont_eventf = Matcont_options.Events;
        Matcont_eventdata = Matcont_eventf(tpart, xpart, varargin{:});
        Matcont_options.OutputFcn = [];
    end
    
    
    sOutput.setStatus('Computing ...');
end


end


function processEvents(t, x, arglist)
    global sOutput Matcont_eventf Matcont_eventdata Matcont_options Matcont_call
    if isempty(Matcont_eventf); return; end
    signchange = 0;
    for k = 1:length(t)
        out = Matcont_eventf(t(k), x(:, k), arglist{:});
        signchange =  signchange || any(Matcont_eventdata .* out <= 0);
        Matcont_eventdata = out;
    end
    
    if signchange
       [~, ~, tE, yE, iE] = Matcont_call([t(1) t(end)], x(:, 1), Matcont_options);
       if ~isempty(tE)
          sOutput.outputODEEvents(tE, yE, iE);
       end
        
    end

end
