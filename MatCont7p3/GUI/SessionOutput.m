classdef SessionOutput
    
    properties
        enabled = 1;
        
        mainhandle = figure('Visible','off','UserData', 0); %HACK
        
        pauseoption;
        SPECIAL = 1;
        ALWAYS = 2;
        NEVER = 3;
        
        currentstate = 0;
        READY = 0;
        COMPUTING = 1;
        PAUSED = 2;
        ENDED = 3;
        
        doStop;
        doPause;
        
        plotpointsperiod;
        
        
        
        statushandle = [];
        outputhandle = [];
        
        
        outputlist = {};
        
        model;
        
        %msgs = {};
    end
    
    
    methods
        
        function obj = SessionOutput(session)
            
            if ~isempty(strfind(session.settings.option_pause, 'Special'))
                obj.pauseoption = obj.SPECIAL;
            elseif ~isempty(strfind(session.settings.option_pause, 'Each'))
                obj.pauseoption = obj.ALWAYS;
            else
                obj.pauseoption = obj.NEVER;
            end
                
            obj.doStop = 0;
            obj.doPause = 0;
            
            obj.plotpointsperiod = session.settings.option_output;
            
            obj.model = SessionOutputModel();
        end
        
        
        function output(obj,xout,s,hout,fout,ind)
            global sOutput
            if (sOutput.enabled)
                try
                    for i = 1:length(obj.outputlist)
                        obj.outputlist{i}.output({xout, hout, fout}, s, ind)
                    end
                catch error
                    fprintf(getReport(error));
                    fprintf(2, '\n\n*** Error detected, stopping interactive output...\n\n');
                    sOutput.outputlist = {};
                    sOutput.enabled = 0;
                end
            end
        end
        
        function outputODE(obj, tpart, ypart)
            global sOutput
            if (sOutput.enabled)
                try
                    for i = 1:length(obj.outputlist)
                        obj.outputlist{i}.output({tpart, ypart}, [], 1:length(tpart))
                    end
                catch error
                    fprintf(getReport(error));
                    fprintf(2, '\n\n*** Error detected, stopping interactive output...\n\n');
                    sOutput.outputlist = {};
                    sOutput.enabled = 0;
                end
            end              
        end
        function outputODEEvents(obj, tE, yE, iE)
            global sOutput
            if (sOutput.enabled)
                try
                    for i = 1:length(obj.outputlist)
                        for k = 1:length(tE)
                            obj.outputlist{i}.outputPoint(tE(k), yE(k, :), ['E', num2str(iE(k))]);
                        end
                    end
                catch error
                    fprintf(getReport(error));
                    fprintf(2, '\n\n*** Error detected, stopping interactive output...\n\n');
                    sOutput.outputlist = {};
                    sOutput.enabled = 0;
                end
            end         
            
        end
        
        function setDuration(obj, string)
            global sOutput
            disp(['duration: ' string]);
        end
        
        
        function setStatus(obj, string)
            global sOutput
            set(obj.statushandle , 'String' ,  string);
        end
        
        function endRun(obj)
            global sOutput
            sOutput.setStatus('Finished');
            sOutput.setCurrentState(sOutput.ENDED);

        end
        
        function doStop = checkPauseResumeStop(obj, s_index , s_msg , it)
            global sOutput
            doStop = sOutput.doStop;
            if (doStop)
                sOutput.doStop = 0;
                return;
            end %doStop = 0;
            
            if (sOutput.doPause)
                sOutput.waitForResume();
            end
            
            if (sOutput.pauseoption == sOutput.SPECIAL)
                if ((it > 1) && ( s_index == (it) )) %first point excluded(it >1)
                    sOutput.pause();
                    sOutput.waitForResume(s_msg);
                end
                
            elseif (sOutput.pauseoption == sOutput.ALWAYS)
                sOutput.pause();
                if ((it > 2) && ( s_index == (it) ))
                    sOutput.waitForResume(s_msg);
                else
                    sOutput.waitForResume();
                end
            end
            
        end
        
        function performStop(obj)
            global sOutput
            sOutput.doStop = 1;
            if (sOutput.currentstate == sOutput.PAUSED)
                sOutput.unpause();
            end
        end
        
        
        function waitForResume(obj,msg)
            global sOutput
            sOutput.setCurrentState(sOutput.PAUSED);
            
            if (nargin == 2)
                
                statusmsg = ['Paused, ' strtrim(msg)];
            else
                statusmsg = 'Paused';
            end
            sOutput.setStatus(statusmsg);
            waitfor(sOutput.mainhandle,'UserData',0);
            sOutput.setCurrentState(sOutput.COMPUTING);
            sOutput.setStatus('Computing ...');
            sOutput.doPause = 0;
        end
        
        
        function pause(obj)
            global sOutput
            set(sOutput.mainhandle, 'UserData' , 1);
            sOutput.doPause = 1;
        end
        
        
        function unpause(obj)
            global sOutput
            sOutput.doPause = 0;
            set(sOutput.mainhandle, 'UserData' , 0);
        end
        
        function n = getPlotPointsInterval(obj)
            n = obj.plotpointsperiod;
        end
        
        
        function setPauseAlways(obj)
            global sOutput
            sOutput.pauseoption = sOutput.ALWAYS;
        end
        function setPauseSpecial(obj)
            global sOutput
            sOutput.pauseoption = sOutput.SPECIAL;
        end
        
        function setPauseNever(obj)
            global sOutput
            sOutput.pauseoption = sOutput.NEVER;
        end
        
        function  b = isPauseNever(obj)
            global sOutput
            b = sOutput.pauseoption == sOutput.NEVER;
        end
        function  b = isPauseAlways(obj)
            global sOutput
            b = sOutput.pauseoption == sOutput.ALWAYS;
        end
        function  b = isPauseSpecial(obj)
            global sOutput
            b = sOutput.pauseoption == sOutput.SPECIAL;
        end
        
        function setCurrentState(obj,state)
            global sOutput
            sOutput.currentstate = state;
            sOutput.model.notify('stateChanged');
        end
        
        function b = hasEnded(obj)
            global sOutput
            b = (sOutput.currentstate == obj.ENDED);
        end
        
        function setStatusHandle(obj, statushandle)
            global sOutput
            sOutput.statushandle = statushandle;
        end
        function setOutputHandle(obj, outputhandle)
            global sOutput
            sOutput.outputhandle = outputhandle;
        end
        function delete(obj)
           global sOutput;
           %delete(obj.mainhandle);
           sOutput = [];
        end
    end
    
    methods(Static)
        function s = vector2str(varargin)
            s = ['( ' sprintf('%f ',varargin{1}) ')'];
        end
    end
end


