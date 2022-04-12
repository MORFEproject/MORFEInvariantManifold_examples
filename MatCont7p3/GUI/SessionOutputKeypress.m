function SessionOutputKeypress(handle, event)
global sOutput
if (isempty(sOutput) || ~sOutput.enabled); return; end


if strcmp(event.Key, 'escape')
    sOutput.performStop();
elseif strcmp(event.Key, 'return')
    sOutput.pause()
elseif strcmp(event.Key, 'space')
    sOutput.unpause();
end

end