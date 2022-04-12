function h = GUISettingLabel(parent, settings, settingname, varargin)

    setting = settings.getSetting(settingname);
    h = uicontrol(parent, 'style', 'text', 'String', setting.displayname,'HorizontalAlignment', 'left','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, varargin{:});

    
    helpstr = CLSettingsHelp.getHelp(settingname);
    if ~isempty(helpstr)
        set(h, 'TooltipString', helpstr);
    end

end