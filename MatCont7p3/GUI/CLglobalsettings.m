function settings = CLglobalsettings(session)
    settings = CLSettings();
    settings.installSetting('windowpositions', CLSettingWindowPositions());
    settings.installSetting('lineoptions', GUILineOptions(session.branchmanager.getAllCurveLabels()));
    settings.versionnumber = session.VERSION;
    
end