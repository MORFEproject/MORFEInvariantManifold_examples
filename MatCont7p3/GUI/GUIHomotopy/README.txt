Homotopy was first implemented in an older version of the MatcontGUI (version 6.xx and earlier).

Due to the nature of the Homotopy-method, the algorithms are not clearly separated from the graphical aspect. There is no
command-line version of the Homotopy-method. The Homotopy (HT) code is interwoven in the GUI and makes integration into the new GUI a hard task. 

In order to avoid a complete rewrite, the HT-code is embedded into the New-GUI where the HT-code acts like it is acting as if in the old version, a kind of 'emulation'.  This GUI sets up the situation as if in the old gui, runs the HT-code and then translates back into the this GUI.



- The file CLSettingHTDS.m  stores the 'ds' global variable (HTHomds, HTHSNds, HTHetds) used in HT continuation curves. This is stored as a setting in a settings object.  When the starter window is called this file will use the information to render the selection of SParam/UParam/T/eps1.   A 'gds' variable is stored in CLSettingHTDS.m (within the ds variable). This 'gds' variable is a variable that describes the state of the 'old gui', it is used to make HT-code believe it is working in  the older version of the GUI.


- The files SimConf_ConnectBase.m / SimConf_ConnectHom.m / SimConf_ConnectHet.m / SimConf_ConnectHSN.m, are used to register the 'simulation step' as a computation.  These files make use of a 'getPoint' file to compute the starting point of the simulation:

GUIConnect_getPointHet.m
GUIConnect_getPointHom.m
GUIConnect_getPointHSN.m

These files are used to compute the initial point (plus metadata). The code has been transfered from gui_Connect[Hom|Het|Hsn].m (oldgui)




- When a simulation is selected as HT node, the following code is executed:
GUIselectHTHet.m
GUIselectHTHom.m
GUIselectHTHSN.m

This code was migrated from the old gui.



- ContConf_HTHet.m ContConf_HTHom.m and ContConf_HTHSN.m are the HT-continuation configuration files. These files make use of migrated code to execute the init-function:
HTHet_performInit.m
HTHom_performInit.m
HTHSN_performInit.m



- When a HT-continuation is finished, the modified ContConf-solution is created, this modified ContConf-solution makes sure that a 'pointload' function is called whenever a HT-node is selected on that curve:
HTCurve.m
HTHomCurve.m
HTHSNCurve.m
HTHetCurve.m


- The pointloaders contain migrated code:
HTHet_pointload.m
HTHom_pointload.m
HTHSN_pointload.m


- The interpretation Files for HT must be able to predict the next HT-phase (described by 'index' in the ds-struct), this is constructed from the original init-files:
HTHet_nextIndex.m
HTHom_nextIndex.m
HTHSN_nextIndex.m


