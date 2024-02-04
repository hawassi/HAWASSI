
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
jFrame    = get(handle(gcf),'JavaFrame');
try
    jRootPane = jFrame.fHG1Client.getWindow;
catch
    jRootPane = jFrame.fHG2Client.getWindow;    %>=2014b
end
statusbarObj = com.mathworks.mwswing.MJStatusBar;

% Add a progress-bar to left side of standard MJStatusBar container
jProgressBar = javax.swing.JProgressBar;
statusbarObj.add(jProgressBar,'West');  % 'West' => left of text; 'East' => right
% Beware: 'East' also works but doesn't resize automatically

% Set this container as the figure's status-bar
jRootPane.setStatusBar(statusbarObj);