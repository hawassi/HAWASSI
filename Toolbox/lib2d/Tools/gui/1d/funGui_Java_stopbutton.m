function [jb]=funGui_Java_stopbutton(statusbarObj)
global Idstop
Idstop=0;
jb = javax.swing.JButton('Stop');
statusbarObj.add(jb,'East');

