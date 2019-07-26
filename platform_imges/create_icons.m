%% psha
handles.Exit_button.CData     = double(imread('exit.jpg'))/255;
handles.setXYLimits.CData     = double(imresize(imread('Ruler.jpg'),[16 16]))/255;
handles.Distance_button.CData = double(imread('Scale.jpg'))/255;
handles.po_refresh_GE.CData   = double(imread('Refresh.jpg'))/255;
plus   =double(imread('plus.jpg'))/255;
exitmat=double(imread('exit.jpg'))/255;
Rulermat = double(imresize(imread('Ruler.jpg'),[16 16]))/255;
Scalemat = double(imread('Scale.jpg'))/255;
Refreshmat = double(imread('Refresh.jpg'))/255;
Refresh2   = double(imread('Refresh2.jpg'))/255;
Engine = double(imread('engine.jpg'))/255;
Leg    = double(imread('Legend.jpg'))/255;
Limits = double(imread('Limits.jpg'))/255;
form1 = double(imread('form1.jpg'))/255;
form2 = double(imread('form2.jpg'))/255;
die   = double(imread('die.jpg'))/255;
Undock = double(imread('Undock.jpg'))/255;
Play   = double(imread('Play.jpg'))/255;
MagicWand = double(imresize(imread('MagicWand.jpg'),[16 16]))/255;
Legend = double(imread('Legend.jpg'))/255;
save PSHAbuttons exitmat Rulermat Scalemat Refreshmat Engine Leg Limits form1 form2 die Refresh2 plus Undock Play MagicWand Legend

%% scenarios
Undockbutton = double(imread('Undock.jpg'))/255;
sumbuttonbutton = double(imread('Sum.jpg'))/255;
ExitButtonbutton = double(imread('exit.jpg'))/255;
plotTesselbutton=double(imresize(imread('Tessel.jpg'),[20 20]))/255;
GridManagerbutton=double(imread('AxisGrid.jpg'))/255;
ShowNodesbutton=double(imread('Nodes.jpg'))/255;
AddNewScenariobutton=double(imread('AddNew.jpg'))/255;
DeleteButtonbutton=double(imread('selection_delete.jpg'))/255;
DiscretizeMButtonbutton=double(imread('MagnitudeDisc.jpg'))/255;
DiscretizeRButtonbutton=double(imread('Grid.jpg'))/255;
SortButtonbutton =double(imread('Sort.jpg'))/255;
Rulerbutton = double(imread('Scale.jpg'))/255;
invokeWizbutton = double(imresize(imread('MagicWand.jpg'),[20 20]))/255;
save All_Scenario_Buttons