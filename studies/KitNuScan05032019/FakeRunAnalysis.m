FR1=RunAnalysis('DataType','Fake','FakeRunType','FakeKITNuScan','RunNr',4,'exclDataStart',2);
FR1.Fit; 
FR1.PlotFit('Mode','Rate');
%FR1.PlotDataModel_KNM1;
