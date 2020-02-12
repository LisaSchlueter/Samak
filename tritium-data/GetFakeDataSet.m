function DataSet = GetFakeDataSet(FakeStudyName)
if contains(FakeStudyName,'KNM1')
    DataSet = 'Knm1';
elseif contains(FakeStudyName,'KNM2') || contains(FakeStudyName,'Final')
    DataSet = 'Knm2';  
end
end