function DataSet = GetFakeDataSet(FakeStudyName)
if contains(FakeStudyName,'KNM1')
    DataSet = 'Knm1';
elseif contains(FakeStudyName,'KNM2')
    DataSet = 'Knm2';  
elseif contains(FakeStudyName,'DesignReport') || contains(FakeStudyName,'Final') || contains(FakeStudyName,'SensitivityLimit')
    DataSet = 'Fake';
end
end