function saveDataToNN(y,u)

    try
        exportDat = [0 0;u y];
    catch
        error('At least one of the inputs should be transport')
    end
    
    csvwrite('_IDdata\NN\NNdata.csv',exportDat)

end
