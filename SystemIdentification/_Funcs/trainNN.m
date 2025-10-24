function varargout = trainNN(saveFileNN,epochs)

seqLength = 10;

disp('Training Neural Network')
    cd('NeuralNetwork\')
    hist = py.main_NN_V4.train_model('../_IDdata\NN\NNdata.csv', ['../',saveFileNN], int32(seqLength), int32(epochs));
    cd('../')

    if nargout>0
        varargout{1} = double(hist.history{'loss'});
    end

disp('Neural Network training done')

end