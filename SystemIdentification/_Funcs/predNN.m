function [pred,varargout] = predNN(saveFileNN,y_val,u_val)

seq_length = 10;
for i=1:(length(u_val)-seq_length-1)
    X_test(i,:,1) = u_val(i+1:i+seq_length+1);
    X_test(i,:,2) = y_val(i:i+seq_length);
end

    cd('NeuralNetwork\')
    pred = py.main_NN_V3.predict_model(['../',saveFileNN], X_test);
    cd('../')
    pred = [y_val(1:seq_length+1)' double(pred)];

if nargout>1
    varargout{1} = y_val'-pred;
end

end