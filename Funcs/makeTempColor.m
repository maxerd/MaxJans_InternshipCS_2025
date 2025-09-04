function tempVal = makeTempColor(temp,tempMax,tempMin)
    clr = zeros(1,3);

    tempDiff = tempMax-tempMin;
    % tempVal = ceil(((temp-tempMin)/tempDiff)*255);
    tempVal = ceil(((temp-tempMin)/tempDiff)*511);
    
    % clr(3) = 1-tempVal;
    % clr(2) = 0;
    % clr(1) = tempVal;
    tempVal = max(1,tempVal);
    tempVal = min(511,tempVal);

end