


dataDir_OL = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\matlabFiles\measurements\processedData\p__CB_none__ST_ms__SM_24w__SA_12w__DT_250829__MD_8h__WT_no__DS_10.mat'; % Open loop data used for identification

dat = load(dataDir_OL);

exportDat = [0 0;dat.Watt' dat.tempTM'];

csvwrite('file.csv',exportDat)

