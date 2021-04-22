function parsave1(matfilename,CELLS,CHANNELS)
%PARSAVE1(matfilename,CELLS,CHANNELS)
%Allows saving in parallel programming (otherwise transparency violation error).
save(matfilename,'CELLS','CHANNELS','-v6')
end