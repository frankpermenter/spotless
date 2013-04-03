classdef SpotProgPreProc 
    methods (Abstract)
        [sdpout,G,h,log] = preProcess(sdpin);
    end
end