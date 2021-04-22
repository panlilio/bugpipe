classdef cellDataMP
    properties
        label = '0';
        frames = [0 0];
        divisionTime = NaN;
        parent = 0;
        lineID = 0;
        
        cellID = []; %Index in CELLS from segmentationMP for each frame
        intensityRange = [];
        intensityMeanSD = [];
        majorAxis = [];
        minorAxis = [];
        orientation = [];
        area = [];
        weightedCentroid = [];
        regionpropsCentroid = [];
        growthRate = NaN;
        QC = [];
        tAcq = [];
    end
    
    methods
        %Save cell properties to the relevant indices
        function obj = saveProps(obj,R)
            obj.intensityRange = [obj.intensityRange; R.intensityRange];
            obj.intensityMeanSD = [obj.intensityMeanSD; R.intensityMeanSD];
            obj.majorAxis = [obj.majorAxis R.majorAxis];
            obj.minorAxis = [obj.minorAxis R.minorAxis];
            obj.orientation = [obj.orientation R.orientation];
            obj.area = [obj.area R.area];
            obj.weightedCentroid = [obj.weightedCentroid; R.weightedCentroid];
            obj.regionpropsCentroid = [obj.regionpropsCentroid; R.regionpropsCentroid];
            obj.QC = [obj.QC; R.QC];
            obj.tAcq = [obj.tAcq; R.tAcq];
            if numel(obj.tAcq)>1
                obj.growthRate = [obj.growthRate; 1./diff(obj.tAcq([end-1 end]))*log(obj.majorAxis(end)/obj.majorAxis(end-1))];
            end
        end
    end
end