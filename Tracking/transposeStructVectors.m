function str = transposeStructVectors(str,type)
if nargin<2 || isempty(type)
    type = 'row';
end
strFields = fields(str);
nStr = length(str);
for i = 1:nStr
    for sf = 1:length(strFields)
        if isstruct(str(i).(strFields{sf}))
            str(i).(strFields{sf}) = transposeStructVectors(str(i).(strFields{sf}),type);
        else
            if ~ischar(str(i).(strFields{sf})) || ismember(1,size(str(i).(strFields{sf})));
                switch type
                    case 'row'
                        if ~isrow(str(i).(strFields{sf}))
                            str(i).(strFields{sf}) = str(i).(strFields{sf})';
                        end
                    case 'col'
                        if isrow(str(i).(strFields{sf}))
                            str(i).(strFields{sf}) = str(i).(strFields{sf})';
                        end
                end
            end
        end
    end
end
