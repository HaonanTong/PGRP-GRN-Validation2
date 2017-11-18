% f_main([1,2]);clc;close all;
% f_main([1,2,3]);clc;close all;
% f_main([1,2,3,4]);clc;close all;
% f_main([1,2,3,4,5]);clc;close all;
% f_main([5,6]);clc;close all;
%f_main([2,3,4,5]);clc;close all;
%f_main([3,4,5]);clc;close all;
%f_main([4,5]);clc;close all;
% target_12 = f_main([1,2]);clc;close all;

% trgt = f_main([1,2]);clc;close all;
% trgt = f_main([1,2,3]);clc;close all;
% trgt = f_main([1,2,3,4]);clc;close all;
% trgt = f_main([1,2,3,4,5]);clc;close all;
% trgt = f_main([5,6]);clc;close all;

tp_array_List = {[1,2],[1,2,3],[1,2,3,4],[1,2,3,4,5],[5,6]};

for i = 1 : length(tp_array_List)   
    tp_array = tp_array_List{i};
    ntp = length(tp_array);
    outputDir = 'Target_Analysis/Target_Analysis';

    for i = 1 : ntp
        outputDir = strcat(outputDir,'%d');
        outputDir = sprintf(outputDir,tp_array(i));
    end

    % Reconstruct Network
    trgt = f_main(tp_array);

    % Analyze and Record Score
    f_trgt_analysis(trgt,outputDir);

end