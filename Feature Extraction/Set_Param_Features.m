function Params_Features = Set_Param_Features(Feat_Type)
%Setting parameters for the feature extraction.

if strcmp(Feat_Type,'SIFT'),
    Params_Features.Type = 0;
    Params_Features.peak = 0.5;
    Params_Features.edge = 800;
else if strcmp(Feat_Type,'SURF'),
        Params_Features.Type = 1;
        Params_Features.Options.upright = false;
        Params_Features.Options.tresh = 0.0001;
     else
         Params_Features.Type = 2;
         fpath = fileparts(which(mfilename));
         Params_Features.exec_str = ['"' fpath '\demo_ASIFT"'];
    end
end