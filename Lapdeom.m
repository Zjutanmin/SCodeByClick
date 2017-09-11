% save('LapSave', 'setting', 'data_fea', 'data_label')
load('LapSave')

% load('LapSave', 'data_fea');
Lap = [];
% if strcmp(cmethodorg, 'MLR')
    if ~isempty(setting.featparaW)
% setting.PCAenergy = 0;
        if setting.PCAenergy
            data_fea1{1} = GetPCAFeature(setting.PCAMethod, setting.featsize, ...
                setting.LocalPCA, setting.NormFeaNF, setting.NormFea1, ...
                [data_fea, data_label], setting.PCACOEFF, setting.FeatConf);
            
        else
            data_fea1{1} = data_fea;
        end
        PCAenergy{1} =  setting.PCAenergy;pcasstr{1} =  setting.pcasstr;
        data_fea1{2} = data_fea;PCAenergy{2} =  0;pcasstr{2} =  '';
%         for jj = length(PCAenergy):1
        for jj = 1:length(PCAenergy)
            setting.pcasstr = pcasstr{jj};
            setting.PCAenergy =  PCAenergy{jj};
 
        method = {'N','P'};
        c_k_cluster_I = {[10, 20, 50, 100], [5, 10, 50, 100, 500, 1000]};
        Mrange = 1:length(method);
        for i_method = Mrange
            k_cluster_I = c_k_cluster_I{i_method};
            Qrange = 1:length(k_cluster_I);
            for i_k_cluster_I = Qrange;
                setting.featparaW{4} = method{i_method};
                setting.featparaW{5} = k_cluster_I(i_k_cluster_I) ;
                if ~isfield(setting, 'pcasstr')
                    setting.pcasstr = '';
                end
                LapName = [setting.config_file, '-' setting.pcasstr...
                    'W' num2str(setting.featpara{1}) 'N' num2str(setting.NormFea)];
                try type = setting.featparaW{4}; catch type = 'N'; end
                try K = setting.featparaW{5}; catch K = 20; end
                Lap = GetLaplace(data_fea1{jj}, LapName, type, K);
            end
        end
        end
        clear 'data_fea1';
    end
% end