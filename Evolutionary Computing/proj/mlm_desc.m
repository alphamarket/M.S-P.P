classdef mlm_desc
% MLM_DESC
% @brief    The description container of XGroups
    properties(GetAccess = public, SetAccess = protected)
        clss
        desc
    end
    methods
        function this = mlm_desc(clss, desc)
%           MLM_DESC
%           @brief          Construct a description
%           @param  clss    The class of description
%           @param  desc    The description
%           @return         A description handler
            this.clss  = clss;
            this.desc = desc;
        end
        function ins = instantize(this, i, max, genNO)
%           INSTANTIZE
%           @brief      instantize a new instance due to current description
%           @param  i   The index of instance
%           @param  max The max# of instances
%           @return     The generated new instance
            switch this.clss
                case 'nn'
                    ins = this.desc(i * randn / (1 * max * sqrt(genNO))); 
                otherwise
                    error('description class `%s` is not defined!', this.clss); 
            end
        end
    end
end