classdef FacadeMatrix
    
    properties
        vector
        dim
    end
    methods
        function obj = FacadeMatrix(matrix, index)
            
            obj.vector = matrix(:, index);
            obj.dim = size(matrix);
            
        end
        function out = subsref(obj , S)
            out = obj.vector(S.subs{1});
        end
       
        
        function [m, n] = size(obj, varargin)
                if ~isempty(varargin)
                    m = obj.dim(varargin{1});
                else
                    m = obj.dim(1);
                    n = obj.dim(2);
                end
        end
        
    end
    
    
end