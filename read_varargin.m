function [transform,err_est_type,q] = read_varargin(varargin)
transform = [];
err_est_type = 'default';
q = 0;
for i = 1:length(varargin)
    if isa(varargin{i}, 'char') || isa(varargin{i}, 'string')
        if ~strcmp(err_est_type, 'default')
            error('Cannot set error estimate type multiple times')
        end
        if ~strcmp(varargin{i}, 'jack')...
                && ~strcmp(varargin{i}, 'leave-one-out')
            error('Unrecognized estimate type %s', varargin{i})
        end
        err_est_type = varargin{i};
    elseif isa(varargin{i}, 'function_handle')
        if ~isempty(transform)
            error('Cannot set transformation twice')
        end
        transform = varargin{i};
    elseif isnumeric(varargin{i}) && ceil(varargin{i}) == varargin{i}
        q = varargin{i};
    else
        error('Optional argument not recognized')
    end
end

if ~isempty(transform) && strcmp(err_est_type, 'leave-one-out')
    error('Cannot use transform with leave-one-out estimate')
end

if strcmp(err_est_type, 'default')
    if isempty(transform)
        err_est_type = 'leave-one-out';
    else
        err_est_type = 'jack';
    end
end
end