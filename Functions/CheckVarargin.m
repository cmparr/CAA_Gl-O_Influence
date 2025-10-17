function a = CheckVarargin(x, varargin)
a = x;
if any(strcmp(varargin, 'thick'))
    disp('thick')
end
end