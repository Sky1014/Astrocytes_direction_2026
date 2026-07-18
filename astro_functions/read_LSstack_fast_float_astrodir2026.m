function [stack,dim]=read_LSstack_fast_float_astrodir2026(fpath,dim)
%% reads whole stack

info=dir(fpath);
if (isempty(info))
    error('File Does not exist');
end

stack=read_LSstack_fast_float_mex64_astrodir2026(fpath,int32(dim(1:2)));

size(stack)

stack=reshape(stack,dim);
