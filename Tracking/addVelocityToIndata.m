function [indata speed] = addVelocityToIndata(indata)
indata = transposeStructVectors(indata,'col');
n_sess = size(indata,1);
n_seg = size(indata,2);
speed = nan(n_sess,n_seg);

for ii = 1:n_sess
    for kk = 1:n_seg
        [v s] = writeVelocity.getVelocity(indata(ii,kk).x,indata(ii,kk).y,indata(ii,kk).t);
        indata(ii,kk).v = v;
        speed(ii,kk) = s;
        disp('Average speed:'); disp(s);
    end
end
