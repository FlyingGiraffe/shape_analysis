function mu = barycenter(mu1,mu2,alpha1,alpha2,kernel,H0,niter)
mu_shape = size(mu1);
v1 = ones(mu_shape); v2 = ones(mu_shape);
%w1 = ones(mu_shape); w2 = ones(mu_shape);
mu = ones(mu_shape);
threshold = 1e-3;
for k = 1:niter
    mu_old = mu;
    mu = ones(mu_shape);
    
    H_v1 = kernel(v1); H_v2 = kernel(v2);
    w1 = mu1./H_v1; w2 = mu2./H_v2;
    w1(H_v1==0) = 0; w2(H_v2==0) = 0;
    
    d1 = v1.*kernel(w1); d2 = v2.*kernel(w2);
    mu = mu.*(d1.^alpha1).*(d2.^alpha2);
    
    mu = entropic_sharpening(mu,H0);
    
    v1 = v1.*mu./d1; v2 = v2.*mu./d2;
    v1(d1==0) = 0; v2(d2==0) = 0;
    sum(abs(mu_old-mu),'all')
    if sum(abs(mu_old-mu),'all')<threshold && k>2
        break;
    end
end

function mu_new = entropic_sharpening(mu,H0)
if sum(mu,'all')-entropy(mu)>H0+1
    fun = @(b) entropy(mu.^b)-entropy(mu.^b)-(1+H0);
    beta = fzero(fun);
    if beta < 0
        beta = 1;
    end
else
    beta = 1;
end
mu_new = mu.^beta;

function H = entropy(mat)
log_mat = log(mat);
log_mat(mat<=0) = 0;
H = -sum(mat.*log_mat,'all');
