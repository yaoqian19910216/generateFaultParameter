lx = 321;       % number of cells lx by lx
lh = lx/2;      % wavenumber for Nyquist frequency
dx=100;
Nx=321;
Ny=101;
Nz=321;
transf = zeros(lx, lx); 
 
hurst = 1.0; % Hurst exponent, for self-similar stress
exph = hurst + 1.; 
 
% Generate a self-affine random function, 
%   lowest mode being one wavelength in x-dir with unit amplitude
ncut = 2; kcut = lh/ncut; % % high-cut filter at a factor below Nyquist frequency 
npole = 4; % was 4  
fcsq = kcut*kcut;  % square of filter cutoff wavenumber
 
% rand('state', 6444);
% rand('state', 644);
 
% === Initialize array 'transf' (function in frequency domain) ===

rand('state',1);
for j = 1:lx
    for i = 1:lx
        % Generate a pair of independent normal random variables
        % using the Ploar Method
        S = 2.;
        while S > 1.
            a=rand(1,2);
            vtmp = 2*a - 1; S = vtmp(1)^2 + vtmp(2)^2;
        end
        transf(i,j) = sqrt(-2*log(S)/S)*sqrt(2.)/2.*complex(vtmp(1),vtmp(2));
    end
end
 
% === Modify array 'transf' in lowest modes, and specify power spectrum ===
for j = 1:lx
    for i = 1:lx
        ik = i-1;
        jk = j-1;
        if (ik > lh); ik = ik-lx; end
        if (jk > lh); jk = jk-lx; end
        ksq = ik^2 + jk^2; 
        fsq = ksq; 
        if ksq == 0
            % Mode (0.0) has zero amplitude
            transf(i,j) = 0.;  
            pspec = 0.;
        elseif ksq == 1
            if jk  == 0
                % Modes (1,0) and (-1,0) have double amplitude
%                transf(i,j) = -2.;
                pspec = 1./fsq^exph;
            else
                % Modes (0,1) and (0,-1) have zero amplitude
%                transf(i,j) = 0.; 
                pspec = 1./fsq^exph;
            end
        elseif ksq == 2
            % Mode (1,1) has zero amplitude
%            transf(i,j) = 0.;
            pspec = 1/fsq^exph;
        else
            % Highers modes have random phase
            pspec = 1/fsq^exph;
        end
        
        % Multiply transform by desired power spectrum (pspec) 
        %  and apply low-pass filter
        filter = 1./(1.+(fsq/fcsq)^npole);
        transf(i,j) = sqrt(pspec*filter)* transf(i,j);
    end
end
 
% Apply inverse fourier transform
w = real(ifft2(transf)); % 0.25*
w = w/max(abs(w(:)));
 
figure; imagesc(w'); xlabel('Lx (strike)'); ylabel('Lz (depth)');

mean=sum(sum(w))/lx/lx;
h_rms=sqrt(sum(sum((w-mean).^2))/(Nx*Nz));
ratio=log10(h_rms/lx/dx);

for i=1:Nx
    for k=1:Nz
        for j=1:(Ny+1)/2;
            x(i,j,k)=(i-1)*dx;
            y(i,j,k)=(j-1)*dx;
            z(i,j,k)=(k-1)*dx;
        end
        for j=(Ny+1)/2+1:Ny+1;
            x(i,j,k)=(i-1)*dx;
            y(i,j,k)=(j-2)*dx;
            z(i,j,k)=(k-1)*dx;
        end
    end
end

imagesc(w)

w1=zeros(Nx,1);
for i=191:211
    w1(i)=-(sin((i-161)/10*pi+pi/2)-1);
end

% w1=zeros(Nx,Nz);
% gaussian_nx=81;
% for i=1:gaussian_nx
%     for j=1:gaussian_nx
%         w1((Nx+1)/2+i,(Nx+1)/2+j)=-exp(-((i-(gaussian_nx+1)/2)^2+(j-(gaussian_nx+1)/2)^2)/300);
%         w1((Nx+1)/2-i,(Nx+1)/2+j)=exp(-((i-(gaussian_nx+1)/2)^2+(j-(gaussian_nx+1)/2)^2)/300);
%         w1((Nx+1)/2-i,(Nx+1)/2-j)=-exp(-((i-(gaussian_nx+1)/2)^2+(j-(gaussian_nx+1)/2)^2)/300);
%         w1((Nx+1)/2+i,(Nx+1)/2-j)=exp(-((i-(gaussian_nx+1)/2)^2+(j-(gaussian_nx+1)/2)^2)/300);
%     end
% end
    
%fid = fopen('Y','w'); fwrite(fid, y, 'single'); fclose(fid);
for i=1:Nx
    for k=1:Nz
        y(i,(Ny+1)/2-1:(Ny+1)/2+2,k)=y(i,(Ny+1)/2-1:(Ny+1)/2+2,k)+500*w1(k);%first i, rough along x, then k
    end
end

npml=10;
% interpolation
bg = npml+2; ed = (Ny+1)/2-1;
h = 1/(ed-bg+2);
for i = bg:ed
    y(:,i,:) = y(:,bg-1,:)+(y(:,ed+1,:)-y(:,bg-1,:))*(i-bg+1)*h;
end
bg = (Ny+1)/2+2; ed = Ny-npml-1;
h = 1/(ed-bg+2);
for i = bg:ed
    y(:,i,:) = y(:,bg-1,:)+(y(:,ed+1,:)-y(:,bg-1,:))*(i-bg+1)*h;
end

 plot(x(:,:,end),y(:,:,end))
% 
% %shear traction
% ts=zeros(lx,Nz+1,lx);
% for i=1:321
%     for k=1:171
%         ts(i,k)=70*10^6;
%         if((i-161)^2+(k-86)^2<=15^2)
%            ts(i,k)=81.6*10^6;
%         end
%     end
% end
% 
% for i=(lx+1)/2-1500/dx:(lx+1)/2+1500/dx
%     for k=(lx+1)/2-1500/dx:(lx+1)/2+1500/dx
%         ts(i,k)=81.60*10^6;
%     end
% end

%% calculate spectrum
fid = fopen('x','w'); fwrite(fid, x, 'single'); fclose(fid);
fid = fopen('y','w'); fwrite(fid, y, 'single'); fclose(fid);
fid = fopen('z','w'); fwrite(fid, z, 'single'); fclose(fid);
% fid = fopen('ts','w'); fwrite(fid, ts, 'single'); fclose(fid);
a = 100 * w(181,:);
N = length(a);
xdft = fft(a);
xdft = xdft(1:(N+1)/2);
psdx = (1/(N*N)) * abs(xdft).^2;
psdx(2:end-1) = 2 * psdx(2:end-1);
freq = (0:(N-1)/2) / N / dx;
Nnew = 80;
p = polyfit(log10(freq(2:end-Nnew)),log10(psdx(2:end-Nnew)),1);
yfit = polyval(p,log10(freq(2:end-Nnew)));
plot(log10(freq),log10(psdx));hold on;plot(log10(freq(2:end-Nnew)),yfit,'k','linewidth',2);

mean=sum(w(:,200))/321;
h_rms=sqrt(sum((w(:,200)-mean).^2)/Nx);
ratio=log10(h_rms/lx/dx);

%% plot profile 
for i=1:Nx
    for j=1:Nz
        wnew(i,j)=w(100,j);
    end
end
surf(wnew*500);hold on
scatter3(161,161,wnew(161,161)*500,'rp','Linewidth',100)
set(gca,'fontsize',20);title('Rough fault')
view(90,280);
