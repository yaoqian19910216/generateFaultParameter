Nx = 621;
Nz = 621;
 
mud=zeros(Nx,Nz);
mus=zeros(Nx,Nz);
dx = 20;
dz = 20;
 
for i=1:Nx
    for k=1:Nz
        mud(i,k)=0.525;
        mus(i,k)=0.75;
    end
end
 
for i=1:dx+1
    mud(i,:)=0.525+(dx+1-i)*0.01;
    mus(i,:)=0.75+(dx+1-i)*0.01;
end
   
for i=Nx-dx:Nx
    mud(i,:)=0.525+(i-Nx+dx)*0.01;
    mus(i,:)=0.75+(i-Nx+dx)*0.01;
end
 
for i=Nz-dz:Nz
    mud(:,i)=0.525+(i-Nx+dx)*0.01;
    mus(:,i)=0.75+(i-Nx+dx)*0.01;
end
 
% 
% for i=1:Nx
%     for k=1:Nz
%         if((i-161)^2+(k-161)^2>=140^2)
%         mud(i,k)=0.525+(sqrt((i-161)^2+(k-161)^2)-140)*0.01;
%         mus(i,k)=0.75+(sqrt((i-161)^2+(k-161)^2)-140)*0.01;
%         end
%     end
% end
 
 
fid = fopen('mud','w'); fwrite(fid, mud, 'single'); fclose(fid);
fid = fopen('mus','w'); fwrite(fid, mus, 'single'); fclose(fid);