clear all
clc

nodes=dlmread('Job-1.inp');
elements_abaqus=dlmread('Job-2.inp');

elements(:,1:2)=elements_abaqus(:,1:2);
elements(:,3)=elements_abaqus(:,6);
elements(:,4)=elements_abaqus(:,7);
elements(:,5)=elements_abaqus(:,3);
elements(:,6)=elements_abaqus(:,5);
elements(:,7)=elements_abaqus(:,9);
elements(:,8)=elements_abaqus(:,8);
elements(:,9)=elements_abaqus(:,4);

loaded_nodes=0;

list_restrained=zeros(length(nodes(:,1)),4);
k=1;

for i=1:length(nodes(:,1))
    % Restrict the lower face
    if ((nodes(i,3)==min(nodes(:,3)))&&(nodes(i,2)==min(nodes(:,2)))&&(nodes(i,4)==min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=1;
        list_restrained(k,3)=1;
        list_restrained(k,4)=1;
        k=k+1;
    end
    
    if ((nodes(i,3)==min(nodes(:,3)))&&(nodes(i,2)==min(nodes(:,2)))&&(nodes(i,4)~=min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=1;
        list_restrained(k,3)=1;
        list_restrained(k,4)=0;
        k=k+1;
    end
    
    if ((nodes(i,3)==min(nodes(:,3)))&&(nodes(i,2)~=min(nodes(:,2)))&&(nodes(i,4)==min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=0;
        list_restrained(k,3)=1;
        list_restrained(k,4)=1;
        k=k+1;
    end
    
    if ((nodes(i,3)==min(nodes(:,3)))&&(nodes(i,2)~=min(nodes(:,2)))&&(nodes(i,4)~=min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=0;
        list_restrained(k,3)=1;
        list_restrained(k,4)=0;
        k=k+1;
    end
    
    % Restrict the top face
    if ((nodes(i,3)==max(nodes(:,3)))&&(nodes(i,2)==min(nodes(:,2)))&&(nodes(i,4)==min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=1;
        list_restrained(k,3)=1;
        list_restrained(k,4)=1;
        k=k+1;
    end
    
    if ((nodes(i,3)==max(nodes(:,3)))&&(nodes(i,2)==min(nodes(:,2)))&&(nodes(i,4)~=min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=1;
        list_restrained(k,3)=1;
        list_restrained(k,4)=0;
        k=k+1;
    end
    
    if ((nodes(i,3)==max(nodes(:,3)))&&(nodes(i,2)~=min(nodes(:,2)))&&(nodes(i,4)==min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=0;
        list_restrained(k,3)=1;
        list_restrained(k,4)=1;
        k=k+1;
    end
    
    if ((nodes(i,3)==max(nodes(:,3)))&&(nodes(i,2)~=min(nodes(:,2)))&&(nodes(i,4)~=min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=0;
        list_restrained(k,3)=1;
        list_restrained(k,4)=0;
        k=k+1;
    end    
    
    % Restrict the x side
    if ((nodes(i,3)~=max(nodes(:,3)))&&(nodes(i,3)~=min(nodes(:,3)))&&(nodes(i,2)==min(nodes(:,2)))&&(nodes(i,4)~=min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=1;
        list_restrained(k,3)=0;
        list_restrained(k,4)=0;
        k=k+1;
    end
      
    % Restrict the z side
    if ((nodes(i,3)~=max(nodes(:,3)))&&(nodes(i,3)~=min(nodes(:,3)))&&(nodes(i,2)~=min(nodes(:,2)))&&(nodes(i,4)==min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=0;
        list_restrained(k,3)=0;
        list_restrained(k,4)=1;
        k=k+1;
    end
    
    % Restrict the central axis
    if ((nodes(i,3)~=max(nodes(:,3)))&&(nodes(i,3)~=min(nodes(:,3)))&&(nodes(i,2)==min(nodes(:,2)))&&(nodes(i,4)==min(nodes(:,4))))
        list_restrained(k,1)=nodes(i,1);
        list_restrained(k,2)=1;
        list_restrained(k,3)=0;
        list_restrained(k,4)=1;
        k=k+1;
    end
end

nr=k-1;

list_restrained(k:length(list_restrained(:,1)),:)=[];

list_fixed=zeros(length(nodes(:,1)),3);
k=1;

for i=1:length(nodes(:,1))
    if (nodes(i,3)==max(nodes(:,3)))
        list_fixed(k,1)=nodes(i,1);
        list_fixed(k,2)=2;
        list_fixed(k,3)=7;
        k=k+1;
    end
end

fixed_nodes=k-1;

list_fixed(k:length(list_fixed(:,1)),:)=[];

limit=1000;
tol=1E-6;
tol2=1E-6;
output_step=1;
num_load=1;
e=1;
v=0.33;

% Open the .dat file and write in it
fid=fopen('input.dat','w');

fprintf(fid,'%d %d %d %d %d 8\n',length(elements(:,1)),length(nodes(:,1)),...
    nr,loaded_nodes,fixed_nodes);
fprintf(fid,'%d %f %f %f\n',limit,tol,e,v);
fprintf(fid,'8\n');
fprintf(fid,'%d %d\n',num_load,output_step);
fprintf(fid,'%f ',tol2);

fclose(fid);

% Open the .d file and write in it
fid=fopen('input.d','w');
fprintf(fid,'*THREE_DIMENSIONAL\n');
fprintf(fid,'*NODES\n');
for i=1:length(nodes(:,1))
    fprintf(fid,'%d %f %f %f\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
end
fprintf(fid,'*ELEMENTS\n');
for i=1:length(elements(:,1))
    fprintf(fid,'%d 3 8 1 %d %d %d %d %d %d %d %d 1\n',elements(i,1),...
        elements(i,2),elements(i,3),elements(i,4),elements(i,5),...
        elements(i,6),elements(i,7),elements(i,8),elements(i,9));
end

fclose(fid);

% Open the .bnd file and write in it
fid=fopen('input.bnd','w');

fprintf(fid,'%d %d %d %d\n',list_restrained(1,1),list_restrained(1,2),list_restrained(1,3),list_restrained(1,4));
for i=2:length(list_restrained(:,1))
    fprintf(fid,'%d %d %d %d\n',list_restrained(i,1),list_restrained(i,2),...
        list_restrained(i,3),list_restrained(i,4));
end

fclose(fid);

% Open the .fix file and write in it
fid=fopen('input.fix','w');

fprintf(fid,'%d %d %f\n',list_fixed(1,1),list_fixed(1,2),list_fixed(1,3));
for i=2:length(list_fixed(:,1))
    fprintf(fid,'%d %d %f\n',list_fixed(i,1),list_fixed(i,2),list_fixed(i,3));
end

fclose(fid);
