% ft_defaults

%% extrakce pouze bod� na povrchu - nen� hotovo
 [k,v]=boundary(mesh.pos(:,1),mesh.pos(:,2),mesh.pos(:,3));
 plot3(mesh.pos(:,1),mesh.pos(:,2),mesh.pos(:,3),'o')