% 3D figure as example for FLY.M. 
% Source: http://www.mathworks.co.uk/help/matlab/visualize/example-moving-the-camera-through-a-scene.html

load wind
wind_speed = sqrt(u.^2 + v.^2 + w.^2);
figure
hpatch = patch(isosurface(x,y,z,wind_speed,35));
isonormals(x,y,z,wind_speed,hpatch)
set(hpatch,'FaceColor',[0.75,0.25,0.25],'EdgeColor',[0.6,0.4,0.4]);
[f vt] = reducepatch(isosurface(x,y,z,wind_speed,45),0.05); 

daspect([1,1,1]);
hcone = coneplot(x,y,z,u,v,w,vt(:,1),vt(:,2),vt(:,3),2);
set(hcone,'FaceColor','blue','EdgeColor','none');
camproj perspective 
hlight = camlight('headlight'); 
set(hpatch,'AmbientStrength',.1,...
      'SpecularStrength',1,...
      'DiffuseStrength',1);
set(hcone,'SpecularStrength',1);
set(gcf,'Color','k')
set(gca,'Color',[0,0,0.25])
lighting gouraud
set(gcf,'Renderer','opengl')
axis off
