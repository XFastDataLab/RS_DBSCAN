
function drawshapes3d(points,class,ND)
    %颜色变换步长
    colorindent = 100/7;
    %形状
    shapes='*o^+sdph>v*ox<.';
    hold on;
    for i=1:ND        
        if (class(i)>=0)
           if(class(i)==1)
               v = [0 0 1];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
           elseif(class(i)==0)
               v = [0.83 0.81 0.92];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==2)
               v = [0 1 0];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==3)
               v = [1 0 0];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==4)
               v = [0 1 1];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==5)
               v = [0.48235 0.40784 0.9333];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==6)
               v = [1 0 1];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==7)
               v = [0 0.80784 0.81963];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
           elseif(class(i)==8)
               v = [0.13333 0.5451 0.13333];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
           elseif(class(i)==9)
               v = [1 0.07843 0.57647];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
           elseif(class(i)==10)
               v = [0.69804 0.13333 0.13333];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==11)
               v = [0.98039 0.50196 0.44706];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==12)
               v = [0.62745 0.12549 0.94118];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==13)
               v = [0.41176 0.41176 0.41176];
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot3(points(i,1),points(i,2),points(i,3),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
           end
            view(3)
            
%             shape_ids(i)=shapeindex;
%             colors_rgb(i,:)=v;
             
            %colr,colg,colb是画图的RGB颜色值
%             col=colorindent*(class(i)-1)*600;
%             colr= mod(col,10);
%             col=fix(col/10);
%             colg=mod(col,10);
%             col=fix(col/10);
%             colb=mod(col,10);
            %v=0.1*[colr,colg,colb];
            %选画图形状
%             shapeindex= mod(class(i),14)+1;
%             set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15])
%             plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
%             set(gca,'FontSize',15);
%             %plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
        else
            %离群点用红色o画出
            plot3(points(i,1),points(i,2),points(i,3),'.','MarkerSize',3,'MarkerEdgeColor','r');  
            view(3)
        end
    end
    
%     ps= find(type==0);
%     
%     plot(points(ps,1),points(ps,2),'o','MarkerSize',5,'MarkerFaceColor','g','MarkerEdgeColor','g'); 
end
