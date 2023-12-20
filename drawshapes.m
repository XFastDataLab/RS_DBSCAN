% Points 是数据集
% class 是每个数据的类别标签
% ND 是数据大小
% shape_ids  saves the shape id of each point 
% colors_rgb save the rgb infor of each  point
function [shape_ids, colors_rgb, label_idx]=drawshapes(points,class,ND)
    %颜色变换步长
    colorindent = 100/7;
    %形状
    shapes='*o^+sx<dp.^sh>dv';
    shape_ids=zeros(ND,1);
    colors_rgb=zeros(ND,3);
    label_idx=zeros(ND,1);
    
    hold on;
    for i=1:ND
        if (class(i)>0)
            if(class(i)==1)
               v = [0 0 1];
               shapeindex= mod(class(i),14)+1;
               label=1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==2)
               v = [0 1 0];
               shapeindex= mod(class(i),14)+1;
               label=2;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==3)
               v = [1 0 0];
               shapeindex= mod(class(i),14)+1;
               label=3;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==4)
               v = [0 1 1];
               shapeindex= mod(class(i),14)+1;
               label=4;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==5)
               v = [1 0.84314 0];
               shapeindex= mod(class(i),14)+1;
               label=5;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==6)
               v = [1 0 1];
               shapeindex= mod(class(i),14)+1;
               label=6;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==7)
               v = [0.4660 0.6740 0.1880];
               shapeindex= mod(class(i),14)+1;
               label=7;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==8)
                v = [1 0.64706 0];
               shapeindex= mod(class(i),14)+1;
               label=8;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            elseif(class(i)==9)
               v = [0.125 0.698 0.6663];
               label=9;
               shapeindex= mod(class(i),14)+1;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
             elseif(class(i)==10)
               v = [0.803 0.521 0.247];
               shapeindex= mod(class(i),14)+1;
               label=10;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
             elseif(class(i)==11)
               v = [0.803 0.3607 0.36078];
               shapeindex= mod(class(i),14)+1;
               label=11;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
             elseif(class(i)==12)
               v = [1 0.411 0.705];
               shapeindex= mod(class(i),14)+1;
               label=12;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
             elseif(class(i)==13)
               v = [0.7411 0.71765 0.41961];
               shapeindex= mod(class(i),14)+1;
               label=13;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
             elseif(class(i)==14)
               v = [0.094 0.4549 0.803];
               shapeindex= mod(class(i),14)+1;
               label=14;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            else(class(i)==15)
               v = [0 0.7440 1];
               shapeindex= mod(class(i),14)+1;
               label=15;
               %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
               plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
               set(gca,'FontSize',15);
            end

        else
            %离群点用红色o画出
            plot(points(i,1),points(i,2),'.','MarkerSize',3,'MarkerEdgeColor','w');  
            v = [1 0 0];
            shapeindex= 10;
            label=-1;
            %set(gcf,'unit','normalized','position',[0.2,0.2,0.1,0.15]);
            plot(points(i,1),points(i,2),shapes(shapeindex),'MarkerSize',5,'MarkerFaceColor','w','MarkerEdgeColor',v); 
            set(gca,'FontSize',15);
        end
        shape_ids(i)=shapeindex;
        colors_rgb(i,:)=v;
        label_idx(i)=label;
        
    end 
end



