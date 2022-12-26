map=ones(250,250);
disp_map=ones(250,250);
st_pt=[10,240];
end_pt=[200,50];
disp_map=add_obs_disp(disp_map,[45,220],20);
disp_map=add_obs_disp(disp_map,[30,130],20);
disp_map=add_obs_disp(disp_map,[150,100],30);

map=add_obs(map,[45,220],20);
map=add_obs(map,[30,130],20);
map=add_obs(map,[150,100],30);


route=astar(map,st_pt,end_pt);

cmap = [1 1 1; ...
    0 0 0; ...
    1 0 0; ...
    0 0 1; ...
    0 1 0; ...
    1 1 0; ...
    0.5 0.5 0.5];
colormap(cmap);
for k = 2:length(route) - 1  
    value=sub2ind(size(map),route(k,1),route(k,2));
    disp_map(value) = 7;
    pause(0.1);
    image(1.5,1.5, disp_map);
    grid on;
    axis image;
end
function map=add_obs(map,centre,radius)
    [x,y]=size(map);
    i=1;
    j=1;
    while(i<x)
        j=1;
        while(j<y)
            if(norm([i-centre(1),j-centre(2)])<radius)
                map(i,j)=2;
            end
        j=j+1;
        end
        i=i+1;
    end
end
function map=add_obs_disp(map,centre,radius)
    [x,y]=size(map);
    i=1;
    j=1;
    while(i<x)
        j=1;
        while(j<y)
            if(norm([i-centre(1),j-centre(2)])<radius)
                map(j,i)=2;
            end
        j=j+1;
        end
        i=i+1;
    end
end
function route=astar(map,start,end_pt)
    st_point=sub2ind(size(map),start(1),start(2));
    en_point=sub2ind(size(map),end_pt(1),end_pt(2));
    map(st_point)=5;
    map(en_point)=6;
    [nrows,ncols]=size(map);
    parent=zeros(nrows,ncols);
    map_g=Inf(nrows,ncols);
    map_f=Inf(nrows,ncols);
    map_f(st_point)=heuristic(start,end_pt);
    disp(map_f(st_point));
    map_g(st_point)=0;
    while true
        [min_f, current] = min(map_f(:));
        if ((current == en_point) || isinf(min_f))
            break;
        end
        if(map(current)==2)
            continue;
        end
        map(current) = 3;
        map_f(current) = Inf;
        [i, j] = ind2sub(size(map_f), current);
        if(i<nrows)
            next_node=sub2ind(size(map),i+1,j);
            if((map(next_node)==1 || map(next_node)==4 || map(next_node)==6) && map_g(next_node)>map_g(current)+1)
                if(map(next_node)==1)
                    map(next_node)=4;
                end
                map_g(next_node)=map_g(current)+1;
                map_f(next_node)=min_f-heuristic([i,j],end_pt)+heuristic([i+1,j],end_pt);
                parent(next_node)=current;
            end
        end
        
        if(i<ncols)
            next_node=sub2ind(size(map),i,j+1);
            if((map(next_node)==1 || map(next_node)==4 || map(next_node)==6) && map_g(next_node)>map_g(current)+1)
                if(map(next_node)==1)
                    map(next_node)=4;
                end
                map_g(next_node)=map_g(current)+1;
                map_f(next_node)=min_f-heuristic([i,j],end_pt)+heuristic([i,j+1],end_pt);
                parent(next_node)=current;
            end
        end

        if(i>1)
            next_node=sub2ind(size(map),i-1,j);
            if((map(next_node)==1 || map(next_node)==4 || map(next_node)==6) && map_g(next_node)>map_g(current)+1)
                if(map(next_node)==1)
                    map(next_node)=4;
                end
                map_g(next_node)=map_g(current)+1;
                map_f(next_node)=min_f-heuristic([i,j],end_pt)+heuristic([i-1,j],end_pt);
                parent(next_node)=current;
            end
        end
        
        if(j>1)
            next_node=sub2ind(size(map),i,j-1);
            if((map(next_node)==1 || map(next_node)==4 || map(next_node)==6) && map_g(next_node)>map_g(current)+1)
                if(map(next_node)==1)
                    map(next_node)=4;
                end
                map_g(next_node)=map_g(current)+1;
                map_f(next_node)=min_f-heuristic([i,j],end_pt)+heuristic([i,j-1],end_pt);
                parent(next_node)=current;
            end
        end
    end  
    x_route=end_pt(2);
    y_route=end_pt(1);
    x=x_route;
    y=y_route;
    while (parent(y,x)~=0)
        [y,x]=ind2sub(size(map),parent(y,x));
        x_route=[x,x_route];
        y_route=[y,y_route];
    end
    route=[x_route;y_route];
    route=transpose(route);
end
function h_val=heuristic(co_ord,end_point)
    h_val=norm([end_point(1)-co_ord(1),end_point(2)-co_ord(2)]);
    h_val=h_val^0.5;
    %h_val=abs(end_point(1)-co_ord(1))+abs(end_point(2)-co_ord(2));
end