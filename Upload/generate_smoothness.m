function smooth = generate_smoothness(xs,ys,zs,len,wid,strike,num_grid,fault_type,fdip,sseg,sid,ssid,top,bot,side,sideid,plt,mode)
    % this function is used to generate smoothness for patches and
    % between segments.
    % mode == 1 adding 0th - order smoothness
    % mode == 2 adding 1st - order smoothness
    if plt == 1
        figure
    end
    smooth = [];
    Layer = 100;
    dpth(1) = Layer;
    % this loop will be used to differentiate fault and diff patches
    for j = 2:1:length(xs)
        if zs(j) == zs(j-1);
            dpth(j) = dpth(j-1);
        elseif zs(j) > zs(j-1);
            dpth(j) = dpth(j-1)+1;
        else
            Layer = Layer + 100;
            dpth(j) = Layer;
        end
    end
    
    N = Layer/100;  % the num of faults
    for j = 1:1:N
        ii = dpth >= j*100 & dpth < (j+1)*100;
        dpth(ii) = -(dpth(ii)-max(dpth(ii))-1);
    end
    
    id = dpth;  % patch index for each fault
    % get the corresponding index for every path on each fault
    % find the layer/patch index for each patch
    l = 0;
    for j = 1:1:length(num_grid)
        for k = 1:1:num_grid(j)
            l = l+1;
            in(l) = l;          % num index
            ip(l) = j;          % patch index
        end
    end
    
    N = sum(num_grid); % pathes number
    layer = max(dpth);  % top layer index  
    % horizontal smoothness within segment
    if fault_type(1) ~= 0
        for j = 1:1:length(in)-1
            if id(j+1) == id(j) && mode == 1 %&& ip(j+1) == ip(j);
                col = zeros(1,N*sum(fault_type));
                col(1,in(j)*sum(fault_type)+1) = -2/(len(j)+len(j+1));     % adding the first-order smoothing constraint
                col(1,(in(j)-1)*sum(fault_type)+1) = 2/(len(j)+len(j+1));   % on horizontal strike slip distribution 
                smooth = [smooth; col];
                col = zeros(1,N*sum(fault_type));
                col(1,in(j)*sum(fault_type)+2) = -2*fdip/(len(j)+len(j+1));     % adding the first-order smoothing constraint
                col(1,(in(j)-1)*sum(fault_type)+2) = 2*fdip/(len(j)+len(j+1));  % on horizontal dip slip distribution
                smooth = [smooth; col];
                if plt == 1
                    line([xs((in(j)+1));xs(in(j))],[ys((in(j)+1));ys(in(j))],[zs((in(j)+1));zs(in(j))],'Color','k'),hold on
                end
            end
        end
        for j = 1:1:length(in)-2
            if id(j+2) == id(j) && mode == 2 %&& ip(j+1) == ip(j);
                col = zeros(1,N*sum(fault_type));
                col(1,in(j+2)*sum(fault_type)-1) = -2*fdip/(len(j+2)+len(j+1));         % adding the second-order smoothing constraint
                col(1,in(j+1)*sum(fault_type)-1) = 8*fdip/(len(j+2)+2*len(j+1)+len(j));         % on horizontal strike slip distribution
                col(1,in(j)*sum(fault_type)-1) = -2*fdip/(len(j+1)+len(j));
                smooth = [smooth; col];
                col = zeros(1,N*sum(fault_type));
                col(1,in(j+2)*sum(fault_type)) = -2*fdip/(len(j+2)+len(j+1));         % adding the second-order smoothing constraint
                col(1,in(j+1)*sum(fault_type)) = 8*fdip/(len(j+2)+2*len(j+1)+len(j));
                col(1,(in(j))*sum(fault_type)) = -2*fdip/(len(j)+len(j+1));     % on horizontal dip slip distribution
                smooth = [smooth; col];
                if plt == 1
                    line([xs((in(j)+2))+300;xs((in(j)+1))],[ys((in(j)+2));ys((in(j)+1))],[zs((in(j)+2));zs((in(j)+1))],'Color','b'),hold on
                    line([xs((in(j)+1))-300;xs(in(j))],[ys((in(j)+1));ys(in(j))],[zs((in(j)+1));zs(in(j))],'Color','g')
                end
            end
        end
    end
    
    % vertical smoothness within segment
    if fault_type(2) ~= 0 && fdip ~= 0
        for j = 1:1:length(num_grid)
            % judge whether it has a contact with upper layer, since the
            % layer order is from bottom to top
            if mode == 2
                for k = layer:-1:3
                    i1 = find(id==k & ip==j);   % get the patch index with same depth (on the same fault )
                    i2 = find(id==k-1 & ip==j); % the lower patch index
                    i3 = find(id==k-2 & ip==j); % this would be added laplace smoothing
                    for l = 1:1:length(i1)
                        x11 = len(i1(1))*(l-1);
                        x12 = len(i1(1))*l;
                        for m = 1:1:length(i2)
                            x21 = len(i2(1))*(m-1);
                            x22 = len(i2(1))*m;
                            for n = 1:1:length(i3)
                                x31 = len(i3(1))*(n-1);
                                x32 = len(i3(1))*n;
                                if (x21 < x12 && x22 > x11) && (x31 < x22 && x32 > x21)       % has contact
                                    col = zeros(1,N*sum(fault_type));
                                    col(1,(in(i3(n))-1)*sum(fault_type)+1) = -2/(wid(i2(m))+wid(i3(n)));
                                    col(1,(in(i2(m))-1)*sum(fault_type)+1) = 8/(2*wid(i2(m))+wid(i1(l))+wid(i3(n)));
                                    col(1,(in(i1(l))-1)*sum(fault_type)+1) = -2/(wid(i2(m))+wid(i1(l)));
                                    smooth = [smooth; col];
                                    col(1,(in(i3(n))-1)*sum(fault_type)+2) = -2/(wid(i2(m))+wid(i3(n)));
                                    col(1,(in(i2(m))-1)*sum(fault_type)+2) = 8/(2*wid(i2(m))+wid(i1(l))+wid(i3(n)));
                                    col(1,(in(i1(l))-1)*sum(fault_type)+2) = -2/(wid(i2(m))+wid(i1(l)));
                                    smooth = [smooth; col];
                                    if plt == 1
                                        line([xs(in(i3(n)));xs(in(i2(m)))],[ys(in(i3(n)));ys(in(i2(m)))],[zs(in(i3(n)));zs(in(i2(m)))],'Color','c'),hold on
                                        line([xs(in(i2(m)));xs(in(i1(l)))],[ys(in(i2(m)));ys(in(i1(l)))],[zs(in(i2(m)));zs(in(i1(l)))],'Color','y')
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if mode == 1
                for k = layer:-1:2
                    i1 = find(id==k & ip==j);
                    i2 = find(id==k-1 & ip==j);
                    for l = 1:1:length(i1)
                        x11 = len(i1(1))*(l-1);
                        x12 = len(i1(1))*l;
                        for m = 1:1:length(i2)
                            x21 = len(i2(1))*(m-1);
                            x22 = len(i2(1))*m;
                            if x21 < x12 && x22 > x11       % has contact
                                col = zeros(1,N*sum(fault_type));
                                col(1,(in(i2(m))-1)*sum(fault_type)+1) = -2/(wid(i2(m))+wid(i1(l)));
                                col(1,(in(i1(l))-1)*sum(fault_type)+1) = 2/(wid(i2(m))+wid(i1(l)));
                                smooth = [smooth; col];
                                col = zeros(1,N*sum(fault_type));
                                col(1,(in(i2(m))-1)*sum(fault_type)+2) = -2*fdip/(wid(i2(m))+wid(i1(l)));
                                col(1,(in(i1(l))-1)*sum(fault_type)+2) = 2*fdip/(wid(i2(m))+wid(i1(l)));
                                smooth = [smooth; col];
                                if plt == 1
                                    line([xs(in(i2(m)));xs(in(i1(l)))],[ys(in(i2(m)));ys(in(i1(l)))],[zs(in(i2(m)));zs(in(i1(l)))],'Color','b'),hold on
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %grid on
    
    % horizontal smoothness between segments
    % note the sid should be aligned along strike direction
    
    if sseg ~= 0 && ~isempty(sid)
        for k = 1:1:length(sid(:,1))
            for j = 1:1:layer
                i1 = find(ip==sid(k,1)&id==j,1,'last');
                i2 = find(ip==sid(k,2)&id==j,1,'first');
                if ~isempty(i1) && ~isempty(i2)
                    col = zeros(1,N*sum(fault_type));
                    col(1,(in(i1)-1)*sum(fault_type)+1) = -2/sqrt((xs(i1)-xs(i2))^2+(ys(i1)-ys(i2))^2);
                    col(1,(in(i2)-1)*sum(fault_type)+1) = 2/sqrt((xs(i1)-xs(i2))^2+(ys(i1)-ys(i2))^2);
                    smooth = [smooth;col];
                    col = zeros(1,N*sum(fault_type));
                    col(1,(in(i1)-1)*sum(fault_type)+2) = -2*fdip/sqrt((xs(i1)-xs(i2))^2+(ys(i1)-ys(i2))^2);
                    col(1,(in(i2)-1)*sum(fault_type)+2) = 2*fdip/sqrt((xs(i1)-xs(i2))^2+(ys(i1)-ys(i2))^2);
                    smooth = [smooth;col];
                    if plt == 1
                        line([xs(in(i2));xs(in(i1))],[ys(in(i2));ys(in(i1))],[zs(in(i2));zs(in(i1))],'Color','r'),hold on
                    end
                end
            end
        end
        
        % for ssid, the third row means its intersection, 1 means the start
        % of the head of id1 intersects with the segment id2, 2 means the
        % end of id1 intersects with the segment id2.
        if ~isempty(ssid)
            for k = 1:1:length(ssid(:,1))
                for j = 1:1:layer
                    if ssid(k,3) == 2
                        i1 = find(ip==ssid(k,1)&id==j,1,'last');
                    elseif ssid(k,3) == 1
                        i1 = find(ip==ssid(k,1)&id==j,1,'first');
                    else
                        continue
                    end
                    i =  find(ip==ssid(k,2)&id==j);
                    if ~isempty(i1) && ~isempty(i)
                        if ssid(k,3) == 2
                            [~, i2] = min((xs(i)-xs(in(i1))-len(in(i1))*sin(strike(in(i1)))).^2+(ys(i)-ys(in(i1))-len(in(i1))*cos(strike(in(i1)))).^2);
                        elseif ssid(k,3) == 1
                            [~, i2] = min((xs(i)-xs(in(i1))).^2+(ys(i)-ys(in(i1))).^2);
                        end
                        i2 = i(i2);
                        col = zeros(1,N*sum(fault_type));
                        col(1,(in(i1)-1)*sum(fault_type)+1) = -2/sqrt((xs(i1)-xs(i2))^2+(ys(i1)-ys(i2))^2);
                        col(1,(in(i2)-1)*sum(fault_type)+1) = 2/sqrt((xs(i1)-xs(i2))^2+(ys(i1)-ys(i2))^2);
                        smooth = [smooth;col];
                        if plt == 1
                            line([xs(in(i2));xs(in(i1))],[ys(in(i2));ys(in(i1))],[zs(in(i2));zs(in(i1))],'Color','g'),hold on
                        end
                    end
                end
            end
        end
        
    end
    
    ed_weigh = 10;
    % 0 edge constrains
    if top ~= 0
        ii = find(id == 1);
        for j = 1:1:length(ii)
            col = zeros(1,N*sum(fault_type));
            col((ii(j)-1)*sum(fault_type)+1) = 2/wid(ii(j))*ed_weigh;
            smooth = [smooth;col];
            col = zeros(1,N*sum(fault_type));
            col(ii(j)*sum(fault_type)) = 2/wid(ii(j))*ed_weigh;
            smooth = [smooth;col];
            if plt == 1
                plot3(xs(ii(j)),ys(ii(j)),zs(ii(j)),'m*'),hold on
            end
        end
    end    
        
    
    if bot ~= 0
        ii = find(id == layer);
        for j = 1:1:length(ii)
            col = zeros(1,N*sum(fault_type));
            col((ii(j)-1)*sum(fault_type)+1) = 2/wid(ii(j))*ed_weigh;
            smooth = [smooth;col];
            col = zeros(1,N*sum(fault_type));
            col(ii(j)*sum(fault_type)) = 2/wid(ii(j))*ed_weigh;
            smooth = [smooth;col];
            if plt == 1
                plot3(xs(ii(j)),ys(ii(j)),zs(ii(j)),'m*'), hold on
            end
        end
    end
    
    if side ~=0
        for j = 1:1:layer
            for k = 1:1:length(sideid(:,1))
                if sideid(k,2) ~=0
                    ii = find(ip==sideid(k,1)&id==j,1,'first');
                end
                if sideid(k,3) ~=0
                    ii = [ii,find(ip==sideid(k,1)&id==j,1,'last')];
                else
                    continue
                end
                for jj = 1:1:length(ii)
                    col = zeros(1,N*sum(fault_type));
                    col((ii(jj)-1)*sum(fault_type)+1) = 2/len(ii(jj))*ed_weigh;
                    smooth = [smooth;col];
                    col = zeros(1,N*sum(fault_type));
                    col(ii(jj)*sum(fault_type)) = 2/len(ii(jj))*ed_weigh;
                    smooth = [smooth;col];
                end
                if plt == 1
                    plot3(xs(ii),ys(ii),zs(ii),'m*'), hold on
                end
            end
        end
    end          
    
