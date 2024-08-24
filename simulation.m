%%
clear;
clc;

%% Specification of the locations of material points   将物质点位置存入数组
length = 100.0;
% length: Total length of the plate  长度
width = 50.0;
% Total width of the plate  宽度
ndivx = 200;
% ndivx: Number of divisions in x direction - except boundary region
% 非边界点的离散数量--x轴
ndivy = 100;
% ndivx: Number of divisions in y direction - except boundary region
%非边界点的离散数量--y轴
dx = length / ndivx;
% dx: Spacing between material points  物质点间距
thick = dx;
% Total thickness of the plate  厚度
nbnd_x = 6;
nbnd_y = 3;
% nbnd: Number of divisions in the boundary region  
%边界上的点数
totnode = (nbnd_x+ndivx)*ndivy + nbnd_y * ndivx;
% totnode: Total number of material points  
%总物质点数量

coord = zeros(totnode, 2);
% coord: Material point locations  
%物质点位置数组

matnum = zeros(totnode,1);
% matnum: Material number of a material point  
%物质点所属材料分类数组 123456分别表示地层123456

nnum = 0;
% nnum: Material point number  循环用到的当前物质点索引

crlength = 8.0;
crack_angle = 30;  % 裂纹倾斜角度，单位为度
crack_angle_radians = deg2rad(crack_angle);  % 将角度转换为弧度
cos_crack_angle = cos(crack_angle_radians);
sin_crack_angle = sin(crack_angle_radians);


% Material points of the bar
for i = 1:ndivx
    for j = 1:ndivy
        %将物质点位置数据存入coord数组
        coordx = (-1.0*length/2.0) + (dx/2.0) + (i - 1)*dx;
        coordy = (-1.0*width/2.0) + (dx/2.0) + (j - 1)*dx;
        nnum = nnum + 1;
        
        coord(nnum, 1) = coordx;
        coord(nnum, 2) = coordy;
        x_transformed = coordx * cos_crack_angle + coordy * sin_crack_angle;
         if x_transformed <= 0  % 裂纹左侧
            if(j <= 20)
                matnum(nnum ,1) = 6;
            elseif(j <= 32)
                matnum(nnum ,1) = 5;
            elseif(j <= 52)
                matnum(nnum ,1) = 4;
            elseif(j <= 76)
                matnum(nnum ,1) = 3;
            elseif(j <= 88)
                matnum(nnum ,1) = 2;
            elseif(j <= 100)
                matnum(nnum ,1) = 1;
            end
         elseif x_transformed > 0  % 裂纹右侧
             if(j <= 8)
                 matnum(nnum ,1) = 7 ;
             elseif(j <= 24)
                 matnum(nnum ,1) = 5;
             elseif(j <= 44)
                 matnum(nnum ,1) = 4;
             elseif(j <= 72)
                 matnum(nnum ,1) = 3;
             elseif(j <= 88)
                 matnum(nnum ,1) = 2;
             elseif(j <= 100)
                 matnum(nnum ,1) = 1;
             end
         end
        %判断位置，区分不同材料，并用matnum数组记录
    end
end

totint = nnum;
for i = 1:nbnd_x/2
    for j = 1: ndivy
        nnum = nnum + 1;
        coord(nnum, 1) = 1.0/2.0*length + (dx/2.0) + (i - 1)*dx;
        coord(nnum, 2) = -1.0/2.0*width + (dx/2.0) + (j - 1)*dx;
        if(j <= 8)
            matnum(nnum ,1) = 7;
        elseif(j <= 20)
            matnum(nnum ,1) = 5;
        elseif(j <= 44)
            matnum(nnum ,1) = 4;
        elseif(j <= 68)
            matnum(nnum ,1) = 3;
        elseif(j <= 88)
            matnum(nnum ,1) = 2;
        elseif(j <= 100)
            matnum(nnum ,1) = 1;
        end
    end
end

totright = nnum;

% Material points of the boundary region - left
for i = 1: nbnd_x/2
    for j = 1: ndivy
        nnum = nnum + 1;
        coord(nnum, 1) = -1.0/2.0*length - (dx/2.0) - (i - 1)*dx;
        coord(nnum, 2) = -1.0/2.0*width + (dx/2.0) + (j - 1)*dx;
        if(j <= 24)
            matnum(nnum ,1) = 6;
        elseif(j <= 32)
            matnum(nnum ,1) = 5;
        elseif(j <= 44)
            matnum(nnum ,1) = 4;
        elseif(j <= 72)
            matnum(nnum ,1) = 3;
        elseif(j <= 88)
            matnum(nnum ,1) = 2;
        elseif(j <= 100)
            matnum(nnum ,1) = 1;
        end
    end
end
totleft = nnum;

for i = 1: nbnd_y
    for j = 1: ndivx
        
        coordx = (-1.0*length/2.0) + (dx/2.0) + (j - 1)*dx;
        coordy = (-1.0*width/2.0) - (dx/2.0) - (i - 1)*dx;
        nnum = nnum + 1;
        
        coord(nnum, 1) = coordx;
        coord(nnum, 2) = coordy;
        
        x_transformed = coordx * cos_crack_angle + coordy * sin_crack_angle;
        if(x_transformed <= 0)
            matnum(nnum ,1) = 6;
        elseif(x_transformed > 0)
            matnum(nnum ,1) = 7;
        end
    end
end

totbottom =nnum;
totnode = nnum;

%% Initialization of fail flag array  初始化失效数组
maxfam = 100;
% maxfam: Maximum number of material points inside a horizon of a material
% point  物质点近场范围最大物质点数量

fail = ones(totnode, maxfam);

%% Determination of material points inside the horizon of each material point
delta = 3.015 * dx;
% delta: Horizon  近场范围

pointfam = int32(zeros(totnode, 1));
% pointfam: index array to find the family members in nodefam array  
% 索引数组，用来查询nodefam中的成员
numfam = int32(zeros(totnode, 1));
% numfam: Number of family members of each material point
% 每个物质点的近场范围所包含的物质点数量
nodefam = int32(zeros(1000000, 1));
% nodefam: array containing family members of all material points
% 包含了所有物质点的族成员的数组

for i = 1:totnode
    if (i == 1)
        pointfam(i, 1) = 1;
    else
        pointfam(i, 1) = pointfam(i - 1, 1) + numfam(i - 1, 1);
    end
    for j = 1:totnode
        idist = sqrt((coord(j, 1) - coord(i, 1))^2+(coord(j, 2) - coord(i, 2))^2);
        if (i ~= j)
            if (idist <= delta)
                numfam(i, 1) = numfam(i, 1) + 1;
                nodefam(pointfam(i, 1)+numfam(i, 1)-1, 1) = j;
            end
        end
    end
end

%% Definition of the crack surface
% crlength = 44.0;
% % crlength: Crack length
% 
% % PD bonds penetrating through the crack surface are broken
% for i = 1:totnode
%     for j = 1:numfam(i, 1)
%         cnode = nodefam(pointfam(i, 1)+j-1, 1);
%         if ((coord(cnode, 1) > 0.0) && (coord(i, 1) < 0.0))
%             if ((abs(coord(i, 2)) - ((crlength / 2.0)-3)) <= 1.0e-10)
%                 fail(i, j) = 0;
%             elseif ((abs(coord(cnode, 2)) - ((crlength / 2.0)-3)) <= 1.0e-10)
%                 fail(i, j) = 0;
%             end
%         elseif ((coord(i, 1) > 0.0) && (coord(cnode, 1) < 0.0))
%             if ((abs(coord(i, 2)) - ((crlength / 2.0)-3)) <= 1.0e-10)
%                 fail(i, j) = 0;
%             elseif ((abs(coord(cnode, 2)) - ((crlength / 2.0)-3)) <= 1.0e-10)
%                 fail(i, j) = 0;
%             end
%         end
%     end
% end




% PD bonds penetrating through the crack surface are broken
for i = 1:totnode
        for j = 1:numfam(i, 1)
            cnode = nodefam(pointfam(i, 1) + j - 1, 1);

            % 对节点坐标进行坐标变换，将裂纹旋转到垂直坐标轴
            x_transformed_i = coord(i, 1) * cos_crack_angle + coord(i, 2) * sin_crack_angle;
            x_transformed_cnode = coord(cnode, 1) * cos_crack_angle + coord(cnode, 2) * sin_crack_angle;

            % 计算裂纹位置的新x坐标
            crack_center_x = 0;  % 这里假设裂纹中心位于坐标原点
            new_x_i = x_transformed_i - crack_center_x;
            new_x_cnode = x_transformed_cnode - crack_center_x;

            % 添加条件，仅在顶部6米以下的节点上执行断裂检查
            if coord(i, 2) <= 1.0  % 顶部6米以下
                if ((new_x_i > 0.0) && (new_x_cnode < 0.0))
                    if ((abs(new_x_i) - ((crlength / 2.0) - 3)) <= 1.0e-10)
                        fail(i, j) = 0;
                    elseif ((abs(new_x_cnode) - ((crlength / 2.0) - 3)) <= 1.0e-10)
                        fail(i, j) = 0;
                    end
                elseif ((new_x_i < 0.0) && (new_x_cnode > 0.0))
                    if ((abs(new_x_i) - ((crlength / 2.0) - 3)) <= 1.0e-10)
                        fail(i, j) = 0;
                    elseif ((abs(new_x_cnode) - ((crlength / 2.0) - 3)) <= 1.0e-10)
                        fail(i, j) = 0;
                    end
                end
            end
        end
end

%% Determination of surface correction factors  表面修正因子
radij = dx / 2.0d0;
%radij: Material point radius  物质点半径
area = dx * dx;
% area: Cross-sectional area  物质点横截面面积
vol = area * dx;
% vol: Volume of a material point  物质点体积

dens_1 = 1700.0;
dens_2 = 1750.0;
dens_3 = 1800.0;
dens_4 = 1850.0;
dens_5 = 1950.0;
dens_6 = 2500.0;
% dens: Density  各层密度
emod_1 = 4.0e6;
emod_2 = 8.5e6;
emod_3 = 9.0e6;
emod_4 = 1.05e7;
emod_5 = 1.2e7;
emod_6 = 10.0e9;
% emod: Elastic modulus  弹性模量
bc_1 = 9.0 * emod_1 / (pi * thick * (delta^3));
bc_2 = 9.0 * emod_2 / (pi * thick * (delta^3));
bc_3 = 9.0 * emod_3 / (pi * thick * (delta^3));
bc_4 = 9.0 * emod_4 / (pi * thick * (delta^3));
bc_5 = 9.0 * emod_5 / (pi * thick * (delta^3));
bc_6 = 9.0 * emod_6 / (pi * thick * (delta^3));
% bc: Bond constant  键常数

disp = zeros(totnode, 2);
% disp: displacement of a material point  点位移数组
stendens = zeros(totnode, 2);
% stendens: strain energy of a material point  物质点应变能数组
fncst = ones(totnode, 2);
% fncst: surface correction factors of a material point, 
% 1:loading 1,2:loading 2  
% 物质点的表面修正因子

sedload1 = 0;

for i = 1:totnode
    disp(i, 1) = 0.001 * coord(i, 1);
    disp(i, 2) = 0.0;
end

for i = 1:totnode
    stendens(i, 1) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2);

        if (idist <= delta - radij)
            fac = 1.0;
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij);
        else
            fac = 0.0;
        end
        
        if (matnum(i,1) == 1)
            bc = bc_1;
        elseif (matnum(i,1) == 2)
            bc = bc_2;  
        elseif (matnum(i,1) == 3)
            bc = bc_3;  
        elseif (matnum(i,1) == 4)
            bc = bc_4;  
        elseif (matnum(i,1) == 5)
            bc = bc_5;  
        elseif (matnum(i,1)==6 || matnum(i,1)==7)
            bc = bc_6;  
        end
        
        stendens(i, 1) = stendens(i, 1) + 0.5 * 0.5 * bc * ((nlength - idist) / idist)^2 * idist * vol * fac;
    end
    % Calculation of surface correction factor in x direction
    % by finding the ratio of the analytical strain energy density value
    % to the strain energy density value obtained from PD Theory
    if (matnum(i,1)==1)
        sedload1 = 9.0 / 16.0 * emod_1 * 1.0e-6;
    elseif(matnum(i,1)==2)
        sedload1 = 9.0 / 16.0 * emod_2 * 1.0e-6;
    elseif(matnum(i,1)==3)
        sedload1 = 9.0 / 16.0 * emod_3 * 1.0e-6;
    elseif(matnum(i,1)==4)
        sedload1 = 9.0 / 16.0 * emod_4 * 1.0e-6;
    elseif(matnum(i,1)==5)
        sedload1 = 9.0 / 16.0 * emod_5 * 1.0e-6;
    elseif(matnum(i,1)==6 || matnum(i,1)==7)
        sedload1 = 9.0 / 16.0 * emod_6 * 1.0e-6;
    end
    
    fncst(i, 1) = sedload1 / stendens(i, 1);
end


sedload2 = 0;

for i = 1:totnode
    disp(i, 1) = 0.0;
    disp(i, 2) = 0.001 * coord(i,2);
end

for i = 1:totnode
    stendens(i, 2) = 0.0;
    for j = 1:numfam(i, 1)
        cnode = nodefam(pointfam(i, 1)+j-1, 1);
        idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2);
        nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2);

        if (idist <= delta - radij)
            fac = 1.0;
        elseif (idist <= delta + radij)
            fac = (delta + radij - idist) / (2.0 * radij);
        else
            fac = 0.0;
        end
        
        if (matnum(i,1) == 1)
            bc = bc_1;
        elseif (matnum(i,1) == 2)
            bc = bc_2;  
        elseif (matnum(i,1) == 3)
            bc = bc_3;  
        elseif (matnum(i,1) == 4)
            bc = bc_4;  
        elseif (matnum(i,1) == 5)
            bc = bc_5;  
        elseif (matnum(i,1)==6 || matnum(i,1)==7)
            bc = bc_6;  
        end
        stendens(i, 2) = stendens(i,2) + 0.5 * 0.5 * bc * ((nlength - idist) / idist)^2 * idist * vol * fac;
    end
    % Calculation of surface correction factor in x direction
    % by finding the ratio of the analytical strain energy density value
    % to the strain energy density value obtained from PD Theory
    if (matnum(i,1)==1)
        sedload2 = 9.0 / 16.0 * emod_1 * 1.0e-6;
    elseif(matnum(i,1)==2)
        sedload2 = 9.0 / 16.0 * emod_2 * 1.0e-6;
    elseif(matnum(i,1)==3)
        sedload2 = 9.0 / 16.0 * emod_3 * 1.0e-6;
    elseif(matnum(i,1)==4)
        sedload2 = 9.0 / 16.0 * emod_4 * 1.0e-6;
    elseif(matnum(i,1)==5)
        sedload2 = 9.0 / 16.0 * emod_5 * 1.0e-6;
    elseif(matnum(i,1)==6 || matnum(i,1)==7)
        sedload2 = 9.0 / 16.0 * emod_6 * 1.0e-6;
    end
    
    fncst(i, 2) = sedload2 / stendens(i, 2);
end


%% Initialization of displacements and velocities
% 初始化位移和速度数组
vel = zeros(totnode, 2);
disp = zeros(totnode, 2);


%% Stable mass vector computation
% 稳定质量矢量计算
dt = 1.0;
% dt: Time interval  时间步长

massvec = zeros(totnode, 2);
% massvec: massvector for adaptive dynamic relaxation  用于自适应动态松弛的质量向量

for i = 1:totnode
    if(matnum(i,1)==1)
        massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_1 / dx ;%  * 5.0;
        massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_1 / dx ;%  * 5.0;
    elseif(matnum(i,1)==2)
        massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_2 / dx ;%  * 5.0;
        massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_2 / dx ;%  * 5.0;
    elseif(matnum(i,1)==3)
        massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_3 / dx ;%  * 5.0;
        massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_3 / dx ;%  * 5.0;
    elseif(matnum(i,1)==4)
        massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_4 / dx ;%  * 5.0;
        massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_4 / dx ;%  * 5.0;
    elseif(matnum(i,1)==5)
        massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_5 / dx ;%  * 5.0;
        massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_5 / dx ;%  * 5.0;
    elseif(matnum(i,1)==6 || matnum(i,1)==7)
        massvec(i, 1) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_6 / dx ;%  * 5.0;
        massvec(i, 2) = 0.25 * dt * dt * (pi * (delta)^2 * thick) * bc_6 / dx ;%  * 5.0;
    end
end

%% Applied loading
appres = 1.3e4;
% appres: Applied pressure  施加载荷

bforce = zeros(totnode, 2);
% bforce: body load acting on a material point  体力数组

for i =1: totint
     %if(matnum(i,1) == 7||matnum(i,1) == 6||matnum(i,1) == 5||matnum(i,1) == 4)
        bforce(i,1) = 0;
        bforce(i,2) =  -1.0 * appres; 
     %end     
end

% for i =1: totint
%     if matnum(i,1) == 1 || matnum(i,1)== 3 ||matnum(i,1) == 5 ||matnum(i,1) ==7 || matnum(i,1) ==2
%         bforce(i,1) = 0;
%         bforce(i,2) =  -1.0 * appres; 
%     else
%         bforce(i,1) = 0;
%         bforce(i,2) = 0; 
%     end
% end

% index = 0;
% for i =1 :ndivx
%     for j =1:ndivy
%         index = index +1;
%         if(matnum(index,1)==6)
%             bforce(index,1) = 0;
%             bforce(index,2) = 0;
%         else
%             bforce(index,1) = 0;
%            % bforce(index,2) =  -1.0 * appres / (dx); 
%             bforce(index,2) =  -1.0 * appres; 
%         end
%     end
% end    


%% Time integration  时间积分
nt =330;
% nt: Total number of time step  总时间步
scr0 = 0.04472;

dmg = zeros(totint, 1);
maxdmg = 0.5;
maxxc = crlength / 2.0;

pforce = zeros(totnode, 2);
% pforce: total peridynamic force acting on a material point  物质点力密度数组
pforceold = zeros(totnode, 2);
% pforce: total peridynamic force acting on a material point in the
% previous time step 1:x-coord, 2:y-coord  物质点前一时间步的物质点力密度数组
acc = zeros(totnode, 2);
% acc: acceleration of a material point  加速度

velhalf = zeros(totnode, 2);
velhalfold = zeros(totnode, 2);
% vel: velocity of a material point  物质点速度

%coord_disp_pd_nt_675 = zeros(totnode, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
%coord_disp_pd_nt_750 = zeros(totnode, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
%coord_disp_pd_nt_825 = zeros(totnode, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
%coord_disp_pd_nt_1000 = zeros(totnode, 5);
% Peridynamic displacement and Analytical displacement of all points at time step of nt

coord_disp_pd_nt = zeros(totnode, 4);
% Peridynamic displacement and Analytical displacement of all points at time step of nt
coord_disp_pd_12 = zeros(totnode, 5);

horiCnt = 0;
horizontal_disps = zeros(ndivx, 4);
% Peridynamic displacement and Analytical displacement of points at y = 0;
vertiCnt = 0;
vertical_disps = zeros(ndivy, 4);


cn = 0.0;
cn1 = 0.0;
cn2 = 0.0;
for tt = 1:nt
    fprintf("%d/%d\n", tt, nt);
    ctime = tt * dt;
    %if(tt == 1)
    for i = (totint + 1): totright
        disp(i, 1) = 0*tt*dt;     
    end
    
    for i = (totright + 1): totleft
        disp(i, 1) = 0*tt*dt;
    end

    for i = (totleft + 1): totnode
        if(matnum(i,1) == 6)
            vel(i, 2) = 0;
            disp(i, 2) = 0*tt*dt;
        end
    end

    for i = 1:totnode
        dmgpar1 = 0.0;
        dmgpar2 = 0.0;
        pforce(i, 1) = 0.0;
        pforce(i, 2) = 0.0;
        for j = 1:numfam(i, 1)
            cnode = nodefam(pointfam(i, 1)+j-1, 1);
            idist = sqrt((coord(cnode, 1) - coord(i, 1))^2+(coord(cnode, 2) - coord(i, 2))^2);
            nlength = sqrt((coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))^2+(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))^2);
 
            % Volume correction
            if (idist <= delta - radij)
                fac = 1.0;
            elseif (idist <= delta + radij)
                fac = (delta + radij - idist) / (2.0 * radij);
            else
                fac = 0.0;
            end
            
            if (abs(coord(cnode, 2)-coord(i, 2)) <= 1.0e-10)
                theta = 0.0;
            elseif (abs(coord(cnode, 1)-coord(i, 1)) <= 1.0e-10)
                theta = 90.0 * pi / 180.0; 
            else
                theta = atan(abs(coord(cnode, 2)-coord(i, 2))/abs(coord(cnode, 1)-coord(i, 1)));
            end
            % Determination of the surface correction between two material points
            scx = (fncst(i, 1) + fncst(cnode, 1)) / 2.0;
            scy = (fncst(i, 2) + fncst(cnode, 2)) / 2.0;
            scr = 1.0 / (((cos(theta))^2.0 / (scx)^2.0) + ((sin(theta))^2.0 / (scy)^2.0));
            scr = sqrt(scr);
            
            if (matnum(i,1) == 1)
                bc = bc_1;
            elseif(matnum(i,1)==2)
                bc = bc_2;
            elseif(matnum(i,1)==3)
                bc = bc_3;
            elseif(matnum(i,1)==4)
                bc = bc_4;
            elseif(matnum(i,1)==5)
                bc = bc_5;
            elseif(matnum(i,1)==6 || matnum(i,1)==7)
                bc = bc_6;
            end
             x_transformed = coord(i,1) * cos_crack_angle + coord(i,2) * sin_crack_angle;
            % Calculation of the peridynamic force in x direction
            % acting on a material point i due to a material point j
            if (fail(i,j)==1)
                dforce1 = bc*(nlength - idist)/idist*vol*scr*fac*(coord(cnode, 1) + disp(cnode, 1) - coord(i, 1) - disp(i, 1))/nlength;
                dforce2 = bc*(nlength - idist)/idist*vol*scr*fac*(coord(cnode, 2) + disp(cnode, 2) - coord(i, 2) - disp(i, 2))/nlength;
            else
                dforce1 = 0;
                dforce2 = 0;
            end
            pforce(i, 1) = pforce(i, 1) + dforce1;
            pforce(i, 2) = pforce(i, 2) + dforce2;
            
            
            % 定义失效区
            if ((nlength - idist)/idist > scr0)
                %if (abs(coord(i, 1)) <= (length/4.0))
                %if (x_transformed>0.0)
                    fail(i, j) = 0;
                %end
                %end
            end
            
            dmgpar1 = dmgpar1 + fail(i, j)*vol*fac;
            dmgpar2 = dmgpar2 + vol*fac;
        end
        % Calculation of the damage parameter  计算损伤系数
        dmg(i, 1) = 1.0 - dmgpar1/dmgpar2;
         if ((dmg(i, 1) > maxdmg) && (coord(i, 2) > maxxc))
             maxxc = coord(i, 2);
         end
    end

    newcrlength = maxxc - crlength / 2.0;
    
    % Adaptive dynamic relaxation  自适应松弛

    for i = 1:totnode
        if (velhalfold(i, 1) ~= 0.0)
            cn1 = cn1 - disp(i, 1) * disp(i, 1) * (pforce(i, 1) / massvec(i, 1) - pforceold(i, 1) / massvec(i, 1)) / (dt * velhalfold(i, 1));
        end
        if (velhalfold(i, 2) ~= 0.0)
            cn1 = cn1 - disp(i, 2) * disp(i, 2) * (pforce(i, 2) / massvec(i, 2) - pforceold(i, 2) / massvec(i, 2)) / (dt * velhalfold(i, 2));
        end
        cn2 = cn2 + disp(i, 1) * disp(i, 1);
        cn2 = cn2 + disp(i, 2) * disp(i, 2);
    end

    if (cn2 ~= 0.0)
        if ((cn1 / cn2) > 0.0)
            cn = 2.0 * sqrt(cn1/cn2);
        else
            cn = 0.0;
        end
    else
        cn = 0.0;
    end

    if (cn > 2.0)
        cn = 1.9;
    end

    for i = 1:totint
        % Integrate acceleration over time.
            if (tt == 1)
                velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
                velhalf(i, 2) = 1.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2)) / 2.0;
            else
                velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
                velhalf(i, 2) = ((2.0 - cn * dt) * velhalfold(i, 2) + 2.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2))) / (2.0 + cn * dt);
            end
            
            vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
            vel(i, 2) = 0.5 * (velhalfold(i, 2) + velhalf(i, 2));
            disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;
            disp(i, 2) = disp(i, 2) + velhalf(i, 2) * dt;
            
            velhalfold(i, 1) = velhalf(i, 1);
            velhalfold(i, 2) = velhalf(i, 2);
            pforceold(i, 1) = pforce(i, 1);
            pforceold(i, 2) = pforce(i, 2); 
    end
    
    for i = (totint + 1):totright
        % Integrate acceleration over time.
%         if (matnum(i,1) == 6)
%             vel(i, 2) = 0;
%             disp(i, 2) = 0;
%         else
            if (tt == 1)
                velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
                velhalf(i, 2) = 1.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2)) / 2.0;
            else
                velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
                velhalf(i, 2) = ((2.0 - cn * dt) * velhalfold(i, 2) + 2.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2))) / (2.0 + cn * dt);
            end
            
            vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
            vel(i, 2) = 0.5 * (velhalfold(i, 2) + velhalf(i, 2));
            disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;
            disp(i, 2) = disp(i, 2) + velhalf(i, 2) * dt;
            
            velhalfold(i, 1) = velhalf(i, 1);
            velhalfold(i, 2) = velhalf(i, 2);
            pforceold(i, 1) = pforce(i, 1);
            pforceold(i, 2) = pforce(i, 2);
%         end   
    end
    
    for i = (totright + 1):totleft
        % Integrate acceleration over time.
%         if (matnum(i,1) == 6)
%             vel(i, 2) = 0;
%             disp(i, 2) = 0;
%         else
            if (tt == 1)
                velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
                velhalf(i, 2) = 1.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2)) / 2.0;
            else
                velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
                velhalf(i, 2) = ((2.0 - cn * dt) * velhalfold(i, 2) + 2.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2))) / (2.0 + cn * dt);
            end
            
            vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
            vel(i, 2) = 0.5 * (velhalfold(i, 2) + velhalf(i, 2));
            disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;
            disp(i, 2) = disp(i, 2) + velhalf(i, 2) * dt;
            
            velhalfold(i, 1) = velhalf(i, 1);
            velhalfold(i, 2) = velhalf(i, 2);
            pforceold(i, 1) = pforce(i, 1);
            pforceold(i, 2) = pforce(i, 2);
%         end
    end
    
    for i = (totleft + 1):totnode
        % Integrate acceleration over time.
        if (matnum(i,1) ==6)
            vel(i, 2) = 0;
            disp(i, 2) = 0;
        else
            if (tt == 1)
                velhalf(i, 1) = 1.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1)) / 2.0;
                velhalf(i, 2) = 1.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2)) / 2.0;
            else
                velhalf(i, 1) = ((2.0 - cn * dt) * velhalfold(i, 1) + 2.0 * dt / massvec(i, 1) * (pforce(i, 1) + bforce(i, 1))) / (2.0 + cn * dt);
                velhalf(i, 2) = ((2.0 - cn * dt) * velhalfold(i, 2) + 2.0 * dt / massvec(i, 2) * (pforce(i, 2) + bforce(i, 2))) / (2.0 + cn * dt);
            end
            
            vel(i, 1) = 0.5 * (velhalfold(i, 1) + velhalf(i, 1));
            vel(i, 2) = 0.5 * (velhalfold(i, 2) + velhalf(i, 2));
            disp(i, 1) = disp(i, 1) + velhalf(i, 1) * dt;
            disp(i, 2) = disp(i, 2) + velhalf(i, 2) * dt;
            
            velhalfold(i, 1) = velhalf(i, 1);
            velhalfold(i, 2) = velhalf(i, 2);
            pforceold(i, 1) = pforce(i, 1);
            pforceold(i, 2) = pforce(i, 2);
        end
    end
    
%     if (tt == nt)
%         for i = 1:totnode
%             coord_disp_pd_nt_1000(i, 1:5) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), dmg(i, 1)];
%         end
%     elseif (tt == 675)
%         for i = 1:totnode
%              coord_disp_pd_nt_675(i, 1:5) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), dmg(i, 1)];
%         end       
%     elseif (tt == 750)
%         for i = 1:totnode
%              coord_disp_pd_nt_750(i, 1:5) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), dmg(i, 1)];
%         end       
%     elseif (tt == 825)
%         for i = 1:totnode
%             coord_disp_pd_nt_825(i, 1:5) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), dmg(i, 1)];
%         end    
%     end
    
        if (tt == nt)
        for i = 1:totint
            coord_disp_pd_nt(i, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
            coord_disp_pd_12(i, 1:5) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2), dmg(i, 1)];
            if (abs(coord(i, 2)-(dx / 2.0)) <= 1.0e-8)
                horiCnt = horiCnt + 1;
                horizontal_disps(horiCnt, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
            end
            if (abs(coord(i, 1)-(dx / 2.0)) <= 1.0e-8)
                vertiCnt = vertiCnt + 1;
                vertical_disps(vertiCnt, 1:4) = [coord(i, 1), coord(i, 2), disp(i, 1), disp(i, 2)];
            end
        end
        end
end

%% plot
% colormap jet;
% subplot(221)

% scatter(coord_disp_pd_nt_675(:, 1)+scale*coord_disp_pd_nt_675(:, 3), coord_disp_pd_nt_675(:, 2)+scale*coord_disp_pd_nt_675(:, 4), [], coord_disp_pd_nt_675(:, 5), "filled")
% subplot(222)
% scale = 1;
% scatter(coord_disp_pd_nt_750(:, 1)+scale*coord_disp_pd_nt_750(:, 3), coord_disp_pd_nt_750(:, 2)+scale*coord_disp_pd_nt_750(:, 4), [], coord_disp_pd_nt_750(:, 5), "filled")
% subplot(223)
% scale = 1;
% scatter(coord_disp_pd_nt_825(:, 1)+scale*coord_disp_pd_nt_825(:, 3), coord_disp_pd_nt_825(:, 2)+scale*coord_disp_pd_nt_825(:, 4), [], coord_disp_pd_nt_825(:, 5), "filled")
% subplot(224)
% scale = 1;
% scatter(coord_disp_pd_nt_1000(:, 1)+scale*coord_disp_pd_nt_1000(:, 3), coord_disp_pd_nt_1000(:, 2)+scale*coord_disp_pd_nt_1000(:, 4), [], coord_disp_pd_nt_1000(:, 5), "filled")

colormap jet;

% Scatter plot for coord_disp_pd_nt
figure;
scale = 1;
scatter(coord_disp_pd_nt(:, 1)+scale*coord_disp_pd_nt(:, 3), coord_disp_pd_nt(:, 2)+scale*coord_disp_pd_nt(:, 4), [], coord_disp_pd_nt(:, 4), 'filled');
title('nt=330 disp');
axis off;
% xlabel('X');
% ylabel('Y');

% Scatter plot for coord_disp_pd_12
colormap jet
figure;
scale = 1;
scatter(coord_disp_pd_12(:, 1)+scale*coord_disp_pd_12(:, 3), coord_disp_pd_12(:, 2)+scale*coord_disp_pd_12(:, 4), [], coord_disp_pd_12(:, 5), 'filled');
title('nt=330 Fail');
axis off;
% xlabel('X');
% ylabel('Y');

% Plot for horizontal_disps
colormap jet
figure;
plot(horizontal_disps(1:horiCnt, 1), horizontal_disps(1:horiCnt, 4));
title('horizontal_disps');
axis off;
% xlabel('X');
% ylabel('Y');

% Plot for vertical_disps
colormap jet
figure;
plot(vertical_disps(1:vertiCnt, 2), vertical_disps(1:vertiCnt, 4));
title('vertical_disps');
axis off;
% xlabel('X');
% ylabel('Y');

%sqrt(coord_disp_pd_nt(:, 3).^2+coord_disp_pd_nt(:, 4).^2)







