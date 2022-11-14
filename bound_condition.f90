subroutine set_BC_fluent

    !-----------------------------------------------------------------------------
    !
    !   方法：根据fluent网格文件设置边界条件
    !   描述：考虑网格不同边界条件，边界条件由pointweise设置.
    !         根据边界条件设定外加的虚拟单元的指向位置。设置指向目的：便于检索
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    use bc_module
    implicit none
    
    integer i,j,k,l
    integer :: index1 = 1,index2 = 1,index3 = 1,index4 = 1,index0 = 1,index_bd1,index_bd2
    real(prec) :: dis1,dis2,d_D,d_R,d_U,d_L,p_c,p_c1,p_c2,dis_p_c,x1,y1,x2,y2
    
    integer :: i_D,i_R,i_U,i_L,O_sideth
    integer :: indexCell,indexCell2,indexNearCell,indexNearCell2,p1,p2
    integer :: count_temp,indexSide,indexSide2
    
    index1 = 1;index2 = 1;index3 = 1;index4 = 1;index0 = 1;
    n_BCells_D = 0;n_BCells_R = 0;n_BCells_U = 0;n_BCells_L = 0 ;
         
    !-----------------------------------------------------------------------------
    !
    !   描述：如果是自生成网格，先添加边界名称
    !
    !-----------------------------------------------------------------------------
     
    if(grid_set == self)then
        call get_bound_Name
    endif
    
    ! 记录侧边单元的索引
    allocate(BoundCells_index_set(nbdsides))      
    count_temp = 1
    
             
    !-----------------------------------------------------------------------------
    !
    !   描述：给虚拟边界单元编号。四条边循环一遍，拐角处有两条侧边在边界亦可
    !
    !-----------------------------------------------------------------------------
    
    do i = 1,ncells
        do j = 1,4
            indexSide = cellset(i).sides(j)
            if(sideset(indexSide).bc /= Interior)then            
                BoundCells_index_set(count_temp) = i
                cellset(i).nearcells(j) = ncells + count_temp
                cellset(ncells+count_temp).index = ncells + count_temp
                count_temp = count_temp + 1
            end if            
        end do    
    end do
    
    !-----------------------------------------------------------------------------
    !
    !   描述：设置虚拟边界单元的指向，也即index。循环边界指向循环对应位置，其余指向相邻边界单元
    !
    !-----------------------------------------------------------------------------
    
    do i = 1,nbdsides
        indexCell = BoundCells_index_set(i)
        do j = 1,4
            indexNearCell = cellset(indexCell).nearcells(j)
            indexSide = cellset(indexCell).sides(j)
             
            if(sideset(indexSide).bc == Periodic)then
                
                !-----------------------------------------------------------------------------
                !
                !   方法：寻找循环侧对应边。在此使用边的中心点比对  ##上下循环##    
                !
                !-----------------------------------------------------------------------------
             
                ! 本单元边界边中心点  
                x1 = (xy_coor(sideset(indexSide).nodes(1),1)+xy_coor(sideset(indexSide).nodes(2),1))*0.5_prec  
                y1 = (xy_coor(sideset(indexSide).nodes(1),2)+xy_coor(sideset(indexSide).nodes(2),2))*0.5_prec  
                
                ! 对应边中心点
                Loop03:do k = 1,nbdsides
                    indexCell2 = BoundCells_index_set(k)
                    Loop04:do l = 1,4                            
                        indexNearCell2 = cellset(indexCell2).nearcells(l)
                        indexSide2 = cellset(indexCell2).sides(l)

                        ! 对应单元边界侧边中心点    
                        if(sideset(indexSide2).bc == Periodic .and. indexSide2 /= indexSide)then
                            
                            x2 = (xy_coor(sideset(indexSide2).nodes(1),1)+xy_coor(sideset(indexSide2).nodes(2),1))*0.5_prec
                            y2 = (xy_coor(sideset(indexSide2).nodes(1),2)+xy_coor(sideset(indexSide2).nodes(2),2))*0.5_prec
                            
                            dis_p_c = sqrt((x1-x2)**2 + (y1-y2)**2)

                            ! 若中点横坐标重合，则对应。计算域若倾斜一定角度无法计算
                            if((abs(dis_p_c-ylong)<1.0e-10 .and. abs(y1-y2)>1.0e-10) .OR. (abs(dis_p_c-xlong)<1.0e-10 .and. abs(x1-x2)>1.0e-10))then
                                cellset(indexNearCell).index = indexCell2
                                cellset(indexNearCell2).index = indexCell
                                exit Loop03     
                            end if
                        end if
                    end do Loop04
                end do Loop03                         
            elseif(sideset(indexSide).bc /= Periodic .and. indexNearCell > ncells)then       
                cellset(indexNearCell).index = indexCell                   
            end if
        end do
    end do

end subroutine set_BC_fluent
    
subroutine get_bound_Name

    !-----------------------------------------------------------------------------
    !
    !   方法：自生成网格边界条件设置
    !   描述：根据算例进行边界条件设置。自生成网格一般计算矩形计算域
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision 
    use global_var
    use type_module
    use bc_module
    use parameter_setting
    implicit none
    
    integer :: i,j,k,l,vertexth1,vertexth2
    character(len = 16) :: bound

    do i = 1,nsides
        if(sideset(i).nearcells(1) == 0 .or. sideset(i).nearcells(2) == 0)then
            vertexth1 = sideset(i).nodes(1)
            vertexth2 = sideset(i).nodes(2)
            if(abs(xy_coor(vertexth1,2)-yl) < 1.0e-10 .and. abs(xy_coor(vertexth2,2)-yl) < 1.0e-10)then
                bound = 'D'
                call case_bc(case_comp,bound,sideset(i).bc)     
            end if
            if(abs(xy_coor(vertexth1,1)-xr) < 1.0e-10 .and. abs(xy_coor(vertexth2,1)-xr) < 1.0e-10)then
                bound = 'R'
                call case_bc(case_comp,bound,sideset(i).bc) 
                
            end if
            if(abs(xy_coor(vertexth1,2)-yr) < 1.0e-10 .and. abs(xy_coor(vertexth2,2)-yr) < 1.0e-10)then
                bound = 'U'
                call case_bc(case_comp,bound,sideset(i).bc)                 
            end if
            if(abs(xy_coor(vertexth1,1)-xl) < 1.0e-10 .and. abs(xy_coor(vertexth2,1)-xl) < 1.0e-10)then
                bound = 'L'
                call case_bc(case_comp,bound,sideset(i).bc)               
            end if
        else
            sideset(i).bc = Interior
        end if
        
    end do
end subroutine get_bound_Name
    
subroutine update_Boundary_fluent
 
    !-----------------------------------------------------------------------------
    !
    !   方法：更新边界
    !   描述：每一时间步，每一时间层都需要更新边界，不同的算例可能有不同的需求
    !         Wall边界设置比较麻烦，不确定对非结构网格而言，设置是否正确。 
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    use bc_module
    implicit none
    
    integer :: i,j,k,l,m,indexCell,indexNearCell,indexSide,index_temp,p1,p2,index,th1,nearcells(4),sides(4)
    integer :: ikesi,ieta,LC_sideth,RC_sideth,RC_k
    real(prec) :: x,y,ss,tt,ruvp(4),coor(2)
    real(prec) :: x1,y1,x2,y2,xa,ya,vertexth1,vertexth2
    
    if(case_comp == DoubleMach_case)then
        
        ! 双马赫算例
        do i = 1,nbdsides    
            indexCell = BoundCells_index_set(i)
            do j = 1,4     
                indexSide = cellset(indexCell).sides(j)
                indexNearCell = cellset(indexCell).nearcells(j)   
                if(sideset(indexSide).bc == Inlet_vent)then
                    
                    ! 左边界
                    ruvp(1)  = 8.0_prec
                    ruvp(2)  = 7.145_prec
                    ruvp(3)  = -4.125_prec
                    ruvp(4)  = 116.5_prec    
                    do k = 1,nsp
                        do l = 1,nsp
                            cellset(indexNearCell).spvalue_ori(l,k,:) = ruvp
                        end do
                    end do       
                elseif(sideset(indexSide).bc == Outlet_vent)then
                    
                    ! 上边界 和 右边界
                    do k = 1, nsp
                        do l = 1,nsp
                            x = cellset(indexNearCell).sp_coor(k,l,1)
                            y = cellset(indexNearCell).sp_coor(k,l,2)
                            ss = 10.0_prec
                            tt = nt_temp*dt
                            if(y - sqrt(3.0_prec)*(x - 1.0_prec/6.0_prec) > -ss*tt*2.0_prec  )then
                                ruvp(1) = 8.0_prec
                                ruvp(2)  = 7.145_prec
                                ruvp(3)  = -4.125_prec
                                ruvp(4)  = 116.5_prec                         
                            else
                                ruvp(1)  = 1.4_prec
                                ruvp(2)  = 0.0_prec
                                ruvp(3)  = 0.0_prec
                                ruvp(4)  = 1.0_prec
                            endif
                            cellset(indexNearCell).spvalue_ori(k,l,:) = ruvp
                        end do 
                    end do

                elseif(sideset(indexSide).bc == Wall)then
                    
                    ! 下边界
                    vertexth1 = j
                    vertexth2 = j+1
                    if(vertexth1 == 4) vertexth2 = 1
                    x1 = xy_coor(cellset(indexCell).nodes(vertexth1),1)
                    y1 = xy_coor(cellset(indexCell).nodes(vertexth1),2)
                    x2 = xy_coor(cellset(indexCell).nodes(vertexth2),1)
                    y2 = xy_coor(cellset(indexCell).nodes(vertexth2),2)
                    do k = 1,nsp
                        do l = 1,nsp
                            
                            ! (x1,y1),(x2,y2)是边界边的端点。确定方向
                            xa = cellset(indexCell).sp_coor(k,l,1)
                            ya = cellset(indexCell).sp_coor(k,l,2)
                            call sovle_sym_point(x1,y1,x2,y2,xa,ya,indexNearCell,cellset(indexCell).spvalue_ori(k,l,:),cellset(indexNearCell).spvalue_ori(:,:,:))
                        end do
                    end do
 
                    do k = 1,nsp                       
                        do l = 1,nsp                
                            x = cellset(indexNearCell).sp_coor(l,k,1)
                            y = cellset(indexNearCell).sp_coor(l,k,2)            
                            if(x < 1.0_prec/6.0_prec)then
                                ruvp(1)  = 8.0_prec
                                ruvp(2)  = 7.145_prec
                                ruvp(3)  = -4.125_prec
                                ruvp(4)  = 116.5_prec    
                                cellset(indexNearCell).spvalue_ori(l,k,:) = ruvp
                            end if
                        end do        
                    end do                                        
                end if
            end do
        end do

    else     
        ! 其余算例按照边界条件设置
        do i = 1,nbdsides    
            indexCell = BoundCells_index_set(i)
            do j = 1,4     
                indexSide = cellset(indexCell).sides(j)
                indexNearCell = cellset(indexCell).nearcells(j) 
                if(sideset(indexSide).bc == Inlet_vent )then
                    if(case_comp == VortexShock_case)then
                        do k = 1,nsp
                            do l = 1,nsp
                                call VortexShock_init(cellset(indexNearCell).sp_coor(k,l,:),cellset(indexNearCell).spvalue_ori(k,l,:))
                            end do                
                        end do  
                    else 
                        LC_sideth = j
                        do k = 1,4
                        if(cellset(indexNearCell).nearcells(k) == indexCell)then
                               RC_sideth = k 
                            end if
                        end do                         
                        do k = 1,nsp
                            call get_RC_l(LC_sideth,RC_sideth,k,RC_k)                 
                            
                            !边界附近的求解点，确定虚拟单元内的与之连线和边界垂直的点
                            if(j==1)then
                                ikesi = 1
                                ieta = k                           
                            elseif(j==2)then
                                ikesi = k
                                ieta =  nsp   
                            elseif(j==3)then
                                ikesi = nsp
                                ieta = k
                            elseif(j==4)then
                                ikesi = k
                                ieta =  1     
                            end if                                                               
                            call sovle_outflow_point(RC_sideth,RC_k,cellset(indexCell).spvalue_ori(ikesi,ieta,:),cellset(indexNearCell).spvalue_ori(:,:,:))                            
                        end do                    
                    end if
                elseif(sideset(indexSide).bc == Outlet_vent)then 
                    if(case_comp == ShuOsher_case)then 
                        ! nothing
                    else
                        LC_sideth = j
                        do k = 1,4
                        if(cellset(indexNearCell).nearcells(k) == indexCell)then
                               RC_sideth = k 
                            end if
                        end do                         
                        do k = 1,nsp
                            call get_RC_l(LC_sideth,RC_sideth,k,RC_k)               
                            
                            !边界附近的求解点，确定虚拟单元内的与之连线和边界垂直的点
                            if(j==1)then
                                ikesi = 1
                                ieta = k                           
                            elseif(j==2)then
                                ikesi = k
                                ieta =  nsp   
                            elseif(j==3)then
                                ikesi = nsp
                                ieta = k
                            elseif(j==4)then
                                ikesi = k
                                ieta =  1     
                            end if                                                               
                            call sovle_outflow_point(RC_sideth,RC_k,cellset(indexCell).spvalue_ori(ikesi,ieta,:),cellset(indexNearCell).spvalue_ori(:,:,:))                            
                        end do    
                    end if               
                elseif(sideset(indexSide).bc == Periodic)then
                    
                    !虚拟边界单元index与指向相同。周期边界
                    nearcells = cellset(indexNearCell).nearcells
                    sides = cellset(indexNearCell).sides
                    cellset(indexNearCell).spvalue_ori(:,:,:) = cellset(cellset(indexNearCell).index).spvalue_ori(:,:,:)
                    !cellset(indexNearCell) = cellset(cellset(indexNearCell).index)
                    
                    cellset(indexNearCell).nearcells = nearcells
                    cellset(indexNearCell).sides = sides

                elseif(sideset(indexSide).bc == Wall)then
                    vertexth1 = j
                    vertexth2 = j+1
                    if(vertexth1 == 4) vertexth2 = 1
                    x1 = xy_coor(cellset(indexCell).nodes(vertexth1),1)
                    y1 = xy_coor(cellset(indexCell).nodes(vertexth1),2)
                    x2 = xy_coor(cellset(indexCell).nodes(vertexth2),1)
                    y2 = xy_coor(cellset(indexCell).nodes(vertexth2),2)
                    do k = 1,nsp
                        do l = 1,nsp
                            ! (x1,y1),(x2,y2)是边界边的端点。确定方向
                            xa = cellset(indexCell).sp_coor(k,l,1)
                            ya = cellset(indexCell).sp_coor(k,l,2)
                            call sovle_sym_point(x1,y1,x2,y2,xa,ya,indexNearCell,cellset(indexCell).spvalue_ori(k,l,:),cellset(indexNearCell).spvalue_ori(:,:,:))   
                            
                        end do
                    end do
                end if
            end do
        end do
    endif

end subroutine update_Boundary_fluent
    
subroutine case_bc(case_th,bound,bound_condition)  
 
    !-----------------------------------------------------------------------------
    !
    !   方法：根据算例设置自生成网格边界条件
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision 
    use parameter_setting
    use bc_module
    implicit none
    
    integer :: case_th
    character(len = 16) :: bound,bound_condition

    !#  Interior = 'Interior',Wall = '3',Inlet_vent = '4',Outlet_vent = '5',Periodic = '20'
    
    if(case_th == equEntropy_case .or. case_th == SinWave_2D_case)then
        bound_condition = Periodic
    end if
    
    if(case_th == SodShockTube_case .or. case_th == LaxShockTube_case .or. case_th == ShuOsher_case )then
        select case(dire_shock)
        case(direct_x)
            if(bound == 'D' .or. bound == 'U')then
                bound_condition = Periodic
            elseif(bound == 'L')then
                bound_condition = Inlet_vent
            elseif(bound == 'R')then
                bound_condition = Outlet_vent
            end if
        case(direct_y)
             if(bound == 'L' .or. bound == 'R')then
                bound_condition = Periodic
            elseif(bound == 'D')then
                bound_condition = Inlet_vent
            elseif(bound == 'U')then
                bound_condition = Outlet_vent
            end if
        end select  
    end if
    
    if(case_th == steadyShock_case)then
            if(bound == 'D')then
                bound_condition = Inlet_vent
            elseif(bound == 'U')then
                bound_condition = Outlet_vent
            elseif(bound == 'L')then
                bound_condition = Inlet_vent
            elseif(bound == 'R')then
                bound_condition = Outlet_vent
            end if
    end if
    
    if(case_th == Riemann2D_case)then
        if(bound == 'U' .or. bound == 'R')then
            bound_condition = Inlet_vent
        elseif(bound == 'L' .or. bound == 'D')then
            bound_condition = Outlet_vent
        end if    
    end if
    
    if(case_th == DoubleMach_case)then
        if(bound == 'R' .or. bound == 'U')then
            bound_condition = Outlet_vent
        elseif(bound == 'D')then
            bound_condition = Wall
        elseif(bound == 'L')then
            bound_condition = Inlet_vent
        end if    
    end if

    if(case_th == VortexShock_case .or. case_th == CompositeVortexShock_case)then
        if(bound == 'D' .or. bound == 'U')then
            bound_condition = Wall
        elseif(bound == 'L')then
            bound_condition = Inlet_vent
        elseif(bound == 'R')then
            bound_condition = Outlet_vent
        end if    
    end if
    
end subroutine case_bc    
    
subroutine pre_BoundaryCell_fluent
 
    !-----------------------------------------------------------------------------
    !
    !   方法：求解虚拟边界单元信息
    !   描述：设置单元坐标，求解点真实坐标，Jacobi等不随时间变化的值.单元形状与相邻单元沿边界线对称。
    !         对于循环边界条件的case之后可以重新赋值网格相关的变量，不必担心单元形状。先试试再说。
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer :: i,j,k,l,indexCell,indexNearCell,index_temp,p1,p2,index,th1,nearcells(4),pCounter,sideth1,vertexth(4)
    real(prec),dimension(:,:)   :: quad_vertex(4,2)
    real(prec),dimension(:,:,:) :: weight_sps(nsp,nsp,4)
    real(prec) :: coor_C(2),M_Jacobi(4),M_direct(4)
    real(prec) :: x1,y1,x2,y2,xa,ya,x,y,xc1,yc1,xc2,yc2
    integer    :: indexNode,indexSide,temp,vertexth1,vertexth2
    real(prec) :: coor_vertex(4,2)
    
    
         
    !-----------------------------------------------------------------------------
    !
    !   描述：生成虚拟单元。形状关于边界与边界单元对称
    !
    !-----------------------------------------------------------------------------
     
    pCounter = nnodes+1
    Loop1:do i = 1,nbdsides
        
        ! 边界侧边所在单元
        indexCell = BoundCells_index_set(i)
        do j = 1,i-1
            if(BoundCells_index_set(j) == indexCell)cycle Loop1
        end do
        
        do j = 1,4                     
            ! 边界侧边所在单元的邻单元
            indexNearCell = cellset(indexCell).nearcells(j)
            if(indexNearCell > ncells)then    
                
                ! 先录入重合的两点.最后按照之前的规则重新编号
                vertexth(1) = j
                vertexth(2) = j+1
                if(j == 4) vertexth(2) = 1
                cellset(indexNearCell).nodes(1) = cellset(indexCell).nodes(vertexth(1))
                cellset(indexNearCell).nodes(2) = cellset(indexCell).nodes(vertexth(2))
                indexNode = 3
                
                ! 求关于直线的对称点
                do k = 1,4
                    x1 = xy_coor(cellset(indexCell).nodes(vertexth(1)),1)
                    y1 = xy_coor(cellset(indexCell).nodes(vertexth(1)),2)
                    x2 = xy_coor(cellset(indexCell).nodes(vertexth(2)),1)
                    y2 = xy_coor(cellset(indexCell).nodes(vertexth(2)),2)
                    if(k /= vertexth(1) .and. k /= vertexth(2))then                        
                        xa = xy_coor(cellset(indexCell).nodes(k),1)
                        ya = xy_coor(cellset(indexCell).nodes(k),2)
                        
                        call sym_Point(x1,y1,x2,y2,xa,ya,x,y)    ! Symmetry point 
                        
                        xy_coor(pCounter,1) = x
                        xy_coor(pCounter,2) = y                
                        cellset(indexNearCell).nodes(indexNode) = pCounter
                        
                        indexNode = indexNode + 1
                        pCounter = pCounter + 1

                    end if                    
                end do                   
            end if      
        end do 
    end do Loop1
    
    !-----------------------------------------------------------------------------
    !
    !   描述：求解虚拟边界单元的信息。和内部单元信息一致
    !
    !-----------------------------------------------------------------------------
    
    do i = 1,nbdsides
        indexCell = ncells + i
        do j =1,4
            coor_vertex(j,:) = xy_coor(cellset(indexCell).nodes(j),:)
        end do

        ! 对值排序
        call sort_qua(coor_vertex(:,:))

        ! 按值调整编号,逆时针
        do j = 1,4
            do k = 1,4
                if(coor_vertex(j,1)==xy_coor(cellset(indexCell).nodes(k),1).and. coor_vertex(j,2)==xy_coor(cellset(indexCell).nodes(k),2))then
                    temp = cellset(indexCell).nodes(j)   
                    cellset(indexCell).nodes(j) = cellset(indexCell).nodes(k) 
                    cellset(indexCell).nodes(k) = temp
                end if
            end do
        end do 
    end do
    
    do i = 1,nbdsides
        indexCell = BoundCells_index_set(i)
        
        ! 设置虚拟边界单元靠近边界单元的邻单元和边序号
        do j = 1,4
            indexNearCell = cellset(indexCell).nearcells(j)
            indexSide = cellset(indexCell).sides(j)
            if(indexNearCell > ncells)then
                vertexth1 = j
                vertexth2 = j+1
                if(vertexth1 == 4)vertexth2 = 1
                xc1 = xy_coor(cellset(indexCell).nodes(vertexth1),1) + xy_coor(cellset(indexCell).nodes(vertexth2),1)
                yc1 = xy_coor(cellset(indexCell).nodes(vertexth1),2) + xy_coor(cellset(indexCell).nodes(vertexth2),2)
                do k = 1,4
                    vertexth1 = k
                    vertexth2 = k+1
                    if(vertexth1 == 4)vertexth2 = 1
                    xc2 = xy_coor(cellset(indexNearCell).nodes(vertexth1),1) + xy_coor(cellset(indexNearCell).nodes(vertexth2),1)
                    yc2 = xy_coor(cellset(indexNearCell).nodes(vertexth1),2) + xy_coor(cellset(indexNearCell).nodes(vertexth2),2)
                    if(xc1 == xc2 .and. yc1 == yc2)then
                        cellset(indexNearCell).nearcells(k) = indexCell
                        cellset(indexNearCell).sides(k) = indexSide
                    end if                    
                end do                 
            end if           
        end do
    end do

    ! 循环求解虚拟边界单元Jacobi和解点全局坐标
    do i = 1,nbdsides   
        indexCell = ncells + i

        ! 取出四边形单元顶点坐标值
        do j = 1,4
            do k = 1,2
                quad_vertex(j,k) = xy_coor(cellset(indexCell).nodes(j),k)
            end do
        end do
    
        ! 求内部解点坐标
        allocate(cellset(indexCell).sp_coor(nsp,nsp,2))  
        do j = 1,nsp
            do k =1,nsp
                call quad_C2Phy(quad_vertex,SPs_local(j,k,:),cellset(indexCell).sp_coor(j,k,:))
            end do
        end do        
   
        ! 求解点Jacobi 矩阵  [xkesi,ykesi;xeta,yeta]
        allocate(cellset(indexCell).MJacobi(nsp,nsp,4),cellset(indexCell).Mdirect(nsp,nsp,4),cellset(indexCell).det_J(nsp,nsp))
        do j = 1,nsp
            do k = 1,nsp
                call solve_Jacobi(quad_vertex,SPs_local(j,k,:),cellset(indexCell).MJacobi(j,k,:),cellset(indexCell).Mdirect(j,k,:),cellset(indexCell).det_J(j,k))  
            end do
        end do
        
        ! 通量点Jacobi 矩阵  [xkesi,ykesi;xeta,yeta]
        allocate(cellset(indexCell).fpMdirect_F(nsp,nsp+1,4),cellset(indexCell).fpMdirect_G(nsp,nsp+1,4),cellset(indexCell).fpdet_J_F(nsp,nsp+1),cellset(indexCell).fpdet_J_G(nsp,nsp+1))
        do j = 1,nsp
            coor_C(2) = SPs(j)
            do k = 1,nsp+1
                coor_C(1) = FPs(k)               
                call solve_Jacobi(quad_vertex,coor_C(:),M_Jacobi,cellset(indexCell).fpMdirect_F(j,k,:),cellset(indexCell).fpdet_J_F(j,k))  
            end do 
        end do
        do j = 1,nsp
            coor_C(1) = SPs(j)
            do k = 1,nsp+1
                coor_C(2) = FPs(k)               
                call solve_Jacobi(quad_vertex,coor_C(:),M_Jacobi,cellset(indexCell).fpMdirect_G(j,k,:),cellset(indexCell).fpdet_J_G(j,k))         
            end do 
        end do  
    end do
    
    !-----------------------------------------------------------------------------
    !
    !   描述：打印虚拟单元网格，FEPOINT形式
    !
    !-----------------------------------------------------------------------------

    open(200,status='REPLACE',file = '.\mesh\boundCell.plt')
    write(200,*)'title="boundMesh"'
    write(200,*) 'variables =x,y'
    write(200,*) 'ZONE, N=',nbdsides*4,', E=',nbdsides,', F=FEPOINT, ET=QUADRILATERAL'
    
    do i = 1,nbdsides   
        indexCell = ncells + i
        do j = 1,4
            write(200,*) xy_coor(cellset(indexCell).nodes(j),:)
        end do        
    end do
    
    do i = 1,4*nbdsides,4   
        indexCell = ncells + i        
        write(200,*) i,i+1,i+2,i+3     
    end do
    
    close(200)
    
end subroutine pre_BoundaryCell_fluent
    
subroutine sym_Point(x1,y1,x2,y2,xa,ya,x,y)    
 
    !-----------------------------------------------------------------------------
    !
    !   方法：求解对称点
    !   描述：1.方向垂直; 2.中点在直线上. 解二元一次方程组
    !         (x1,y1),(x2,y2)定义直线
    !         (xa,ya)为已知点,(x,y)为待求点
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none
    
    real(prec) :: x1,y1,x2,y2,xa,ya,x,y
    real(prec) :: vector12
    
    if(abs(x1-x2)<1.0e-14)then
        x = 2.0_prec*x1 - xa
        y = ya
    elseif(abs(y1-y2)<1.0e-14)then
        x = xa
        y = 2.0_prec*y1 - ya
    else
        vector12 = (y2-y1)/(x2-x1)
        x = (vector12**2*(2.0_prec*x2-xa) + xa + vector12*2.0_prec*(ya-y2))/(vector12**2+1.0_prec)
        y = -((x-xa)*(x2-x1))/(y2-y1) + ya     
    end if    
    
end subroutine
    
subroutine vertical_if(x1,y1,x2,y2,xa,ya,xb,yb,res)    
 
    !-----------------------------------------------------------------------------
    !
    !   方法：判断直线是否相互垂直
    !   描述：1.方向垂直; 2.中点在直线上. 解二元一次方程组
    !         (x1,y1),(x2,y2)定义直线
    !         (xa,ya)为已知点,(x,y)为待求点
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none
    
    real(prec) :: x1,y1,x2,y2,xa,ya,xb,yb
    real(prec) :: dot_take
    integer :: res
    
    dot_take = (x1-x2)*(xa-xb)+(y1-y2)*(ya-yb)
    if(abs(dot_take) <1.0e-14)then
        res = 1
    else
        res = 0
    end if
    
end subroutine
    
subroutine sovle_sym_point(x1,y1,x2,y2,xa,ya,indexNearCell,varCellNode,varNearCell)
 
    !-----------------------------------------------------------------------------
    !
    !   方法：滑移边界条件
    !   描述：计算对称点处参数。沿垂直壁面方向速度大小相等，方向相反
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use global_var 
    use type_module
    use parameter_setting
    implicit none
    
    real(prec) :: x1,y1,x2,y2,xa,ya,xa_sym,ya_sym,x,y
    integer :: k,l,indexNearCell
    real(prec) :: varCellNode(4),varNearCell(nsp,nsp,4)
    real(prec) :: theta_x12,theta_VaX,theta_symVaX      ! 点1 2与x轴 的夹角，点A的速度与x轴的夹角，对称点A'与x轴的夹角
    real(prec) :: r,u,v,p,r_sym,u_sym,v_sym,p_sym,VV
    
    r = varCellNode(1)
    u = varCellNode(2)
    v = varCellNode(3)
    p = varCellNode(4)
    
    !-----------------------------------------------------------------------------
    !
    !   描述：夹角。从几何关系推导
    !
    !-----------------------------------------------------------------------------
    
    VV = sqrt(u**2 + v**2)
    theta_x12 = atan((y2-y1),(x2-x1))
    theta_VaX = atan(v,u)
    theta_symVaX = 2.0_prec*theta_x12 - theta_VaX       ! theta_symVaX = -(theta_VaX - theta_x12) + theta_x12

    r_sym = r
    u_sym = VV * cos(theta_symVaX)
    v_sym = VV * sin(theta_symVaX)
    p_sym = p

    !-----------------------------------------------------------------------------
    !
    !   描述：寻找相邻的虚拟边界单元内对应的求解点。中点在直线上，垂直向量乘积为0
    !
    !-----------------------------------------------------------------------------
     
    call sym_Point(x1,y1,x2,y2,xa,ya,x,y)
    loop1:do k = 1,nsp
        loop2:do l = 1,nsp
            xa_sym = cellset(indexNearCell).sp_coor(k,l,1)
            ya_sym = cellset(indexNearCell).sp_coor(k,l,2)         
 
            if(abs(xa_sym-x)<1.0e-14 .and. abs(ya_sym-y)<1.0e-14)then
                varNearCell(k,l,1) = r_sym
                varNearCell(k,l,2) = u_sym
                varNearCell(k,l,3) = v_sym
                varNearCell(k,l,4) = p_sym   
                exit loop1
            end if       
        end do loop2
    end do loop1

end subroutine sovle_sym_point
    
subroutine sovle_sym_point2(x1,y1,x2,y2,xa,ya,indexNearCell,varCellNode,varNearCell)

    !-----------------------------------------------------------------------------
    !
    !   方法：另一种滑移边界条件处理方法--不确定对不对，最后未完全验证
    !   描述：计算对称点处参数。沿垂直壁面方向速度大小相等，方向相反
    !   参考：刘巍，张理论，王勇献，邓小刚.计算空气动力学的并行编程基础[M].北京：国防工业出版社，2013: p116
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use global_var 
    use type_module
    use parameter_setting
    implicit none
    
    real(prec) :: x1,y1,x2,y2,xa,ya,xa_sym,ya_sym,x,y
    integer :: k,l,indexNearCell
    real(prec) :: varCellNode(4),varNearCell(nsp,nsp,4)
    real(prec) :: theta_x12,theta_VaX,theta_symVaX  ! 点1 2与x轴 的夹角，点A的速度与x轴的夹角，对称点A'与x轴的夹角
    real(prec) :: r,u,v,p,r_sym,u_sym,v_sym,p_sym,VV
    
    ! 计算对称点处参数。沿垂直壁面方向速度大小相等，方向相反
    r = varCellNode(1)
    u = varCellNode(2)
    v = varCellNode(3)
    p = varCellNode(4)
    VV = sqrt(u**2 + v**2)
    theta_x12 = atan((y2-y1),(x2-x1))
    theta_VaX = atan(v,u)
    theta_symVaX = 2.0_prec*theta_x12 - theta_VaX  !!!! -(theta_VaX - theta_x12) + theta_x12
    !write(*,*) theta_x12,theta_VaX,theta_symVaX
    r_sym = r
    u_sym = VV * cos(theta_symVaX)
    v_sym = VV * sin(theta_symVaX)
    p_sym = p
    !write(*,*)r,u,v,p
    !write(*,*)r_sym,u_sym,v_sym,p_sym
    ! 寻找相邻的虚拟边界单元内对应的求解点
    ! 中点在直线上，垂直向量乘积为0

    call sym_Point(x1,y1,x2,y2,xa,ya,x,y)
    loop1:do k = 1,nsp
        loop2:do l = 1,nsp
            xa_sym = cellset(indexNearCell).sp_coor(k,l,1)
            ya_sym = cellset(indexNearCell).sp_coor(k,l,2)         
            !write(*,*)xa_sym,ya_sym,abs(xa_sym-x),abs(ya_sym-y)
            if(abs(xa_sym-x)<1.0e-14 .and. abs(ya_sym-y)<1.0e-14)then
                !write(*,*)xa_sym,ya_sym,abs(xa_sym-x),abs(ya_sym-y)
                varNearCell(k,l,1) = r_sym
                varNearCell(k,l,2) = u_sym
                varNearCell(k,l,3) = v_sym
                varNearCell(k,l,4) = p_sym   
                !write(*,*)r_sym,u_sym,v_sym,p_sym
                !stop
                exit loop1
            end if       
        end do loop2
    end do loop1
    !varNearCell(k,l,1) = r
    !varNearCell(k,l,2) = -u
    !varNearCell(k,l,3) = v
    !varNearCell(k,l,4) = p
    !stop
end subroutine sovle_sym_point2    
    
subroutine sovle_outflow_point(nearSideth,near_k,varCellNode,varNearCell)

    use real_precision
    use type_module
    use parameter_setting
    implicit none
    
    integer :: k,l,nearSideth,near_k
    real(prec) :: varCellNode(4),varNearCell(nsp,nsp,4)
                
    if(nearsideth == 1 .or. nearsideth == 3)then
        do l = 1,nsp
            varNearCell(l,near_k,:) = varCellNode(:)                   
        end do 
    elseif(nearsideth == 2 .or. nearsideth == 4)then
        do l = 1,nsp
            varNearCell(near_k,l,:) = varCellNode(:)                   
        end do       
    end if
   
end subroutine sovle_outflow_point

subroutine update_Boundary_Detect
 
    !-----------------------------------------------------------------------------
    !
    !   方法：虚拟边界单元是否为问题单元的更新
    !   描述：双马赫反射的所有虚拟边界单元都设置为问题单元，其他的算例是设置问题边界单元的相邻虚拟单元
    !         边界处理目前使用的是虚拟单元，虚拟单元的设置根据边界条件确定
    !         边界处理也可以考虑直接对边界特殊处理
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i,j,k,l,indexCell,indexNearCell,index_temp,p1,p2,index

    if(case_comp == DoubleMach_case )then! .or. case_comp == HyperCylinder_case .or. case_comp == WingFlow_case
        do i = ncells+1,ncells+nbdsides     
            cellset(i).Beta = 1  
            cellset(i).Beta_line = 1
        end do
    else  
        do i = 1,nbdsides
            indexCell = BoundCells_index_set(i)
            if(cellset(indexCell).Beta==0)cycle
            do j = 1,4              
                indexNearCell = cellset(indexCell).nearcells(j)     
                if(indexNearCell>ncells .and.cellset(indexCell).Beta == 1 )then  
                    cellset(indexNearCell).Beta = 1  
                    cellset(indexNearCell).Beta_line = 1
                end if
            end do
        end do
    end if
    
end subroutine update_Boundary_Detect