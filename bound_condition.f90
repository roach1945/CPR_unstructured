subroutine set_BC_fluent

    !-----------------------------------------------------------------------------
    !
    !   ����������fluent�����ļ����ñ߽�����
    !   ��������������ͬ�߽��������߽�������pointweise����.
    !         ���ݱ߽������趨��ӵ����ⵥԪ��ָ��λ�á�����ָ��Ŀ�ģ����ڼ���
    !   ���ߣ�gqShi 
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
    !   �������������������������ӱ߽�����
    !
    !-----------------------------------------------------------------------------
     
    if(grid_set == self)then
        call get_bound_Name
    endif
    
    ! ��¼��ߵ�Ԫ������
    allocate(BoundCells_index_set(nbdsides))      
    count_temp = 1
    
             
    !-----------------------------------------------------------------------------
    !
    !   ������������߽絥Ԫ��š�������ѭ��һ�飬�սǴ�����������ڱ߽����
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
    !   ��������������߽絥Ԫ��ָ��Ҳ��index��ѭ���߽�ָ��ѭ����Ӧλ�ã�����ָ�����ڱ߽絥Ԫ
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
                !   ������Ѱ��ѭ�����Ӧ�ߡ��ڴ�ʹ�ñߵ����ĵ�ȶ�  ##����ѭ��##    
                !
                !-----------------------------------------------------------------------------
             
                ! ����Ԫ�߽�����ĵ�  
                x1 = (xy_coor(sideset(indexSide).nodes(1),1)+xy_coor(sideset(indexSide).nodes(2),1))*0.5_prec  
                y1 = (xy_coor(sideset(indexSide).nodes(1),2)+xy_coor(sideset(indexSide).nodes(2),2))*0.5_prec  
                
                ! ��Ӧ�����ĵ�
                Loop03:do k = 1,nbdsides
                    indexCell2 = BoundCells_index_set(k)
                    Loop04:do l = 1,4                            
                        indexNearCell2 = cellset(indexCell2).nearcells(l)
                        indexSide2 = cellset(indexCell2).sides(l)

                        ! ��Ӧ��Ԫ�߽������ĵ�    
                        if(sideset(indexSide2).bc == Periodic .and. indexSide2 /= indexSide)then
                            
                            x2 = (xy_coor(sideset(indexSide2).nodes(1),1)+xy_coor(sideset(indexSide2).nodes(2),1))*0.5_prec
                            y2 = (xy_coor(sideset(indexSide2).nodes(1),2)+xy_coor(sideset(indexSide2).nodes(2),2))*0.5_prec
                            
                            dis_p_c = sqrt((x1-x2)**2 + (y1-y2)**2)

                            ! ���е�������غϣ����Ӧ������������бһ���Ƕ��޷�����
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
    !   ����������������߽���������
    !   �����������������б߽��������á�����������һ�������μ�����
    !   ���ߣ�gqShi 
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
    !   ���������±߽�
    !   ������ÿһʱ�䲽��ÿһʱ��㶼��Ҫ���±߽磬��ͬ�����������в�ͬ������
    !         Wall�߽����ñȽ��鷳����ȷ���Էǽṹ������ԣ������Ƿ���ȷ�� 
    !   ���ߣ�gqShi 
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
        
        ! ˫�������
        do i = 1,nbdsides    
            indexCell = BoundCells_index_set(i)
            do j = 1,4     
                indexSide = cellset(indexCell).sides(j)
                indexNearCell = cellset(indexCell).nearcells(j)   
                if(sideset(indexSide).bc == Inlet_vent)then
                    
                    ! ��߽�
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
                    
                    ! �ϱ߽� �� �ұ߽�
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
                    
                    ! �±߽�
                    vertexth1 = j
                    vertexth2 = j+1
                    if(vertexth1 == 4) vertexth2 = 1
                    x1 = xy_coor(cellset(indexCell).nodes(vertexth1),1)
                    y1 = xy_coor(cellset(indexCell).nodes(vertexth1),2)
                    x2 = xy_coor(cellset(indexCell).nodes(vertexth2),1)
                    y2 = xy_coor(cellset(indexCell).nodes(vertexth2),2)
                    do k = 1,nsp
                        do l = 1,nsp
                            
                            ! (x1,y1),(x2,y2)�Ǳ߽�ߵĶ˵㡣ȷ������
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
        ! �����������ձ߽���������
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
                            
                            !�߽總�������㣬ȷ�����ⵥԪ�ڵ���֮���ߺͱ߽紹ֱ�ĵ�
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
                            
                            !�߽總�������㣬ȷ�����ⵥԪ�ڵ���֮���ߺͱ߽紹ֱ�ĵ�
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
                    
                    !����߽絥Ԫindex��ָ����ͬ�����ڱ߽�
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
                            ! (x1,y1),(x2,y2)�Ǳ߽�ߵĶ˵㡣ȷ������
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
    !   ����������������������������߽�����
    !   ���ߣ�gqShi 
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
    !   �������������߽絥Ԫ��Ϣ
    !   ���������õ�Ԫ���꣬������ʵ���꣬Jacobi�Ȳ���ʱ��仯��ֵ.��Ԫ��״�����ڵ�Ԫ�ر߽��߶Գơ�
    !         ����ѭ���߽�������case֮��������¸�ֵ������صı��������ص��ĵ�Ԫ��״����������˵��
    !   ���ߣ�gqShi 
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
    !   �������������ⵥԪ����״���ڱ߽���߽絥Ԫ�Գ�
    !
    !-----------------------------------------------------------------------------
     
    pCounter = nnodes+1
    Loop1:do i = 1,nbdsides
        
        ! �߽������ڵ�Ԫ
        indexCell = BoundCells_index_set(i)
        do j = 1,i-1
            if(BoundCells_index_set(j) == indexCell)cycle Loop1
        end do
        
        do j = 1,4                     
            ! �߽������ڵ�Ԫ���ڵ�Ԫ
            indexNearCell = cellset(indexCell).nearcells(j)
            if(indexNearCell > ncells)then    
                
                ! ��¼���غϵ�����.�����֮ǰ�Ĺ������±��
                vertexth(1) = j
                vertexth(2) = j+1
                if(j == 4) vertexth(2) = 1
                cellset(indexNearCell).nodes(1) = cellset(indexCell).nodes(vertexth(1))
                cellset(indexNearCell).nodes(2) = cellset(indexCell).nodes(vertexth(2))
                indexNode = 3
                
                ! �����ֱ�ߵĶԳƵ�
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
    !   �������������߽絥Ԫ����Ϣ�����ڲ���Ԫ��Ϣһ��
    !
    !-----------------------------------------------------------------------------
    
    do i = 1,nbdsides
        indexCell = ncells + i
        do j =1,4
            coor_vertex(j,:) = xy_coor(cellset(indexCell).nodes(j),:)
        end do

        ! ��ֵ����
        call sort_qua(coor_vertex(:,:))

        ! ��ֵ�������,��ʱ��
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
        
        ! ��������߽絥Ԫ�����߽絥Ԫ���ڵ�Ԫ�ͱ����
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

    ! ѭ���������߽絥ԪJacobi�ͽ��ȫ������
    do i = 1,nbdsides   
        indexCell = ncells + i

        ! ȡ���ı��ε�Ԫ��������ֵ
        do j = 1,4
            do k = 1,2
                quad_vertex(j,k) = xy_coor(cellset(indexCell).nodes(j),k)
            end do
        end do
    
        ! ���ڲ��������
        allocate(cellset(indexCell).sp_coor(nsp,nsp,2))  
        do j = 1,nsp
            do k =1,nsp
                call quad_C2Phy(quad_vertex,SPs_local(j,k,:),cellset(indexCell).sp_coor(j,k,:))
            end do
        end do        
   
        ! ����Jacobi ����  [xkesi,ykesi;xeta,yeta]
        allocate(cellset(indexCell).MJacobi(nsp,nsp,4),cellset(indexCell).Mdirect(nsp,nsp,4),cellset(indexCell).det_J(nsp,nsp))
        do j = 1,nsp
            do k = 1,nsp
                call solve_Jacobi(quad_vertex,SPs_local(j,k,:),cellset(indexCell).MJacobi(j,k,:),cellset(indexCell).Mdirect(j,k,:),cellset(indexCell).det_J(j,k))  
            end do
        end do
        
        ! ͨ����Jacobi ����  [xkesi,ykesi;xeta,yeta]
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
    !   ��������ӡ���ⵥԪ����FEPOINT��ʽ
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
    !   ���������ԳƵ�
    !   ������1.����ֱ; 2.�е���ֱ����. ���Ԫһ�η�����
    !         (x1,y1),(x2,y2)����ֱ��
    !         (xa,ya)Ϊ��֪��,(x,y)Ϊ�����
    !   ���ߣ�gqShi 
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
    !   �������ж�ֱ���Ƿ��໥��ֱ
    !   ������1.����ֱ; 2.�е���ֱ����. ���Ԫһ�η�����
    !         (x1,y1),(x2,y2)����ֱ��
    !         (xa,ya)Ϊ��֪��,(x,y)Ϊ�����
    !   ���ߣ�gqShi 
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
    !   ���������Ʊ߽�����
    !   ����������ԳƵ㴦�������ش�ֱ���淽���ٶȴ�С��ȣ������෴
    !   ���ߣ�gqShi 
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
    real(prec) :: theta_x12,theta_VaX,theta_symVaX      ! ��1 2��x�� �ļнǣ���A���ٶ���x��ļнǣ��ԳƵ�A'��x��ļн�
    real(prec) :: r,u,v,p,r_sym,u_sym,v_sym,p_sym,VV
    
    r = varCellNode(1)
    u = varCellNode(2)
    v = varCellNode(3)
    p = varCellNode(4)
    
    !-----------------------------------------------------------------------------
    !
    !   �������нǡ��Ӽ��ι�ϵ�Ƶ�
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
    !   ������Ѱ�����ڵ�����߽絥Ԫ�ڶ�Ӧ�����㡣�е���ֱ���ϣ���ֱ�����˻�Ϊ0
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
    !   ��������һ�ֻ��Ʊ߽�����������--��ȷ���Բ��ԣ����δ��ȫ��֤
    !   ����������ԳƵ㴦�������ش�ֱ���淽���ٶȴ�С��ȣ������෴
    !   �ο�����Ρ�������ۣ������ף���С��.�����������ѧ�Ĳ��б�̻���[M].������������ҵ�����磬2013: p116
    !   ���ߣ�gqShi 
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
    real(prec) :: theta_x12,theta_VaX,theta_symVaX  ! ��1 2��x�� �ļнǣ���A���ٶ���x��ļнǣ��ԳƵ�A'��x��ļн�
    real(prec) :: r,u,v,p,r_sym,u_sym,v_sym,p_sym,VV
    
    ! ����ԳƵ㴦�������ش�ֱ���淽���ٶȴ�С��ȣ������෴
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
    ! Ѱ�����ڵ�����߽絥Ԫ�ڶ�Ӧ������
    ! �е���ֱ���ϣ���ֱ�����˻�Ϊ0

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
    !   ����������߽絥Ԫ�Ƿ�Ϊ���ⵥԪ�ĸ���
    !   ������˫��շ������������߽絥Ԫ������Ϊ���ⵥԪ����������������������߽絥Ԫ���������ⵥԪ
    !         �߽紦��Ŀǰʹ�õ������ⵥԪ�����ⵥԪ�����ø��ݱ߽�����ȷ��
    !         �߽紦��Ҳ���Կ���ֱ�ӶԱ߽����⴦��
    !   ���ߣ�gqShi 
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