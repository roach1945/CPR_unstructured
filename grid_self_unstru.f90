subroutine self_unstru_input

    !-----------------------------------------------------------------------------
    !
    !   方法：网格生成程序
    !   描述：生成简单的结构网格，计算域为矩形，网格单元为矩形，梯形或扰动网格。
    !   作者：gqShi 
    !   历史：添加注释 //2021.12.16
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use parameter_setting
    use global_var
    implicit none
    
    ! 曲网格先不考虑，朱老师那复制的代码，曲网格雅可比矩阵还没算。留待以后补充
    
    select case(grid_type)
        case(grid_straight)
            call straightGrid_generate
        case(grid_curve)
            call curveGrid_generate 
        case default
            write(*,*)'Grid generate error!--grid_self_unstru.f90'
    end select   
    
end subroutine self_unstru_input

subroutine straightGrid_generate

    !-----------------------------------------------------------------------------
    !
    !   方法：直边网格生成模块
    !   描述：矩形网格：直接生成
    !         梯形网格：在矩形网格上从中间划条斜杠分成两个梯形
    !         扰动网格：在矩形网格上给每个顶点一个随机扰动（h/4）
    !   作者：gqShi 
    !   历史：其他类型非结构以后补充 //20201226 
    !         修改注释  //2021.12.16
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    implicit none
         
    !-----------------------------------------------------------------------------
    !
    !   描述：标准单元索引貌似 i j反了          
    !         在此用的是第i行，第j列
    !         (3,1) (3,2) (3,3)          *     *     *   
    !         (2,1) (2,2) (2,3)     -->  *     *     *
    !         (1,1) (1,2) (1,3)          *     *     *
    !
    !-----------------------------------------------------------------------------
     
    select case(cell_shape) !----对单元顶点全局坐标赋值
    case(cell_rect)
        call get_nodes_coor_rect
    case(cell_trap)
        call get_nodes_coor_trap 
    case(cell_dist)
        call get_nodes_coor_rect
        call add_disturbance
    end select
     
    call get_cellset        !----对cell 结构体数组赋值
    call get_sideset        !----对side 结构体数组赋值  
    
end subroutine straightGrid_generate

subroutine get_nodes_coor_rect
 
    !-----------------------------------------------------------------------------
    !
    !   方法：矩形网格生成
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer :: i,j,xindex,yindex
    real(prec) :: dx0,dy0,x0,y0 

    nsdpx = nsdx+1          !----x方向节点数
    nsdpy = nsdy+1          !----y方向节点数

    ! 节点、边、单元数
    nnodes = nsdpx*nsdpy
    nsides = (nsdpx-1)*nsdpy+nsdy*nsdpx
    ncells = (nsdpx-1)*(nsdpy-1)
    nbdsides = (nsdpx + nsdpy - 2)*2
    
    write(*,"(' nodes:'I8,'  sides:'I8,'  cells:'I4 '  *'I4,'=',I8)")nnodes,nsides,nsdpx-1,nsdpy-1,ncells
    allocate(xx(1:nsdpx,1:nsdpy),yy(1:nsdpx,1:nsdpy))   
    allocate(xy_coor(nnodes+nbdsides*2,2))
        
    ! 初始化 0
    xy_coor = 0.0_prec
    xlong = xr - xl
    ylong = yr - yl
    dx0 = xlong/real((nsdpx-1),prec)
    dy0 = ylong/real((nsdpy-1),prec)
    
    ! 对单元顶点全局坐标赋值
    do i = 1,nnodes
        xindex = mod(i,nsdpx)
        if(xindex == 0)then
            xindex = nsdpx
        end if
        yindex = int((i+nsdpx-1)/nsdpx)

        xx(xindex,yindex) = dx0*(xindex-1) + xl
        yy(xindex,yindex) = dy0*(yindex-1) + yl 
        
        xy_coor(i,1) = xx(xindex,yindex)
        xy_coor(i,2) = yy(xindex,yindex)
    end do

end subroutine get_nodes_coor_rect
    
subroutine add_disturbance

    !-----------------------------------------------------------------------------
    !
    !   方法：扰动四边形网格添加随机扰动
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer :: i,j,xindex,yindex
    real(prec) :: dx0,dy0,x0,y0,randomNum,r_min,r_max,r_len,stren_disturb
    
    dx0 = xlong/real((nsdpx-1),prec)
    dy0 = ylong/real((nsdpy-1),prec)
    r_min = -1.0_prec
    r_max =  1.0_prec
    r_len = r_max-r_min
    stren_disturb = 0.25_prec
    
    ! 对单元顶点坐标添加随机扰动
    call RANDOM_SEED()
    do i = 1,nnodes
        call RANDOM_NUMBER(randomNum) 

        if(abs(xy_coor(i,1)-xl)<1.0e-14 .OR. abs(xy_coor(i,1)-xr)<1.0e-14 .OR. abs(xy_coor(i,2)-yl)<1.0e-14 .OR. abs(xy_coor(i,2)-yr)<1.0e-14)cycle
        if(abs(xy_coor(i,1)-xl-dx0)<1.0e-14 .OR. abs(xy_coor(i,1)-xr-dx0)<1.0e-14 .OR. abs(xy_coor(i,2)-yl-dy0)<1.0e-14 .OR. abs(xy_coor(i,2)-yr-dy0)<1.0e-14)cycle
        randomNum = r_min + r_len*randomNum
        xy_coor(i,1) = xy_coor(i,1) + dx0*stren_disturb*randomNum
        xy_coor(i,2) = xy_coor(i,2) + dy0*stren_disturb*randomNum

    end do

end subroutine add_disturbance
    
subroutine get_nodes_coor_trap

    !-----------------------------------------------------------------------------
    !
    !   方法：梯形网格生成
    !   描述：梯形网格用处不大，之前是想测试一下不规则网格的表现
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer :: i,j,xindex,yindex
    real(prec) :: dx0,dy0,x0,y0

    nsdx = nsdx/2
    nsdpx = nsdx*2+1      !----x方向节点数
    nsdpy = nsdy+1        !----y方向节点数

    ! 节点、边、单元数
    nnodes = nsdpx*nsdpy
    nsides = (nsdpx-1)*nsdpy+nsdy*nsdpx
    ncells = (nsdpx-1)*(nsdpy-1)
    nbdsides = (nsdpx + nsdpy - 2)*2

    write(*,"(' nodes:'I8,'  sides:'I8,'  cells:'I8)")nnodes,nsides,ncells
    allocate(xx(1:nsdpx,1:nsdpy),yy(1:nsdpx,1:nsdpy))     
    allocate(xy_coor(nnodes+nbdsides*2,2))
        
    ! 初始化 0
    xy_coor = 0.0_prec
    xlong = xr - xl
    ylong = yr - yl
    dx0 = xlong/real((nsdx),prec)
    dy0 = ylong/real((nsdy),prec)
    
    ! 对单元顶点全局坐标赋值 
    do i = 1,nnodes
        
        xindex = mod(i,nsdpx)
        if(xindex == 0)then
            xindex = nsdpx
        end if       
        yindex = int((i+nsdpx-1)/nsdpx)
        
        if(mod(xindex+1,2) == 0)then
            xx(xindex,yindex) = dx0*((xindex+1)/2-1) + xl
            yy(xindex,yindex) = dy0*(yindex-1) + yl 
        else
            if(mod(yindex,2) == 0)then
                xx(xindex,yindex) = dx0*(xindex/2-1) + dx0*1.0_prec/3.0_prec + xl
            else
                xx(xindex,yindex) = dx0*(xindex/2-1) + dx0*2.0_prec/3.0_prec + xl
            end if  
            yy(xindex,yindex) = dy0*(yindex-1) + yl 
        end if   
        
        xy_coor(i,1) = xx(xindex,yindex)
        xy_coor(i,2) = yy(xindex,yindex)

    end do

    nsdx = nsdx*2

end subroutine get_nodes_coor_trap

subroutine get_cellset

    !-----------------------------------------------------------------------------
    !
    !   方法：获取单元结构体cellset的数据 
    !   描述：包括索引、顶点、侧边、相邻单元、Jacobi矩阵等
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    integer i,j
    integer :: cell_vertex(4),cell_side(4),cell_near(4)

    call allo_stru_record
    do i = 1,ncells
        
        ! 单元编号
        cellset(i).index = i
        
        ! 单元顶点编号
        cell_vertex(1) = i + int(i/(nsdpx-1)); if(mod(i,nsdpx-1) == 0) cell_vertex(1) = cell_vertex(1)-1;
        cell_vertex(2) = cell_vertex(1)+1
        cell_vertex(3) = cell_vertex(1)+1+nsdpx
        cell_vertex(4) = cell_vertex(3)-1
        cellset(i).nodes(:) = cell_vertex(:)
        
        ! 单元侧边编号
        cell_side(1) = mod(i,nsdpx-1) + int(i/(nsdpx-1))*(2*nsdpx-1);if(mod(i,nsdpx-1) == 0) cell_side(1) = nsdpx-1 + (int(i/(nsdpx-1))-1)*(2*nsdpx-1);
        cell_side(2) = cell_side(1) + nsdpx
        cell_side(3) = cell_side(1) + (2*nsdpx-1)
        cell_side(4) = cell_side(1) + (nsdpx-1)
        cellset(i).sides(:) = cell_side(:)

    end do
    
   
    call record_near_cell   !----录入相邻单元    
    call solve_SPs_local    !----求解点Jacobi矩阵和全局坐标值      
    call solve_FPs_local    !----通量点Jacobi矩阵和全局坐标值
    
end subroutine get_cellset

subroutine record_near_cell

    !-----------------------------------------------------------------------------
    !
    !   方法：获取相邻单元
    !   描述：主要思想是查询有相同侧边的单元，若有则记录为相邻单元，若无，则说明是边界处的单元，相邻单元编号记为0
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i,j,k,l
    integer :: cell_vertex(4),cell_side(4),cell_near(4)
    
    ! 遍历，搜索，记录(record) C--cell S--side
    recordedC: do i = 1,ncells
        recordedCS: do k = 1,4
            searchedC: do j = 1,ncells
                if(j == i) cycle searchedC
                searchedCS: do l = 1,4
                    if(cellset(i).sides(k) == cellset(j).sides(l))then
                        cellset(i).nearcells(k) = cellset(j).index
                        cycle recordedCS                !记录i单元k侧边后跳出，进行k+1侧边记录
                    else
                        cellset(i).nearcells(k) = 0     !表示在边界处，无相邻单元，需要进行边界条件设置
                    end if
                end do searchedCS
            end do searchedC
        end do recordedCS

    end do recordedC
    
end subroutine record_near_cell

subroutine solve_SPs_local
 
    !-----------------------------------------------------------------------------
    !
    !   方法：求解点计算
    !   描述：单元顶点坐标、全局坐标、Jacobi矩阵
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i,j,k
    real(prec),dimension(:,:) :: quad_vertex(4,2)
    real(prec),dimension(:,:,:) :: weight_sps(nsp,nsp,4)
    
    
    !循环求解每个单元Jacobi和解点全局坐标
    do i = 1,ncells
    
        ! 取出四边形单元顶点坐标值
        do j = 1,4
            do k = 1,2
                quad_vertex(j,k) = xy_coor(cellset(i).nodes(j),k)
            end do
        end do

        ! 求内部解点坐标
        allocate(cellset(i).sp_coor(nsp,nsp,2)) 
        do j = 1,nsp
            do k =1,nsp
                call quad_C2Phy(quad_vertex,SPs_local(j,k,:),cellset(i).sp_coor(j,k,:))
            end do
        end do        
        
        ! 求解点Jacobi 矩阵  [xkesi,ykesi;xeta,yeta]
        allocate(cellset(i).MJacobi(nsp,nsp,4),cellset(i).Mdirect(nsp,nsp,4),cellset(i).det_J(nsp,nsp))
        do j = 1,nsp
            do k = 1,nsp
                call solve_Jacobi(quad_vertex,SPs_local(j,k,:),cellset(i).MJacobi(j,k,:),cellset(i).Mdirect(j,k,:),cellset(i).det_J(j,k))                   
            end do
        end do
           
    end do 
end subroutine solve_SPs_local
    
subroutine solve_FPs_local
 
    !-----------------------------------------------------------------------------
    !
    !   方法：通量点计算
    !   描述：单元顶点坐标、全局坐标、Jacobi矩阵
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i,j,k
    real(prec),dimension(:,:) :: quad_vertex(4,2)
    real(prec),dimension(:) :: coor_C(2),M_Jacobi(4),M_direct(4)
    real(prec) :: detJ
    
    ! 循环求解每个单元Jacobi和通量点全局坐标
    do i = 1,ncells 
        
        ! 取出四边形单元顶点坐标值
        do j = 1,4
            do k = 1,2
                quad_vertex(j,k) = xy_coor(cellset(i).nodes(j),k)
            end do
        end do
        
        ! 求通量点Jacobi 矩阵  [xkesi,ykesi;xeta,yeta]
        allocate(cellset(i).fpMdirect_F(nsp,nsp+1,4),cellset(i).fpMdirect_G(nsp,nsp+1,4),cellset(i).fpdet_J_F(nsp,nsp+1),cellset(i).fpdet_J_G(nsp,nsp+1))!分配内存给MJacobi
        do j = 1,nsp
            coor_C(2) = SPs(j)
            do k = 1,nsp+1
                coor_C(1) = FPs(k)               
                call solve_Jacobi(quad_vertex,coor_C(:),M_Jacobi,cellset(i).fpMdirect_F(j,k,:),cellset(i).fpdet_J_F(j,k))  
            end do 
        end do
        do j = 1,nsp
            coor_C(1) = SPs(j)
            do k = 1,nsp+1
                coor_C(2) = FPs(k)               
                call solve_Jacobi(quad_vertex,coor_C(:),M_Jacobi,cellset(i).fpMdirect_G(j,k,:),cellset(i).fpdet_J_G(j,k))         
            end do 
        end do
        
    end do 

end subroutine solve_FPs_local
    
subroutine quad_C2Phy(cell_vertex,coor_C,coor_P)
 
    !-----------------------------------------------------------------------------
    !
    !   方法：映射求解全局坐标
    !   描述：根据计算域坐标，可求解直边网格中对应的求解点、通量点的物理域坐标
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    implicit none

    real(prec),dimension(:)   :: weight(4)
    real(prec),dimension(:,:) :: cell_vertex(4,2)   !----实际单元顶点，标准单元已知
    real(prec),dimension(:)   :: coor_C(2),coor_P(2)!----标准单元，实际单元
    
    ! 根据标准单元(-1,-1)(1,1)求解点处的权重
    weight(1) = 1.0_prec - 0.5_prec*(coor_C(2)+1)
    weight(2) = 0.5_prec*(coor_C(2)+1)
    weight(3) = 1.0_prec - 0.5_prec*(coor_C(1)+1)
    weight(4) = 0.5_prec*(coor_C(1)+1) 
    
    ! x
    coor_P(1) = 0.5_prec*(weight(1)*((cell_vertex(2,1)-cell_vertex(1,1))*coor_C(1)+((cell_vertex(2,1)+cell_vertex(1,1))))&
              + weight(2)*((cell_vertex(3,1)-cell_vertex(4,1))*coor_C(1)+((cell_vertex(3,1)+cell_vertex(4,1)))))    
    ! y
    coor_P(2) = 0.5_prec*(weight(3)*((cell_vertex(4,2)-cell_vertex(1,2))*coor_C(2)+((cell_vertex(4,2)+cell_vertex(1,2))))&
              + weight(4)*((cell_vertex(3,2)-cell_vertex(2,2))*coor_C(2)+((cell_vertex(3,2)+cell_vertex(2,2)))))   
                                
end subroutine quad_C2Phy

subroutine solve_Jacobi(cell_vertex,coor_C,M_Jacobi,M_direct,detJ)
 
    !-----------------------------------------------------------------------------
    !
    !   方法：Jacobi 矩阵计算
    !   描述：单元变换时按照侧边伸缩比进行变换。双线性变换
    !         具体实施方案：
    !         xkesi = Lx/2; yeta = Ly/2
    !         x = w1((x1-x0)/2*kesi - (x0+x1)/2) + w2((x2-x3)/2*kesi - (x2+x3)/2)
    !         y = w3((y3-y0)/2*eta - (y0+y3)/2) + w4((y2-y1)/2*eta - (y2+y1)/2)
    !         均匀分布，wi可从 标准四边形A(-1,-1)(1,-1)(1,1)(-1,1)得到,w1 = 1-(y'+1)/2;w2 = (y'+1)/2;w3 = 1-(x'+1)/2;w4 = (x'+1)/2
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use real_precision
    use parameter_setting,only:nsp
    implicit none

    real(prec),dimension(:,:) :: cell_vertex(4,2)                   !----实际单元顶点，标准单元已知
    real(prec),dimension(:)   :: coor_C(2),M_Jacobi(4),M_direct(4)  !----标准单元,Jacobi矩阵
    real(prec) :: detJ
    
    !!xkesi
    M_Jacobi(1) = (1.0_prec-0.5_prec*(coor_C(2)+1.0_prec))*(cell_vertex(2,1)-cell_vertex(1,1))*0.5_prec&
                + 0.25_prec*(coor_C(2)+1.0_prec)*(cell_vertex(3,1)-cell_vertex(4,1))
    !!xeta
    M_Jacobi(2) = -0.25_prec*((cell_vertex(2,1)-cell_vertex(1,1))*coor_C(1)+(cell_vertex(2,1)+cell_vertex(1,1)))&
                + 0.25_prec*((cell_vertex(3,1)-cell_vertex(4,1))*coor_C(1)+(cell_vertex(3,1)+cell_vertex(4,1)))   
    !!ykesi
    M_Jacobi(3) = -0.25_prec*((cell_vertex(4,2)-cell_vertex(1,2))*coor_C(2)+(cell_vertex(4,2)+cell_vertex(1,2)))&
                + 0.25_prec*((cell_vertex(3,2)-cell_vertex(2,2))*coor_C(2)+(cell_vertex(3,2)+cell_vertex(2,2)))         
    !!yeta
    M_Jacobi(4) = (1.0_prec-0.5_prec*(coor_C(1)+1.0_prec))*(cell_vertex(4,2)-cell_vertex(1,2))*0.5_prec&
                + 0.25_prec*(coor_C(1)+1.0_prec)*(cell_vertex(3,2)-cell_vertex(2,2))    
    
    ! detJ
    detJ = M_Jacobi(1)*M_Jacobi(4)-M_Jacobi(2)*M_Jacobi(3)
    
    ! kesi_x
    M_direct(1) = M_Jacobi(4)/detJ            
    ! kesi_y
    M_direct(2) = -M_Jacobi(2)/detJ               
    ! eta_x
    M_direct(3) = -M_Jacobi(3)/detJ               
    ! eta_y
    M_direct(4) = M_Jacobi(1)/detJ               
                
end subroutine solve_Jacobi

subroutine solve_Jacobi2(cell_vertex,coor_C,M_Jacobi,M_direct,detJ)
    !。另一种计算方法
     
    !-----------------------------------------------------------------------------
    !
    !   方法：Jacobi 矩阵计算
    !   描述：等参单元方法。有限元方法里的内容，自行搜索“等参单元”
    !         或参考 2009，Wang Z J，LCP; sgq 硕士论文
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit none

    real(prec),dimension(:,:) :: cell_vertex(4,2)                           !实际单元顶点，标准单元已知
    real(prec),dimension(:)   :: coor_C(2),M_Jacobi(4),M_direct(4)          !标准单元,Jacobi矩阵
    real(prec),dimension(:)   :: cell_vertexL(4,2)
    real(prec) :: detJ
    
    cell_vertexL(1,1) = -1.0_prec;cell_vertexL(1,2) = -1.0_prec;
    cell_vertexL(2,1) =  1.0_prec;cell_vertexL(2,2) = -1.0_prec;
    cell_vertexL(3,1) =  1.0_prec;cell_vertexL(3,2) =  1.0_prec;
    cell_vertexL(4,1) = -1.0_prec;cell_vertexL(4,2) =  1.0_prec;

    !xkesi
    call solve_M_Jacobi(coor_C(2),cell_vertexL(:,1),cell_vertexL(:,2),cell_vertex(:,1),M_Jacobi(1))
    
    !xeta
    call solve_M_Jacobi(coor_C(1),cell_vertexL(:,2),cell_vertexL(:,1),cell_vertex(:,1),M_Jacobi(2))
    
    !ykesi
    call solve_M_Jacobi(coor_C(2),cell_vertexL(:,1),cell_vertexL(:,2),cell_vertex(:,2),M_Jacobi(3))

    !yeta
    call solve_M_Jacobi(coor_C(1),cell_vertexL(:,2),cell_vertexL(:,1),cell_vertex(:,2),M_Jacobi(4))

    !detJ
    detJ = M_Jacobi(1)*M_Jacobi(4)-M_Jacobi(2)*M_Jacobi(3)
    
    !kesi_x
    M_direct(1) = M_Jacobi(4)/detJ            
    !kesi_y
    M_direct(2) = -M_Jacobi(2)/detJ               
    !eta_x
    M_direct(3) = -M_Jacobi(3)/detJ               
    !eta_y
    M_direct(4) = M_Jacobi(1)/detJ               
                
end subroutine solve_Jacobi2
    
subroutine solve_M_Jacobi(coorL_xy,vertexL1_xy,vertexL2_xy,vertexP_xy,M_Jacobi_th)

    use real_precision
    implicit none
    
    real(prec) :: coorL_xy,vertexL1_xy(4),vertexL2_xy(4),vertexP_xy(4),M_Jacobi_th
    integer :: i
    
    M_Jacobi_th = 0.0_prec
    do i = 1,4
        M_Jacobi_th = M_Jacobi_th + 0.25_prec*vertexL1_xy(i)*(1.0_prec + vertexL2_xy(i)*coorL_xy)*vertexP_xy(i)         
    end do
    
end subroutine solve_M_Jacobi
    
subroutine get_sideset

    !-----------------------------------------------------------------------------
    !
    !   方法：获取边结构体sideset的数据 
    !   描述：包括顶点、相邻单元等
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i,j,indexSide
    integer :: cell_vertex(4),cell_side(4),cell_near(4)
    
    do i = 1,nsides
        sideset(i).nodes(:) = 0
        sideset(i).nearcells(:) = 0
    end do

    do i = 1,ncells
        do j = 1,4
            indexSide = cellset(i).sides(j)
            if(sideset(indexSide).nearcells(1)==0)then !----首先存入左单元，之后存入右单元
                sideset(indexSide).nearcells(1) = i   
                !起点，终点
                if(j==4)then
                    sideset(indexSide).nodes(1) = cellset(i).nodes(4)
                    sideset(indexSide).nodes(2) = cellset(i).nodes(1)
                else
                    sideset(indexSide).nodes(1) = cellset(i).nodes(j)
                    sideset(indexSide).nodes(2) = cellset(i).nodes(j+1)
                end if
            else
                sideset(indexSide).nearcells(2) = i  
            end if
        end do    
    end do
    
end subroutine get_sideset

subroutine Matrix_Inverse(A,n,B)
     
    !-----------------------------------------------------------------------------
    !
    !   方法：矩阵求逆
    !   描述：修改网上撸的代码：https://blog.csdn.net/zxylv/article/details/87028584 
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    implicit None
    
    integer :: n
    integer :: i
    real(prec) :: A(n,n),B(n,n)
    
    b=inv(a,n)
    
    contains
 
    ! 递归求行列式的值
    recursive function det(A,col,row) result(D)
        Implicit None
        integer row,col
        real(prec)::A(col,row),B(row-1,col-1)
        real(prec)::D
        integer row_now,col_now,k,c,f
        row_now=row 
        col_now=col
        if (row_now>1) then
            D = 0.0;
            do k=1,row_now
                if(k==1)then
                    B=A(2:col_now,2:row_now)
                elseif(k<row_now)then
                    B(1:k-1,:)=A(1:k-1,2:row_now)
                    B(k:col_now-1,:)=A(k+1:col_now,2:row_now)
                else
                    B=A(1:col_now-1,2:row_now)
                endif
                c=col-1
                f=row-1
                D = D + (-1)**(1+k) * A(k,1) *&
                    det(B,c,f)
            end do
        else
            D = A(1,1);
        end if
    end function det
 
    function delete_col_row(A,m,n,k,tar) result(b)
        !tar == 1，在 m*n 的矩阵中删去第 k 行
        !tar != 1，在 m*n 的矩阵中删去第 k 列,方便伴随矩阵的计算
        Implicit None
        integer::m,n,k,tar
        real(prec)::a(m,n)
        real(prec),allocatable::B(:,:)
        if(tar==1)then
            allocate(B(m-1,n))
            if(k==1)then
                B=A(2:m,:)
            elseif(k<m)then
                B(1:k-1,:)=A(1:k-1,:)
                B(k:m-1,:)=A(k+1:m,:)
            else
                B=A(1:m-1,:)
            endif
        else
            allocate(B(m,n-1))
            if(k==1)then
                B=A(:,2:n)
            elseif(k<n)then
                B(:,1:k-1)=A(:,1:k-1)
                B(:,k:n-1)=A(:,k+1:n)
            else
                B=A(:,1:n-1)
            endif
        endif
    end function delete_col_row
 
    !用  A* / |A|  计算逆矩阵
    function inv(A,row) 
        Implicit None
        integer::row,i,j
        real(prec)::A(row,row)
        real(prec)::inv(row,row),b(row-1,row-1)
        do i=1,row
            do j=1,row
                b=delete_col_row(delete_col_row(A,row,row,i,1),row-1,row,j,2)
                inv(i,j)=(-1)**(i+j)*det(b,row-1,row-1)/det(A,row,row)
            end do
        end do
        inv=transpose(inv)
    end function inv
    
end subroutine Matrix_Inverse

subroutine curveGrid_generate
    !结构网格，曲网格，非结构储存；结构网格是特殊的非结构网格，，，
    use real_precision
    use global_var
    implicit none
    integer :: i,j,ex
    real(prec) :: dx0,dy0,x0,x1,x2,y0,y1,y2,ds,aa0
    character(len=20)filename
    !调整曲网格曲率
    ni = nsdx+1
    nj = nsdy+1 
    allocate(xx(1:ni,1:nj),yy(1:ni,1:nj))
    
    xlong = xr - xl
    ylong = yr - yl
    ex = 10
    aa0 = 1.0_prec!
    dx0 = xlong/real((ni-1),prec)
    dy0 = ylong/real((nj-1),prec)
    
    do j=1,nj
        do i=1,ni
            x0 = dx0*real((i-1),prec) + xl
            y0 = dy0*real((j-1),prec) + yl

            x1 = pi*x0 + 0.5_prec*pi
            y1 = 2.0_prec*pi*y0 + pi
            x2 = sin(x1)**2.0_prec
            y2 = sin(y1)**2.0_prec
        
            xx(i,j) = x0 
            ds=8.0_prec
            if ( abs(x0) .LT. ds ) then
                y1 = 5.0_prec*pi*x0/(10.0_prec*ds) + 0.5_prec*pi
                y2 = sin(y1)**ex
            else 
                y2 = 0.0_prec
            end if

            y1 = pi*y0/ylong + 0.5_prec*pi
            x2 = sin(y1)**2.0_prec
            
            yy(i,j) = y0 + aa0*x2*y2
        end do
    end do
 
    !输出曲边网格 mesh.plt查看效果
    filename = 'mesh_curve.plt'
    open(15,file = filename)
    write(15,*)'variables =x,y'
    write(15,*)'zone i=',ni,'j=',nj,'   ZONETYPE = Ordered' 
    write(15,*) ((xx(i,j),i=1,ni),j=1,nj)
    write(15,*) ((yy(i,j),i=1,ni),j=1,nj)
    close(15)
                     
    !之后的编号以后进行
    !call 
end subroutine curveGrid_generate