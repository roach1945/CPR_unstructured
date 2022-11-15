subroutine pre_processing

    write(*,*)'------------------------------------------------------'  
    write(*,*)'------------------------------------------------------'
    
    call read_parameter         !----读取控制参数

    call pre_setting            !----根据读取的参数计算其他参数 

    call grid_generate          !----网格生成   
    
    call fluid_initialize       !----初始化流场 
    
    call get_min_dis            !----求最小全局横纵坐标变换值，求dt
        
    call print_submesh_data     !----网格

    call print_mesh_data        !----输出plt格式网格
    
    open(unit=1010,file='error.dat')
    
end subroutine pre_processing 

subroutine read_parameter

    !-----------------------------------------------------------------------------
    !
    !   描述：读取预设参数. 输入文件 input.in 
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use preParameter
    use global_var,only: T
    use parameter_setting
    implicit none
    
    character(len = 20)::filename   
    character(len = 80)::msg       
    integer    :: status1           
    real(prec) :: nsdx_temp     
    
    filename = 'input.in'  
    write(*,*) 'Read parameters.'     
    write(*,"('     parameters in file: 'A)") filename
    open(UNIT = 10 , FILE = filename , STATUS = 'OLD' , ACTION = 'read' , IOSTAT = status1 , IOMSG = msg)
    
    nsdx_temp = nsdx
    
    if(status1 == 0 ) then                                           
 
        read(10,NML=PARA)
        read(10,NML=PARA_debug)
        
        ! 测试收敛阶时，避免从输入文件中读取相同的单元划分数
        if(compuPurpose == order_solve .or. compuPurpose == order_solve2)then            
            nsdx = nsdx_temp
            nsdy = nsdx
        end if
    else                                                                
        ! 返回错误信息
        write(*,"('Error opening file : IOSTAT = ',I6)") status1
        write(*,*) trim(msg)         
        stop  
    end if
    
    
    close(10)     
end subroutine read_parameter

subroutine pre_setting 

    !-----------------------------------------------------------------------------
    !
    !   描述：根据输入文件，预处理提前求解程序所用常量
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use parameter_setting   
    use global_var
    use type_module
    implicit none
 
    
    call nsp_setting            !----确定空间离散阶数
    
    call Gauss_points_set       !----Gauss点
    
    call Lobatto_points_set     !----Lobatto点
    
    call SPs_setting            !----解点坐标赋值
       
    call FPs_setting            !----通量点坐标赋值
    
    call get_collocation        !----求解矫正函数在解点的值
    
    call get_dis_sp_fp          !----求解子单元解点与通量点的距离
    
    call get_Indicator_para     !----Modal Energy侦测所需,Lagrange->Legend转化矩阵

end subroutine pre_setting 
    
subroutine grid_generate

    !-----------------------------------------------------------------------------
    !
    !   描述：网格生成或网格读取
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use parameter_setting 
    implicit none
    
    write(*,*) 'Read meshes.'
    select case(grid_set)!
    case(self)
        call self_unstru_input
    case(fluent_cas)
        call fluent_unstru_input
    case default
        write(*,"('Error! No Defined grid generate type')")
        stop
    end select

end subroutine grid_generate
    
subroutine fluid_initialize

    !-----------------------------------------------------------------------------
    !
    !   描述：流场初始化
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    implicit none
    
    write(*,*) 'Fluid initialization.'
    
    call set_BC_fluent              !----设置边界条件(根据cell所在的位置)
    
    call pre_BoundaryCell_fluent    !----设置虚拟边界单元和求解点顶点坐标     
    
    call set_IC                     !----设置初始条件

end subroutine fluid_initialize

subroutine get_min_dis 

    !-----------------------------------------------------------------------------
    !
    !   方法：求解最小距离
    !   描述：求解网格单元的顶点的x方向最小距离，顶点的y方向最小距离，中心点的最小距离。根据CFL数确定dt，目前使用的是min_dis
    !   返回：mindx_global, mindy_global, min_dis
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use global_var
    use type_module
    implicit none
    
    integer :: i,j,nearCellIndex
    real(prec) :: dx_temp,dy_temp,xc_x0,xc_x1,yc_y0,yc_y1,dis
    mindx_global = 1.0e5
    mindy_global = 1.0e5
    min_dis = 1.0e5
    
    ! 顶点的x方向最小距离，顶点的y方向最小距离
    
    do i = 1,ncells
        dx_temp = abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(i).nodes(2),1))
        mindx_global = min(mindx_global,dx_temp)
        dx_temp = abs(xy_coor(cellset(i).nodes(3),1)-xy_coor(cellset(i).nodes(4),1))
        mindx_global = min(mindx_global,dx_temp)
        
        dy_temp = abs(xy_coor(cellset(i).nodes(1),2)-xy_coor(cellset(i).nodes(4),2))
        mindy_global = min(mindy_global,dy_temp)
        dy_temp = abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(i).nodes(3),2))
        mindy_global = min(mindy_global,dy_temp)  
    end do
    
    ! 单元中心点的最小距离
    
    do i = 1,ncells
        xc_x0 = xy_coor(cellset(i).nodes(1),1)+xy_coor(cellset(i).nodes(3),1)
        yc_y0 = xy_coor(cellset(i).nodes(1),2)+xy_coor(cellset(i).nodes(3),2)
        
        do j = 1,4
            nearCellIndex = cellset(i).nearcells(j)
            xc_x1 = xy_coor(cellset(nearCellIndex).nodes(1),1)+xy_coor(cellset(nearCellIndex).nodes(3),1)
            yc_y1 = xy_coor(cellset(nearCellIndex).nodes(1),2)+xy_coor(cellset(nearCellIndex).nodes(3),2)
            
            dis = sqrt((xc_x0 - xc_x1)**2+(yc_y0 - yc_y1)**2)
            min_dis = min(min_dis,dis)
        end do  
    end do
    min_dis = min_dis*0.5_prec

end subroutine get_min_dis 
    
    
subroutine nsp_setting

    !-----------------------------------------------------------------------------
    !
    !   方法：设置求解点数目
    !   描述：根据CPR阶数，确定解点分布。选择Gauss点分布，Larange多项式，单边的解点数等于阶数
    !   返回：nsp
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use parameter_setting
    implicit none
    
    nsp = nsp_set
    !select case(nsp_set)!
    !case(nsp_one)
    !    nsp = nsp_one
    !case(nsp_two)
    !    nsp = nsp_two
    !case(nsp_t)
    !    nsp = nsp_t
    !case(nsp_fo)
    !    nsp = nsp_fo
    !case(nsp_fi)
    !    nsp = nsp_fi
    !case default
    !    write(*,*)'Error! No defined order input, default: nsp = 3!'
    !    stop
    !end select

end subroutine nsp_setting

subroutine SPs_setting     

    !-----------------------------------------------------------------------------
    !
    !   方法：设置求解点类型
    !   描述：内部解点采取Gauss Points， Lobatto求解点 ，等距点.            
    !         等距求解点： 把单元均分为nsp个的子单元中心点. 涉及到积分的时候就不能直接用这些点值直接计算
    !   返回：SPs(:), SPs_local(:,:)
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use parameter_setting
    use global_var
    implicit none
    integer :: i,j
    real(prec) :: dh
    
    ! 1D 标准单元求解点坐标
    
    allocate(SPs(nsp)) 
    select case(sp_type)
    case(Gauss_p)           
        SPs = GPs                   
    case(Lobatto_p)
        ! TODO
        SPs=LPs 
    case(Equidistant_p)    
        dh = 2.0_prec / real(nsp)
        do i = 1,nsp
            SPs(i) = -1.0_prec + dh*(i-1) + 0.5_prec*dh
        end do                       
    case default
        write(*,*)'Error! Occur in subroutine solution_points_set'
    end select
        
    ! 2D 标准单元求解点坐标
        
    allocate(SPs_local(nsp,nsp,2))
    do j = 1,nsp
        SPs_local(j,:,1) = SPs(:)   ! x
        SPs_local(:,j,2) = SPs(:)   ! y
    end do
end subroutine SPs_setting

subroutine FPs_setting    

    !-----------------------------------------------------------------------------
    !
    !   方法：设置通量点类型
    !   描述：内部通量点采取Lobatto Points， Gauss Weight Points求解点 ，等距点. 通量点比求解点多一个            
    !   返回：nsp
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use parameter_setting
    use global_var
    implicit none
    integer j
    allocate(FPs(nsp+1))
    select case(fp_type)
    case(Lobatto_FPs)
        FPs = LPs 
    case(GaussW_FPs)        
        FPs(1) = -1.0_prec
        if(nsp>1)then
            do j = 2,nsp
                FPs(j) = FPs(j-1)+GCoe(j-1)
            end do           
        end if
        FPs(nsp+1) = 1.0_prec
    case(Equidistant_FPs)   
        FPs(1) = -1.0_prec
        if(nsp>1)then
            do j = 2,nsp
                FPs(j) = FPs(j-1) + 2.0_prec/nsp
            end do
        end if
        FPs(nsp+1) = 1.0_prec
    case default
        write(*,*)'Error! Occur in subroutine solution_points_set'
    end select
   
end subroutine FPs_setting

subroutine Gauss_points_set

    !-----------------------------------------------------------------------------
    !
    !   方法：设置Gauss点和Gauss权
    !   描述：          
    !   返回：
    !   作者：gqShi 
    !   历史：修改求解Gauss点和Gauss权方法，可以生成任意点的点值和权值（当然实际上也有限制，和程序实现有关，参见程序get_Gauss_Point_Weight）//2021.12.15
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use Gauss_P_W
    implicit none
    real(prec_dd) :: GPs_tmp(nsp),GCoe_tmp(nsp)     !----四精度，保证计算准确
    
    allocate(GPs(nsp),GCoe(nsp))
    
    Gps = 0.0_prec
    GCoe = 0.0_prec
    
    call get_Gauss_Point_Weight(GPs_tmp,GCoe_tmp)   !----根据点数任意生成Gauss点值和高斯权值
    
    GPs(:) = GPs_tmp(:)
    GCoe(:) = GCoe_tmp(:)

    !旧的求解方法
    !select case(nsp)
    !case(nsp_one)
    !    
    !    GPs(1) = 0.0_prec
    !
    !    GCoe(1) = 2.0_prec
    !    
    !case(nsp_two)
    !    
    !    GPs(1) = -0.5773502691896257645092_prec
    !    GPs(2) = 0.5773502691896257645092_prec
    !     
    !    GCoe(1) = 1.0_prec
    !    GCoe(2) = 1.0_prec
    !    
    !case(nsp_t)
    !    
    !    GPs(1) = -sqrt(0.6_prec)
    !    GPs(2) = 0.0_prec
    !    GPs(3) = sqrt(0.6_prec)
    !     
    !    GCoe(1) = 5.0_prec/9.0_prec
    !    GCoe(2) = 8.0_prec/9.0_prec
    !    GCoe(3) = 5.0_prec/9.0_prec
    !    
    !case(nsp_fo)
    !    
    !    GPs(1) = -0.861136311594052575224_prec
    !    GPs(2) = -0.3399810435848562648027_prec
    !    GPs(3) = 0.3399810435848562648027_prec
    !    GPs(4) = 0.861136311594052575224_prec
    !    
    !    GCoe(1) = 0.3478548451374538573731_prec
    !    GCoe(2) = 0.6521451548625461426269_prec
    !    GCoe(3) = 0.6521451548625461426269_prec
    !    GCoe(4) = 0.3478548451374538573731_prec
    !
    !case(nsp_fi)
    !     
    !    GPs(1) = -0.9061798459386639927976_prec
    !    GPs(2) = -0.5384693101056830910363_prec
    !    GPs(3) =  0.0_prec
    !    GPs(4) =  0.5384693101056830910363_prec
    !    GPs(5) =  0.9061798459386639927976_prec
    !    
    !    GCoe(1) = 0.2369268850561890875143_prec
    !    GCoe(2) = 0.4786286704993664680413_prec
    !    GCoe(3) = 0.5688888888888888888889_prec
    !    GCoe(4) = 0.4786286704993664680413_prec
    !    GCoe(5) = 0.2369268850561890875143_prec
    !           
    !case default
    !    write(*,"('Error! No Defined Gauss Points Number,nsp ='I3 )")nsp
    !    stop
    !end select
    
end subroutine Gauss_points_set

subroutine Lobatto_points_set

    !-----------------------------------------------------------------------------
    !
    !   方法：设置Lobatto点坐标
    !   描述：          
    !   返回：
    !   作者：gqShi 
    !   历史：修改求解Lobatto点，可以生成任意点的点值和权值（当然实际上也有限制，和程序实现有关，参见程序get_Lobatto_Point_Weight）//2021.12.15
    !   这里面不知道为什么Lobatto点的个数比预先设定多1个，先接受这个设定
    !   观察了斐然师兄的程序里面Lobatto点个数和解点数一致的，先进行修改
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use Lobatto_P_W
    implicit none
    
    real(prec_dd) :: LPs_tmp(nsp+1),LCoe_tmp(nsp+1)     !----四精度，保证计算准确
    
    allocate(LPs(nsp+1))
    
    call get_Lobatto_Point_Weight(LPs_tmp,LCoe_tmp) !----根据点数任意生成Lobatto点值和权值
    
    LPs(:) = LPs_tmp(:)        

    
    !
    !select case(nsp)
    !case(nsp_one)
    !         
    !    LPs(1) = -1.000000000000000000000_prec 
    !    LPs(2) =  1.000000000000000000000_prec 
    !        
    !case(nsp_two)
    !         
    !    LPs(1) = -1.000000000000000000000_prec
    !    LPs(2) =  0.0_prec
    !    LPs(3) =  1.000000000000000000000_prec
    !case(nsp_t)
    !         
    !    LPs(1) = -1.000000000000000000000_prec 
    !    LPs(2) = -0.447213595499957939282_prec
    !    LPs(3) =  0.447213595499957939282_prec
    !    LPs(4) =  1.000000000000000000000_prec 
    !
    !case(nsp_fo)         
    !         
    !    LPs(1) = -1.0000000000000000000000_prec 
    !    LPs(2) = -0.6546536707079771437983_prec 
    !    LPs(3) =  0.0000000000000000000000_prec
    !    LPs(4) =  0.6546536707079771437983_prec
    !    LPs(5) =  1.0000000000000000000000_prec
    !        
    !case(nsp_fi)
    !         
    !    LPs(1) = -1.000000000000000000000_prec 
    !    LPs(2) = -0.765055323929464692851_prec
    !    LPs(3) = -0.2852315164806450963142_prec
    !    LPs(4) =  0.2852315164806450963142_prec
    !    LPs(5) =  0.765055323929464692851_prec            
    !    LPs(6) =  1.000000000000000000000_prec
    !        
    !case default
    !    write(*,"('Error! No Defined Gauss Points Number,nsp ='I3 )")nsp
    !    stop
    !end select

end subroutine Lobatto_points_set
    
subroutine get_dis_sp_fp

    !-----------------------------------------------------------------------------
    !
    !   方法：求解点和通量点之间的距离（计算空间）
    !   描述：          
    !   返回：
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    implicit none
    
    integer :: i,j
    
    allocate(dis_sp_fp(2*nsp))
    do i = 1,nsp
        dis_sp_fp(2*i-1) = SPs(i) - FPs(i)
        dis_sp_fp(2*i)   = FPs(i+1) - SPs(i)
    end do

end subroutine get_dis_sp_fp

subroutine get_Indicator_para

    !-----------------------------------------------------------------------------
    !
    !   方法：预设侦测所需变量值
    !   描述：MDH侦测所需,Lagrange->Legend转化矩阵. K_La_Le  原始形式，K_La_Le2 添加边界点的修改形式
    !         目前适用于5点CPR格式，需要拓展研究任意点的实现程序
    !   返回：
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use parameter_setting   
    use global_var
    use type_module
    
    implicit none
    integer :: i
    
    allocate(K_La_Le(nsp,nsp),K_La_Le_2(nsp+2,nsp+2))
    K_La_Le=reshape( (/  -0.44542134632044532872_prec,  0.44856404312571887327_prec,    1.4079281687625479597_prec,     0.44856404312571887327_prec,    -0.44542134632044532872_prec,&
                         -0.26295072534780090115_prec,  -0.31564963775813707471_prec,   0.0_prec,                       0.31564963775813707471_prec,    0.26295072534780090115_prec, &
                         0.27412134137104513456_prec,   -0.04924826331462704873_prec,   -0.44974615611283617166_prec,   -0.04924826331462704873_prec,   0.27412134137104513456_prec, &
                         -0.222081873567294817185_prec, 0.37373739635316142603_prec,    0.0_prec,                       -0.37373739635316142603_prec,   0.222081873567294817185_prec,&
                         0.123506106334137037799_prec,  -0.349780276313832245607_prec,  0.452548339959390415617_prec,   -0.349780276313832245607_prec,  0.123506106334137037799_prec/),   &
                        shape(K_La_Le), order=(/2,1/) ) 
    K_La_Le_2=reshape( (/  0.0_prec, 0.0_prec, 0.0_prec, 1.00000000000000000000_prec, 0.0_prec, 0.0_prec, 0.0_prec,&
 -0.93750000000000000000_prec, 1.68402696032017526627_prec, -2.02152893587538279044_prec,    0.0_prec,                       2.02152893587538279044_prec, -1.68402696032017526627_prec,   0.93750000000000000000_prec,&
 0.9375000000000000000_prec,  -1.8583805056663784005_prec,   3.7542138389997117338_prec,     -5.6666666666666666667_prec,  3.7542138389997117338_prec, -1.8583805056663784005_prec,   0.9375000000000000000_prec,&
 4.3750000000000000000_prec,  -7.4920339228080194586_prec,   4.4833198487460769383_prec,     0.0_prec,                      -4.4833198487460769383_prec,   7.4920339228080194586_prec, -4.3750000000000000000_prec,&
 -4.3750000000000000000_prec,  8.2677119297962491841_prec, -8.3260452631295825175_prec,      8.8666666666666666667_prec,  -8.3260452631295825175_prec,   8.2677119297962491841_prec, -4.3750000000000000000_prec,&
 -3.9375000000000000000_prec,  5.8080069624878441923_prec, -2.4617909128706941478_prec,      0.0_prec,                       2.4617909128706941478_prec, -5.8080069624878441923_prec,   3.9375000000000000000_prec,&
 3.93750000000000000000_prec, -6.4093314241298707836_prec,   4.57183142412987078361_prec,    -4.20000000000000000000_prec, 4.57183142412987078361_prec, -6.4093314241298707836_prec,   3.93750000000000000000_prec/),&
                        shape(K_La_Le_2), order=(/2,1/) )                   

end subroutine get_Indicator_para 

subroutine get_collocation

    !-----------------------------------------------------------------------------
    !
    !   方法：修正系数
    !   描述：修正函数在固定求解点是定值，预先求出，节省计算量. gDG修正函数          
    !   返回：
    !   作者：gqShi 
    !   历史：修改求解修正项方法，可以生成任意阶的惩罚系数（当然实际上也有限制，和程序实现有关，参见程序get_gDG_collocation）//2021.12.15
    !-----------------------------------------------------------------------------

    use real_precision
    use global_var
    use parameter_setting
    use Correction_Func
    implicit none
    
    integer j
    real(prec),external :: gl_sub_kesi,gr_sub_kesi 
    real(prec_dd):: gl_coe_tmp(nsp),gr_coe_tmp(nsp) 
    
    allocate(gl_coe(nsp),gr_coe(nsp))
 
    call get_gDG_collocation(gl_coe_tmp,gr_coe_tmp)
    
    gl_coe = gl_coe_tmp
    gr_coe = gr_coe_tmp


    ! 旧的求解方法，直接推导出公式，然后编写子程序求解
    !do j = 1,nsp       
    !    gl_coe(j) = gl_sub_kesi(SPs(j))
    !    gr_coe(j) = gr_sub_kesi(SPs(j))
    !end do
end subroutine get_collocation


