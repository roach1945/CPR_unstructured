module real_precision
    
    !-----------------------------------------------------------------------------
    !
    !   模块：实型精度
    !   描述：预设实型变量精度值。包括单精度 (kind=4 / kind=1.0e0)、双精度(kind=8 / kind=1.0d0)
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
    
    implicit none
    integer,parameter :: single_prec = 4
    integer,parameter :: double_prec = 8 
    integer,parameter :: ddouble_prec = 16

    integer,parameter :: prec = double_prec
    integer,parameter :: prec_dd = ddouble_prec
end module real_precision
    
module preParameter
    
    !-----------------------------------------------------------------------------
    !
    !   模块：程序预设功能
    !   描述：参数调节求解目的，输出相应结果
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    implicit none
    
    !! 求解目的 测误差，测精度 
    
    integer :: compuPurpose = 0                                                      
    integer,parameter :: error_solve = 0,order_solve = 1,order_solve2 = 2            
    
    !! 求解目的 测离散守恒  
    
    integer :: compuConsLaw = 0 
    integer,parameter :: consLaw_close = 0,consLaw_open = 1
    
    !! 求解目的 测残差 
    
    integer :: compuResidua = 0                                             
    integer,parameter :: residua_no = 0,residua_yes = 1
    
    !! 求解目的 记录问题单元 
    
    integer :: trouCell_record = 0                                              
    integer,parameter :: trouCell_no = 0,trouCell_yes = 1

    Namelist /prePARA/ compuPurpose,compuConsLaw,compuResidua,trouCell_record 
    
end module preParameter
    
module parameter_setting

    !-----------------------------------------------------------------------------
    !
    !   模块：程序预设参数
    !   描述：参数调节不同方法实现，实现不同功能
    !         预设初始值，在此模块中不做更改
    !         可在pre_processing 读取 input.in 覆盖
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
  
    use real_precision
    implicit none
    
    !----预设值-------------------------------------------------------------------       
    
    real(prec) :: xl = -5.0_prec, xr = 5.0_prec     !! 定义求解区间：左值 xl，右值 xr                                      
    real(prec) :: yl = -5.0_prec, yr = 5.0_prec
    integer :: nsdx = 40,nsdy = 40                  !! Number of space discretize 空间离散数                                                             
    integer :: nsp = 3                              !! Number of solution points 解点数                     
    integer :: nstep = 30000                        !! number of time step 时间步           
    real(prec) :: dt = 0.0001_prec                  !! 时间步长
    real(prec) :: cfl = 0.001_prec                  !! CFL数    
    real(prec) :: T = 1.0_prec                      !! 计算时间   
    real(prec) :: TVB_M = 0.01_prec                 !! TVB侦测器自由参数 M
    real(prec) :: efix = 0.1_prec                   !! Roe Flux 自由参数
    real(prec) :: d_smooth = 4.0_prec               !! 格式过渡，控制MDH侦测器侦测区域参数
    !------------------------------------------------------------------------------
    
    
    !----其他全局变量-------------------------------------------------------------- 
    
    character(len=32) :: mesh_file      ! 网格文件名
    integer :: print_step               ! 输出步数
    real(prec) :: Ms,Mv,attackAngle     ! 算例中shock 强度，vortex 强度，攻角
    !------------------------------------------------------------------------------
    
    
    !----参数控制器定义------------------------------------------------------------
    
    ! 求解点个数. 3 -- 3点；4 -- 4点；5 -- 5点
    integer :: nsp_set = 3                  
    integer,parameter :: nsp_one = 1,nsp_two = 2,nsp_t = 3, nsp_fo = 4, nsp_fi = 5  
    
    ! 通量导数计算方法. Larange Polynomial, Chain Rule
    integer :: fluxD_type = 0
    integer,parameter :: LP = 0,CR = 1  
    
    ! 网格生成方法. self, fluent cas file
    integer :: grid_set = 0         
    integer,parameter :: self = 0, fluent_cas = 1

    ! 自生成网格类型. 直网格, 曲网格(TODO) 
    integer :: grid_type = 0    
    integer,parameter :: grid_straight = 0, grid_curve = 1
    
    ! 直网格类型. 矩形, 梯形, 扰动
    integer :: cell_shape = 0   
    integer,parameter :: cell_rect = 0, cell_trap = 1, cell_dist = 2           
    
    ! 求解点类型. Gauss point, Lobatto point, equidistan point    
    integer :: sp_type = 0      
    integer,parameter :: Gauss_p = 0, Lobatto_p = 1,Equidistant_p = 2
    
    ! 通量点类型. Labotto, Gauss weight, equidistan point    
    integer :: fp_type = 0      
    integer,parameter :: Lobatto_FPs = 0, GaussW_FPs = 1, Equidistant_FPs = 2
    
    ! 插值变量. 原始变量, 特征变量
    integer :: var_type
    integer,parameter :: ori_type = 0, character_type = 1  
    
    ! 标准算例选择. 等熵涡, Sod, Lax, Shu-Osher, 2D Riemann, double Mach, shock-vortex interaction, 
    !               Shock-composite vortex interaction, cylinder, wing, steady shovk,2D sin wave
    integer :: case_comp = 0            
    integer,parameter :: equEntropy_case = 0, SodShockTube_case = 1, LaxShockTube_case = 2,ShuOsher_case = 3,Riemann2D_case = 4,&
                         DoubleMach_case = 5, VortexShock_case  = 6, CompositeVortexShock_case = 7,HyperCylinder_case = 8,test_case = 99,&
                         WingFlow_case   = 9, steadyShock_case = 10, SinWave_2D_case = 11, Doublerarefaction_1D_case=12
    
    ! 激波管问题激波运动方向. x, y
    integer :: dire_shock = 0
    integer,parameter :: direct_x = 0, direct_y = 1
    
    ! 问题单元侦测. TVB, ATV, MV, MDH, MDHm, KXRCF, JST
    integer :: detection_type
    integer,parameter :: detection_TVB = 0,detection_ATV = 1,detection_MV = 2,detection_MDH = 3,detection_MDHm = 4, detection_KXRCF = 5, detection_JST = 6
    
    ! 时间步长的选取方法. CFL 自适应, 固定dt
    integer :: dt_switch = 1
    integer,parameter :: dt_adapt = 0, dt_fixed = 1
    
    ! 求解器选择. 1D, 2D
    integer :: solver_dim = 2
    integer,parameter :: dim_1D = 1, dim_2D = 2
    
    
    Namelist /PARA/ xl,xr,yl,yr,nsdx,nsdy,fluxD_type,T,nstep,dt,cfl,nsp_set,grid_set,grid_type,cell_shape,&
                    sp_type,fp_type,case_comp,var_type,TVB_M,dire_shock,mesh_file,detection_type,print_step,&
                    Ms,Mv,dt_switch,solver_dim,efix,attackAngle,d_smooth
    !------------------------------------------------------------------------------
    
    
    !----在此定义调试程序所需变量--------------------------------------------------
    
    ! 通量计算左右值插值方法. NNW采取非线性插值，Godunov左右值等于解点处值，为一阶方法
    integer :: method_subcell = 0                           
    integer,parameter :: method_NNW = 0, method_Godunov = 1 
    
    ! 针对单元边界通量点处初始值，NNW采取非线性插值的不同方法。// 两点反距离加权，多点反距离加权,两单元两点线性插值再平均
    integer :: method_ori_fp = 0                            
    integer,parameter :: TPs_weight_NNW = 0, MPs_weight_NNW = 1, TTPsL_ave_NNW = 2  
    
    ! NNW求对流项导数方法. 两点差分，三点加权处理
    integer :: method_adv_Der = 0                           
    integer,parameter :: Operator_2Ps = 0, Operator_3Ps = 1
    
    ! 侦测类型. 分维，分单元
    integer :: detect_type = 1                              
    integer,parameter :: ByDim = 0, ByCell = 1          
    
    ! 是否加过渡单元. 
    integer :: buffer_cell_switch = 0                       
    integer,parameter :: buffer_cell_no = 0, buffer_cell_yes = 1    
    
    ! 格式类型. CPR, CNNW2, hybrid CPR-CNNW2 
    integer :: scheme_kind = 2                       
    integer,parameter :: scheme_cpr = 0, scheme_two = 1,scheme_hybrid = 2 
    
    ! 格式过渡. 双阈值调节
    integer :: switch_updown = 0                       
    integer,parameter :: close_updown = 0, open_updown = 1
    
    ! NNW插值 限制器开关
    integer :: limiter_switch = 1                       
    integer,parameter :: limiter_close = 0, limiter_open = 1
    
    ! 通量导数的两种求法. FPs = SPs, FPs = other  
    integer :: fluxDer_switch = 0                       
    integer,parameter :: fluxDer_SPs = 0, fluxDer_FPs = 1
    
    ! Riemann Flux. LLF, Roe  
    integer :: Riemann_flux = 0                     
    integer,parameter :: LLF_flux = 0, Roe_flux = 1
    
    
    Namelist /PARA_debug/ method_subcell,method_ori_fp,method_adv_Der,detect_type,buffer_cell_switch,scheme_kind,limiter_switch,fluxDer_switch,Riemann_flux,switch_updown
    !--------------------------------------------------------------------------------
    
end module parameter_setting

module global_var

    !其他全局变量设置，利用输入参数可控制或本身不需要改变的

    use real_precision
    use parameter_setting
    implicit none
   
    real(prec),parameter :: pi = 3.1415926535897932384626433832795_prec
    real(prec),parameter :: gamma = 1.4_prec,gamma1=gamma-1.0_prec
    
    real(prec) :: xlong,ylong
    real(prec) :: dt_temp,T_temp
    real(prec),parameter :: kesi_l = -1.0_prec, kesi_r = 1.0_prec           !标准单元 ：左值 -1.0，右值 1.0
    
    real(prec),dimension(:),allocatable :: GPs,GCoe,SPs,LPs,FPs             !Gauss points,Gauss cofficient
    real(prec),dimension(:),allocatable :: dis_sp_fp                        !子单元解点与通量点之间的距离，C2NNW2需要(2nsp)
    real(prec),dimension(:,:,:),allocatable :: SPs_local                    !Solution points 标准单元坐标，具体多少视CPR阶数定(nsp,nsp,2) 2指x,y坐标
    real(prec),dimension(:),allocatable :: gl_coe,gr_coe
    integer,parameter :: nghn = 1                                           !Number of Left Boundary of Mould，模板向外扩展数
    
    integer :: ni,nj,nsdpx,nsdpy                                            !nsdpx,nsdpy是x方向，y方向离散点，nsdx,nsdy是x,y方向离散单元
    integer :: nnodes,nsides,ncells,nbdsides 
    integer,dimension(:),allocatable  :: BoundCells_index_set
    integer :: n_BCells_D = 0,n_BCells_R = 0,n_BCells_U = 0,n_BCells_L = 0  !计算区域边界上的侧边数
    real(prec),dimension(:,:),allocatable :: xx,yy                          !顶点的横、纵坐标
    real(prec),dimension(:,:),allocatable :: xy_coor                        !X of solution points，单元顶点的全局坐标，对应顶点编号
    integer :: sum_vertexCell                                               !后处理。包含该顶点的子单元总数，填补单元间空白
    real(prec),dimension(10) :: a_Norm_L1,a_Norm_L2,a_Norm_Linf             !array of Norm Ll,array of Norm Lint
    real(prec) :: Norm_L1,Norm_L2,Norm_Linf                                 !LI范数误差，L2范数误差，Linf范数误差
    integer :: counter = 1                                                  !计数器
    integer :: nt_temp,hour_whole,minite_whole,second_whole                 !计时所用变量
    
    real(prec),dimension(:, :),allocatable :: K_La_Le,K_La_Le_2             !MDH侦测，转换矩阵
    
    real(prec) :: sum_q,sum_q_initial,sum_a_exact                           !离散守恒律测试所用变量
    real(prec) :: maxResiduals,FirstResi,FirstAveResi                       !残差测试
    real(prec),dimension(:) :: CL_xy(2)                                     !翼型绕流，升力系数
    real(prec) :: mindx_global,mindy_global,min_dis                         !全局横纵坐标变化值，全局单元中心最小距离，求dt   
    
    !contains
    !subroutine Allocate_memory_Cell
    !    implicit none
    !    !---分配数组内存-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !    !---结束-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !end subroutine Allocate_memory_Cell 

end module global_var

module type_module

    !-----------------------------------------------------------------------------
    !
    !   方法：定义流场数据储存
    !   描述：单元结构体，边结构体
    !   作者：gqShi 
    !   历史：** 
    !         2021.12.08 注释
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use global_var,only:nnodes,nsides,ncells,nbdsides !bound condition cells
    implicit none
    
    ! 单元结构体类型
    
    type cell_record                    
        
        integer    :: index                             !单元编号
        integer    :: nodes(4)                          !4顶点编号
        integer    :: sides(4)                          !4条边编号
        integer    :: nearcells(4)                      !4个相邻单元编号
        character(len = 16)    :: bc                    !见fluent 边界条件数字代码 Fluent_File_Doc.txt
        real(prec),allocatable :: det_J(:,:)            !Jacobi 行列式值  
        real(prec),allocatable :: MJacobi(:,:,:)        !Matrix Jacobi (nsp,nsp,4)    逆度量    
        real(prec),allocatable :: Mdirect(:,:,:)        !直接度量矩阵 (nsp,nsp,4)  
        real(prec),allocatable :: fpdet_J_F(:,:)        !flux point Jacobi(nsp,nsp+1)按行数
        real(prec),allocatable :: fpdet_J_G(:,:)        !flux point Jacobi(nsp,nsp+1)按列数
        real(prec),allocatable :: fpMdirect_F(:,:,:)    !flux point直接度量矩阵 (nsp,nsp+1,4) nsp行
        real(prec),allocatable :: fpMdirect_G(:,:,:)    !flux point直接度量矩阵 (nsp,nsp+1,4) nsp列
        
        real(prec),allocatable :: sp_coor(:,:,:)        !内部解点全局坐标    (nsp,nsp,2)
        real(prec),allocatable :: spvalue_ori(:,:,:)    !内部解点(nsp,nsp,4)  原始变量( r u v p )
        real(prec),allocatable :: spvalue_con(:,:,:)    !内部解点(nsp,nsp,4)  守恒变量( r ru rv E)
        real(prec),allocatable :: spvalue_fluF(:,:,:)   !内部解点(nsp,nsp,4)  kesi方向通量  
        real(prec),allocatable :: spvalue_fluG(:,:,:)   !内部解点(nsp,nsp,4)  eat 方向通量 
        
        real(prec),allocatable :: spvalue_con_loc(:,:,:)    !计算空间内部解点(nsp,nsp,4)  守恒变量( r ru rv E)
        real(prec),allocatable :: spvalue_fluF_loc(:,:,:)   !计算空间内部解点(nsp,nsp,4)  kesi方向通量  
        real(prec),allocatable :: spvalue_fluG_loc(:,:,:)   !计算空间内部解点(nsp,nsp,4)  eat 方向通量 
        
        real(prec),allocatable :: spvalue_con_tem(:,:,:)    !RK推进暂时储存数组，替代spvalue_con(:,:,:)
        
        real(prec),allocatable :: spvalue_ori_exa(:,:,:)    !准确解 内部解点(nsp,nsp,4)  原始变量( r u v p )
        real(prec),allocatable :: spvalue_ori_old(:,:,:)
        
        ! 关于子单元限制的一些变量
        
        integer :: Beta,Beta_old = 0                                        !标定单元是否Trouble Cell Beta = 1为Trouble Cell 
        integer,allocatable  :: Beta_line(:,:)                              !(nsp,2)分维侦测,(i,1)代表横方向是否问题，(i,2)代表纵方向是否问题 
        real(prec),allocatable :: Smooth_line_x(:,:),Smooth_line_y(:,:)     !(nsp,nsp,2),(:,:,1)代表横方向
        real(prec),allocatable :: fluxF_innerfp(:,:,:),fluxG_innerfp(:,:,:) !储存子单元通量，也即单元内部通量点
        
    end type cell_record
    
    ! 边结构体类型
    
    type side_record                    
  
        integer :: nodes(2)                             !起点，终点
        integer :: nearcells(2)                         !左单元，右单元
        character(len = 16)    :: bc                    !见fluent 边界条件数字代码 Fluent_File_Doc.txt
        real(prec),allocatable :: fp_coor(:,:,:)        !边界通量点全局坐标   (4,nsp,2)
        real(prec),allocatable :: FPdirect(:,:,:)       !直接度量矩阵         (nsp,nsp,4) 
        real(prec),allocatable :: fpvalue_upw(:,:)      !边界迎风通量         (nsp,4)
        
    end type side_record
    
    ! 定义单元结构体和边结构体 
    type(cell_record),allocatable :: cellset(:) 
    type(side_record),allocatable :: sideset(:)
    
    contains
    
    ! 分配内存
    subroutine allo_stru_record
        allocate(cellset(ncells+nbdsides))  
        allocate(sideset(nsides))
    end subroutine allo_stru_record
    
end module type_module
    
module bc_module

    !-----------------------------------------------------------------------------
    !
    !   方法：流场边界类型模块
    !   描述：定义流场边界类型
    !   作者：gqShi 
    !   历史：** 
    !         2021.12.08 添加注释
    !
    !-----------------------------------------------------------------------------

    implicit none
    
    ! 新的边界类型可以在这添加
    
    character(len = 16) :: Interior = 'Interior',Wall = '3',Inlet_vent = '4',Outlet_vent = '5',Periodic = '20'!,Periodic_X = '20X',Periodic_Y = '20Y'
    
    
end module bc_module

module Time_module

    !-----------------------------------------------------------------------------
    !
    !   方法：计算时间模块
    !   描述：定义时间统计所需变量
    !   作者：gqShi 
    !   历史：** 
    !         2021.12.08 添加注释
    !
    !-----------------------------------------------------------------------------

    use real_precision
    implicit none
    real(prec) :: TimeStart,TimeEnd,TimeFrame1,TimeFrame2,TimeFrame3
    real(prec) :: DetectionTime, DetectionTimeStart, DetectionTimeEnd
    
end module Time_module
    
module time_test
    use real_precision
    implicit none
    real(prec) :: Ts,Te,T1,T2,T3,T4,T5
   
end module time_test

module test_module

    !中间测试用的模块，程序写完删，留着也行。。

    type test1
        integer,allocatable :: array(:)
    end type test1
    integer :: m,n
    type(test1),allocatable :: arrayset(:)
    contains
    subroutine test_sub1
        allocate(arrayset(m))
    end subroutine test_sub1

end module test_module
