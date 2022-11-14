module real_precision
    
    !-----------------------------------------------------------------------------
    !
    !   ģ�飺ʵ�;���
    !   ������Ԥ��ʵ�ͱ�������ֵ������������ (kind=4 / kind=1.0e0)��˫����(kind=8 / kind=1.0d0)
    !   ���ߣ�gqShi 
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
    !   ģ�飺����Ԥ�蹦��
    !   �����������������Ŀ�ģ������Ӧ���
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    implicit none
    
    !! ���Ŀ�� �����⾫�� 
    
    integer :: compuPurpose = 0                                                      
    integer,parameter :: error_solve = 0,order_solve = 1,order_solve2 = 2            
    
    !! ���Ŀ�� ����ɢ�غ�  
    
    integer :: compuConsLaw = 0 
    integer,parameter :: consLaw_close = 0,consLaw_open = 1
    
    !! ���Ŀ�� ��в� 
    
    integer :: compuResidua = 0                                             
    integer,parameter :: residua_no = 0,residua_yes = 1
    
    !! ���Ŀ�� ��¼���ⵥԪ 
    
    integer :: trouCell_record = 0                                              
    integer,parameter :: trouCell_no = 0,trouCell_yes = 1

    Namelist /prePARA/ compuPurpose,compuConsLaw,compuResidua,trouCell_record 
    
end module preParameter
    
module parameter_setting

    !-----------------------------------------------------------------------------
    !
    !   ģ�飺����Ԥ�����
    !   �������������ڲ�ͬ����ʵ�֣�ʵ�ֲ�ͬ����
    !         Ԥ���ʼֵ���ڴ�ģ���в�������
    !         ����pre_processing ��ȡ input.in ����
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
  
    use real_precision
    implicit none
    
    !----Ԥ��ֵ-------------------------------------------------------------------       
    
    real(prec) :: xl = -5.0_prec, xr = 5.0_prec     !! ����������䣺��ֵ xl����ֵ xr                                      
    real(prec) :: yl = -5.0_prec, yr = 5.0_prec
    integer :: nsdx = 40,nsdy = 40                  !! Number of space discretize �ռ���ɢ��                                                             
    integer :: nsp = 3                              !! Number of solution points �����                     
    integer :: nstep = 30000                        !! number of time step ʱ�䲽           
    real(prec) :: dt = 0.0001_prec                  !! ʱ�䲽��
    real(prec) :: cfl = 0.001_prec                  !! CFL��    
    real(prec) :: T = 1.0_prec                      !! ����ʱ��   
    real(prec) :: TVB_M = 0.01_prec                 !! TVB��������ɲ��� M
    real(prec) :: efix = 0.1_prec                   !! Roe Flux ���ɲ���
    real(prec) :: d_smooth = 4.0_prec               !! ��ʽ���ɣ�����MDH���������������
    !------------------------------------------------------------------------------
    
    
    !----����ȫ�ֱ���-------------------------------------------------------------- 
    
    character(len=32) :: mesh_file      ! �����ļ���
    integer :: print_step               ! �������
    real(prec) :: Ms,Mv,attackAngle     ! ������shock ǿ�ȣ�vortex ǿ�ȣ�����
    !------------------------------------------------------------------------------
    
    
    !----��������������------------------------------------------------------------
    
    ! �������. 3 -- 3�㣻4 -- 4�㣻5 -- 5��
    integer :: nsp_set = 3                  
    integer,parameter :: nsp_one = 1,nsp_two = 2,nsp_t = 3, nsp_fo = 4, nsp_fi = 5  
    
    ! ͨ���������㷽��. Larange Polynomial, Chain Rule
    integer :: fluxD_type = 0
    integer,parameter :: LP = 0,CR = 1  
    
    ! �������ɷ���. self, fluent cas file
    integer :: grid_set = 0         
    integer,parameter :: self = 0, fluent_cas = 1

    ! ��������������. ֱ����, ������(TODO) 
    integer :: grid_type = 0    
    integer,parameter :: grid_straight = 0, grid_curve = 1
    
    ! ֱ��������. ����, ����, �Ŷ�
    integer :: cell_shape = 0   
    integer,parameter :: cell_rect = 0, cell_trap = 1, cell_dist = 2           
    
    ! ��������. Gauss point, Lobatto point, equidistan point    
    integer :: sp_type = 0      
    integer,parameter :: Gauss_p = 0, Lobatto_p = 1,Equidistant_p = 2
    
    ! ͨ��������. Labotto, Gauss weight, equidistan point    
    integer :: fp_type = 0      
    integer,parameter :: Lobatto_FPs = 0, GaussW_FPs = 1, Equidistant_FPs = 2
    
    ! ��ֵ����. ԭʼ����, ��������
    integer :: var_type
    integer,parameter :: ori_type = 0, character_type = 1  
    
    ! ��׼����ѡ��. ������, Sod, Lax, Shu-Osher, 2D Riemann, double Mach, shock-vortex interaction, 
    !               Shock-composite vortex interaction, cylinder, wing, steady shovk,2D sin wave
    integer :: case_comp = 0            
    integer,parameter :: equEntropy_case = 0, SodShockTube_case = 1, LaxShockTube_case = 2,ShuOsher_case = 3,Riemann2D_case = 4,&
                         DoubleMach_case = 5, VortexShock_case  = 6, CompositeVortexShock_case = 7,HyperCylinder_case = 8,test_case = 99,&
                         WingFlow_case   = 9, steadyShock_case = 10, SinWave_2D_case = 11, Doublerarefaction_1D_case=12
    
    ! ���������⼤���˶�����. x, y
    integer :: dire_shock = 0
    integer,parameter :: direct_x = 0, direct_y = 1
    
    ! ���ⵥԪ���. TVB, ATV, MV, MDH, MDHm, KXRCF, JST
    integer :: detection_type
    integer,parameter :: detection_TVB = 0,detection_ATV = 1,detection_MV = 2,detection_MDH = 3,detection_MDHm = 4, detection_KXRCF = 5, detection_JST = 6
    
    ! ʱ�䲽����ѡȡ����. CFL ����Ӧ, �̶�dt
    integer :: dt_switch = 1
    integer,parameter :: dt_adapt = 0, dt_fixed = 1
    
    ! �����ѡ��. 1D, 2D
    integer :: solver_dim = 2
    integer,parameter :: dim_1D = 1, dim_2D = 2
    
    
    Namelist /PARA/ xl,xr,yl,yr,nsdx,nsdy,fluxD_type,T,nstep,dt,cfl,nsp_set,grid_set,grid_type,cell_shape,&
                    sp_type,fp_type,case_comp,var_type,TVB_M,dire_shock,mesh_file,detection_type,print_step,&
                    Ms,Mv,dt_switch,solver_dim,efix,attackAngle,d_smooth
    !------------------------------------------------------------------------------
    
    
    !----�ڴ˶�����Գ����������--------------------------------------------------
    
    ! ͨ����������ֵ��ֵ����. NNW��ȡ�����Բ�ֵ��Godunov����ֵ���ڽ�㴦ֵ��Ϊһ�׷���
    integer :: method_subcell = 0                           
    integer,parameter :: method_NNW = 0, method_Godunov = 1 
    
    ! ��Ե�Ԫ�߽�ͨ���㴦��ʼֵ��NNW��ȡ�����Բ�ֵ�Ĳ�ͬ������// ���㷴�����Ȩ����㷴�����Ȩ,����Ԫ�������Բ�ֵ��ƽ��
    integer :: method_ori_fp = 0                            
    integer,parameter :: TPs_weight_NNW = 0, MPs_weight_NNW = 1, TTPsL_ave_NNW = 2  
    
    ! NNW������������. �����֣������Ȩ����
    integer :: method_adv_Der = 0                           
    integer,parameter :: Operator_2Ps = 0, Operator_3Ps = 1
    
    ! �������. ��ά���ֵ�Ԫ
    integer :: detect_type = 1                              
    integer,parameter :: ByDim = 0, ByCell = 1          
    
    ! �Ƿ�ӹ��ɵ�Ԫ. 
    integer :: buffer_cell_switch = 0                       
    integer,parameter :: buffer_cell_no = 0, buffer_cell_yes = 1    
    
    ! ��ʽ����. CPR, CNNW2, hybrid CPR-CNNW2 
    integer :: scheme_kind = 2                       
    integer,parameter :: scheme_cpr = 0, scheme_two = 1,scheme_hybrid = 2 
    
    ! ��ʽ����. ˫��ֵ����
    integer :: switch_updown = 0                       
    integer,parameter :: close_updown = 0, open_updown = 1
    
    ! NNW��ֵ ����������
    integer :: limiter_switch = 1                       
    integer,parameter :: limiter_close = 0, limiter_open = 1
    
    ! ͨ��������������. FPs = SPs, FPs = other  
    integer :: fluxDer_switch = 0                       
    integer,parameter :: fluxDer_SPs = 0, fluxDer_FPs = 1
    
    ! Riemann Flux. LLF, Roe  
    integer :: Riemann_flux = 0                     
    integer,parameter :: LLF_flux = 0, Roe_flux = 1
    
    
    Namelist /PARA_debug/ method_subcell,method_ori_fp,method_adv_Der,detect_type,buffer_cell_switch,scheme_kind,limiter_switch,fluxDer_switch,Riemann_flux,switch_updown
    !--------------------------------------------------------------------------------
    
end module parameter_setting

module global_var

    !����ȫ�ֱ������ã�������������ɿ��ƻ�����Ҫ�ı��

    use real_precision
    use parameter_setting
    implicit none
   
    real(prec),parameter :: pi = 3.1415926535897932384626433832795_prec
    real(prec),parameter :: gamma = 1.4_prec,gamma1=gamma-1.0_prec
    
    real(prec) :: xlong,ylong
    real(prec) :: dt_temp,T_temp
    real(prec),parameter :: kesi_l = -1.0_prec, kesi_r = 1.0_prec           !��׼��Ԫ ����ֵ -1.0����ֵ 1.0
    
    real(prec),dimension(:),allocatable :: GPs,GCoe,SPs,LPs,FPs             !Gauss points,Gauss cofficient
    real(prec),dimension(:),allocatable :: dis_sp_fp                        !�ӵ�Ԫ�����ͨ����֮��ľ��룬C2NNW2��Ҫ(2nsp)
    real(prec),dimension(:,:,:),allocatable :: SPs_local                    !Solution points ��׼��Ԫ���꣬���������CPR������(nsp,nsp,2) 2ָx,y����
    real(prec),dimension(:),allocatable :: gl_coe,gr_coe
    integer,parameter :: nghn = 1                                           !Number of Left Boundary of Mould��ģ��������չ��
    
    integer :: ni,nj,nsdpx,nsdpy                                            !nsdpx,nsdpy��x����y������ɢ�㣬nsdx,nsdy��x,y������ɢ��Ԫ
    integer :: nnodes,nsides,ncells,nbdsides 
    integer,dimension(:),allocatable  :: BoundCells_index_set
    integer :: n_BCells_D = 0,n_BCells_R = 0,n_BCells_U = 0,n_BCells_L = 0  !��������߽��ϵĲ����
    real(prec),dimension(:,:),allocatable :: xx,yy                          !����ĺᡢ������
    real(prec),dimension(:,:),allocatable :: xy_coor                        !X of solution points����Ԫ�����ȫ�����꣬��Ӧ������
    integer :: sum_vertexCell                                               !���������ö�����ӵ�Ԫ���������Ԫ��հ�
    real(prec),dimension(10) :: a_Norm_L1,a_Norm_L2,a_Norm_Linf             !array of Norm Ll,array of Norm Lint
    real(prec) :: Norm_L1,Norm_L2,Norm_Linf                                 !LI������L2������Linf�������
    integer :: counter = 1                                                  !������
    integer :: nt_temp,hour_whole,minite_whole,second_whole                 !��ʱ���ñ���
    
    real(prec),dimension(:, :),allocatable :: K_La_Le,K_La_Le_2             !MDH��⣬ת������
    
    real(prec) :: sum_q,sum_q_initial,sum_a_exact                           !��ɢ�غ��ɲ������ñ���
    real(prec) :: maxResiduals,FirstResi,FirstAveResi                       !�в����
    real(prec),dimension(:) :: CL_xy(2)                                     !��������������ϵ��
    real(prec) :: mindx_global,mindy_global,min_dis                         !ȫ�ֺ�������仯ֵ��ȫ�ֵ�Ԫ������С���룬��dt   
    
    !contains
    !subroutine Allocate_memory_Cell
    !    implicit none
    !    !---���������ڴ�-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !    !---����-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !end subroutine Allocate_memory_Cell 

end module global_var

module type_module

    !-----------------------------------------------------------------------------
    !
    !   �����������������ݴ���
    !   ��������Ԫ�ṹ�壬�߽ṹ��
    !   ���ߣ�gqShi 
    !   ��ʷ��** 
    !         2021.12.08 ע��
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use global_var,only:nnodes,nsides,ncells,nbdsides !bound condition cells
    implicit none
    
    ! ��Ԫ�ṹ������
    
    type cell_record                    
        
        integer    :: index                             !��Ԫ���
        integer    :: nodes(4)                          !4������
        integer    :: sides(4)                          !4���߱��
        integer    :: nearcells(4)                      !4�����ڵ�Ԫ���
        character(len = 16)    :: bc                    !��fluent �߽��������ִ��� Fluent_File_Doc.txt
        real(prec),allocatable :: det_J(:,:)            !Jacobi ����ʽֵ  
        real(prec),allocatable :: MJacobi(:,:,:)        !Matrix Jacobi (nsp,nsp,4)    �����    
        real(prec),allocatable :: Mdirect(:,:,:)        !ֱ�Ӷ������� (nsp,nsp,4)  
        real(prec),allocatable :: fpdet_J_F(:,:)        !flux point Jacobi(nsp,nsp+1)������
        real(prec),allocatable :: fpdet_J_G(:,:)        !flux point Jacobi(nsp,nsp+1)������
        real(prec),allocatable :: fpMdirect_F(:,:,:)    !flux pointֱ�Ӷ������� (nsp,nsp+1,4) nsp��
        real(prec),allocatable :: fpMdirect_G(:,:,:)    !flux pointֱ�Ӷ������� (nsp,nsp+1,4) nsp��
        
        real(prec),allocatable :: sp_coor(:,:,:)        !�ڲ����ȫ������    (nsp,nsp,2)
        real(prec),allocatable :: spvalue_ori(:,:,:)    !�ڲ����(nsp,nsp,4)  ԭʼ����( r u v p )
        real(prec),allocatable :: spvalue_con(:,:,:)    !�ڲ����(nsp,nsp,4)  �غ����( r ru rv E)
        real(prec),allocatable :: spvalue_fluF(:,:,:)   !�ڲ����(nsp,nsp,4)  kesi����ͨ��  
        real(prec),allocatable :: spvalue_fluG(:,:,:)   !�ڲ����(nsp,nsp,4)  eat ����ͨ�� 
        
        real(prec),allocatable :: spvalue_con_loc(:,:,:)    !����ռ��ڲ����(nsp,nsp,4)  �غ����( r ru rv E)
        real(prec),allocatable :: spvalue_fluF_loc(:,:,:)   !����ռ��ڲ����(nsp,nsp,4)  kesi����ͨ��  
        real(prec),allocatable :: spvalue_fluG_loc(:,:,:)   !����ռ��ڲ����(nsp,nsp,4)  eat ����ͨ�� 
        
        real(prec),allocatable :: spvalue_con_tem(:,:,:)    !RK�ƽ���ʱ�������飬���spvalue_con(:,:,:)
        
        real(prec),allocatable :: spvalue_ori_exa(:,:,:)    !׼ȷ�� �ڲ����(nsp,nsp,4)  ԭʼ����( r u v p )
        real(prec),allocatable :: spvalue_ori_old(:,:,:)
        
        ! �����ӵ�Ԫ���Ƶ�һЩ����
        
        integer :: Beta,Beta_old = 0                                        !�궨��Ԫ�Ƿ�Trouble Cell Beta = 1ΪTrouble Cell 
        integer,allocatable  :: Beta_line(:,:)                              !(nsp,2)��ά���,(i,1)����᷽���Ƿ����⣬(i,2)�����ݷ����Ƿ����� 
        real(prec),allocatable :: Smooth_line_x(:,:),Smooth_line_y(:,:)     !(nsp,nsp,2),(:,:,1)����᷽��
        real(prec),allocatable :: fluxF_innerfp(:,:,:),fluxG_innerfp(:,:,:) !�����ӵ�Ԫͨ����Ҳ����Ԫ�ڲ�ͨ����
        
    end type cell_record
    
    ! �߽ṹ������
    
    type side_record                    
  
        integer :: nodes(2)                             !��㣬�յ�
        integer :: nearcells(2)                         !��Ԫ���ҵ�Ԫ
        character(len = 16)    :: bc                    !��fluent �߽��������ִ��� Fluent_File_Doc.txt
        real(prec),allocatable :: fp_coor(:,:,:)        !�߽�ͨ����ȫ������   (4,nsp,2)
        real(prec),allocatable :: FPdirect(:,:,:)       !ֱ�Ӷ�������         (nsp,nsp,4) 
        real(prec),allocatable :: fpvalue_upw(:,:)      !�߽�ӭ��ͨ��         (nsp,4)
        
    end type side_record
    
    ! ���嵥Ԫ�ṹ��ͱ߽ṹ�� 
    type(cell_record),allocatable :: cellset(:) 
    type(side_record),allocatable :: sideset(:)
    
    contains
    
    ! �����ڴ�
    subroutine allo_stru_record
        allocate(cellset(ncells+nbdsides))  
        allocate(sideset(nsides))
    end subroutine allo_stru_record
    
end module type_module
    
module bc_module

    !-----------------------------------------------------------------------------
    !
    !   �����������߽�����ģ��
    !   ���������������߽�����
    !   ���ߣ�gqShi 
    !   ��ʷ��** 
    !         2021.12.08 ���ע��
    !
    !-----------------------------------------------------------------------------

    implicit none
    
    ! �µı߽����Ϳ����������
    
    character(len = 16) :: Interior = 'Interior',Wall = '3',Inlet_vent = '4',Outlet_vent = '5',Periodic = '20'!,Periodic_X = '20X',Periodic_Y = '20Y'
    
    
end module bc_module

module Time_module

    !-----------------------------------------------------------------------------
    !
    !   ����������ʱ��ģ��
    !   ����������ʱ��ͳ���������
    !   ���ߣ�gqShi 
    !   ��ʷ��** 
    !         2021.12.08 ���ע��
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

    !�м�����õ�ģ�飬����д��ɾ������Ҳ�С���

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
