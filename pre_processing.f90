subroutine pre_processing

    write(*,*)'------------------------------------------------------'  
    write(*,*)'------------------------------------------------------'
    
    call read_parameter         !----��ȡ���Ʋ���

    call pre_setting            !----���ݶ�ȡ�Ĳ��������������� 

    call grid_generate          !----��������   
    
    call fluid_initialize       !----��ʼ������ 
    
    call get_min_dis            !----����Сȫ�ֺ�������任ֵ����dt
        
    call print_submesh_data     !----����

    call print_mesh_data        !----���plt��ʽ����
    
    open(unit=1010,file='error.dat')
    
end subroutine pre_processing 

subroutine read_parameter

    !-----------------------------------------------------------------------------
    !
    !   ��������ȡԤ�����. �����ļ� input.in 
    !   ���ߣ�gqShi 
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
        
        ! ����������ʱ������������ļ��ж�ȡ��ͬ�ĵ�Ԫ������
        if(compuPurpose == order_solve .or. compuPurpose == order_solve2)then            
            nsdx = nsdx_temp
            nsdy = nsdx
        end if
    else                                                                
        ! ���ش�����Ϣ
        write(*,"('Error opening file : IOSTAT = ',I6)") status1
        write(*,*) trim(msg)         
        stop  
    end if
    
    
    close(10)     
end subroutine read_parameter

subroutine pre_setting 

    !-----------------------------------------------------------------------------
    !
    !   ���������������ļ���Ԥ������ǰ���������ó���
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------

    use real_precision
    use parameter_setting   
    use global_var
    use type_module
    implicit none
 
    
    call nsp_setting            !----ȷ���ռ���ɢ����
    
    call Gauss_points_set       !----Gauss��
    
    call Lobatto_points_set     !----Lobatto��
    
    call SPs_setting            !----������긳ֵ
       
    call FPs_setting            !----ͨ�������긳ֵ
    
    call get_collocation        !----�����������ڽ���ֵ
    
    call get_dis_sp_fp          !----����ӵ�Ԫ�����ͨ����ľ���
    
    call get_Indicator_para     !----Modal Energy�������,Lagrange->Legendת������

end subroutine pre_setting 
    
subroutine grid_generate

    !-----------------------------------------------------------------------------
    !
    !   �������������ɻ������ȡ
    !   ���ߣ�gqShi 
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
    !   ������������ʼ��
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------

    implicit none
    
    write(*,*) 'Fluid initialization.'
    
    call set_BC_fluent              !----���ñ߽�����(����cell���ڵ�λ��)
    
    call pre_BoundaryCell_fluent    !----��������߽絥Ԫ�����㶥������     
    
    call set_IC                     !----���ó�ʼ����

end subroutine fluid_initialize

subroutine get_min_dis 

    !-----------------------------------------------------------------------------
    !
    !   �����������С����
    !   �������������Ԫ�Ķ����x������С���룬�����y������С���룬���ĵ����С���롣����CFL��ȷ��dt��Ŀǰʹ�õ���min_dis
    !   ���أ�mindx_global, mindy_global, min_dis
    !   ���ߣ�gqShi 
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
    
    ! �����x������С���룬�����y������С����
    
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
    
    ! ��Ԫ���ĵ����С����
    
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
    !   ����������������Ŀ
    !   ����������CPR������ȷ�����ֲ���ѡ��Gauss��ֲ���Larange����ʽ�����ߵĽ�������ڽ���
    !   ���أ�nsp
    !   ���ߣ�gqShi 
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
    !   ������������������
    !   �������ڲ�����ȡGauss Points�� Lobatto���� ���Ⱦ��.            
    !         �Ⱦ����㣺 �ѵ�Ԫ����Ϊnsp�����ӵ�Ԫ���ĵ�. �漰�����ֵ�ʱ��Ͳ���ֱ������Щ��ֱֵ�Ӽ���
    !   ���أ�SPs(:), SPs_local(:,:)
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------

    use parameter_setting
    use global_var
    implicit none
    integer :: i,j
    real(prec) :: dh
    
    ! 1D ��׼��Ԫ��������
    
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
        
    ! 2D ��׼��Ԫ��������
        
    allocate(SPs_local(nsp,nsp,2))
    do j = 1,nsp
        SPs_local(j,:,1) = SPs(:)   ! x
        SPs_local(:,j,2) = SPs(:)   ! y
    end do
end subroutine SPs_setting

subroutine FPs_setting    

    !-----------------------------------------------------------------------------
    !
    !   ����������ͨ��������
    !   �������ڲ�ͨ�����ȡLobatto Points�� Gauss Weight Points���� ���Ⱦ��. ͨ����������һ��            
    !   ���أ�nsp
    !   ���ߣ�gqShi 
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
    !   ����������Gauss���GaussȨ
    !   ������          
    !   ���أ�
    !   ���ߣ�gqShi 
    !   ��ʷ���޸����Gauss���GaussȨ�������������������ĵ�ֵ��Ȩֵ����Ȼʵ����Ҳ�����ƣ��ͳ���ʵ���йأ��μ�����get_Gauss_Point_Weight��//2021.12.15
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use Gauss_P_W
    implicit none
    real(prec_dd) :: GPs_tmp(nsp),GCoe_tmp(nsp)     !----�ľ��ȣ���֤����׼ȷ
    
    allocate(GPs(nsp),GCoe(nsp))
    
    Gps = 0.0_prec
    GCoe = 0.0_prec
    
    call get_Gauss_Point_Weight(GPs_tmp,GCoe_tmp)   !----���ݵ�����������Gauss��ֵ�͸�˹Ȩֵ
    
    GPs(:) = GPs_tmp(:)
    GCoe(:) = GCoe_tmp(:)

    !�ɵ���ⷽ��
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
    !   ����������Lobatto������
    !   ������          
    !   ���أ�
    !   ���ߣ�gqShi 
    !   ��ʷ���޸����Lobatto�㣬�������������ĵ�ֵ��Ȩֵ����Ȼʵ����Ҳ�����ƣ��ͳ���ʵ���йأ��μ�����get_Lobatto_Point_Weight��//2021.12.15
    !   �����治֪��ΪʲôLobatto��ĸ�����Ԥ���趨��1�����Ƚ�������趨
    !   �۲����Ȼʦ�ֵĳ�������Lobatto������ͽ����һ�µģ��Ƚ����޸�
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use Lobatto_P_W
    implicit none
    
    real(prec_dd) :: LPs_tmp(nsp+1),LCoe_tmp(nsp+1)     !----�ľ��ȣ���֤����׼ȷ
    
    allocate(LPs(nsp+1))
    
    call get_Lobatto_Point_Weight(LPs_tmp,LCoe_tmp) !----���ݵ�����������Lobatto��ֵ��Ȩֵ
    
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
    !   �����������ͨ����֮��ľ��루����ռ䣩
    !   ������          
    !   ���أ�
    !   ���ߣ�gqShi 
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
    !   ������Ԥ������������ֵ
    !   ������MDH�������,Lagrange->Legendת������. K_La_Le  ԭʼ��ʽ��K_La_Le2 ��ӱ߽����޸���ʽ
    !         Ŀǰ������5��CPR��ʽ����Ҫ��չ�о�������ʵ�ֳ���
    !   ���أ�
    !   ���ߣ�gqShi 
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
    !   ����������ϵ��
    !   ���������������ڹ̶������Ƕ�ֵ��Ԥ���������ʡ������. gDG��������          
    !   ���أ�
    !   ���ߣ�gqShi 
    !   ��ʷ���޸�������������������������׵ĳͷ�ϵ������Ȼʵ����Ҳ�����ƣ��ͳ���ʵ���йأ��μ�����get_gDG_collocation��//2021.12.15
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


    ! �ɵ���ⷽ����ֱ���Ƶ�����ʽ��Ȼ���д�ӳ������
    !do j = 1,nsp       
    !    gl_coe(j) = gl_sub_kesi(SPs(j))
    !    gr_coe(j) = gr_sub_kesi(SPs(j))
    !end do
end subroutine get_collocation


