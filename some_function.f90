!���ļ��ռ���̹����ߴ�·�Ĵ��룬�����������Ӱ�졣
!subroutine Perspective_Mappings
!    !��ע������͸�ӱ任����ʵ��ͶӰ�Ŀ��ǣ��������ڵ�Ԫ�߽��ϵĽ�㲻һ�£���ʽ���غ㣬�ʲ�ѡ�ô��ӳ���ķ���
!    !���ÿ����Ԫ�任����׼��Ԫ��Jacobi����
!    !�Ƶ���ʽ��ֱ�Ӽ���
!    !���崦����: 1.�����A(0,0)(1,0)(1,1)(0,1)���ǽṹ�ı���B�ı任���� B = AM1;                            --���ղο��������
!    !              2.��M1�����,Ҳ��B��A�ı任���� BM1_inv = A;                                               --��ֵ����
!    !              3.C = AM2,C(-1,-1)(1,-1)(1,1)(-1,1),M2=[2 0 0; 0 2 0; -1 -1 1];                            --���ݲο������ܽ�
!    !              4.�õ�C = BM1_invM2 = BM3;                                                                       --���Դ���
!    !              5.��� J                                                                                   --����
!    !              ����Ϊʲô��ֱ���Ƶ�B->C,����A��������ʹ��M1�ܼ򵥣����������ĵ�Ԫ�����̻���鷳����һ��������Ҫ�á�
!    !Refs
!    !   [1] ͸���任���� https://www.cnblogs.com/purehol/p/11734085.html
!    !   [2]��ͼ����͸�ӱ任 Perspective Transformation https://blog.csdn.net/xiaowei_cqu/article/details/26471527
!     
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer i,j,k
!    real(prec) :: dx1,dx2,dx3,dy1,dy2,dy3
!    real(prec),dimension(:,:) :: M1(3,3),M1_inv(3,3),M2(3,3),M3(3,3),M3_inv(3,3),C(3,3)
!    real(prec),dimension(:) :: V_sp_local(3)
!    
!
!    do i = 1,ncells
!        dx1 = xy_coor(cellset(i).nodes(2),1) - xy_coor(cellset(i).nodes(3),1)
!        dx2 = xy_coor(cellset(i).nodes(4),1) - xy_coor(cellset(i).nodes(3),1)
!        dx3 = xy_coor(cellset(i).nodes(1),1) - xy_coor(cellset(i).nodes(2),1) + xy_coor(cellset(i).nodes(3),1) - xy_coor(cellset(i).nodes(4),1)
!        dy1 = xy_coor(cellset(i).nodes(2),2) - xy_coor(cellset(i).nodes(3),2)
!        dy2 = xy_coor(cellset(i).nodes(4),2) - xy_coor(cellset(i).nodes(3),2)
!        dy3 = xy_coor(cellset(i).nodes(1),2) - xy_coor(cellset(i).nodes(2),2) + xy_coor(cellset(i).nodes(3),2) - xy_coor(cellset(i).nodes(4),2)
!        !write(*,*)'ddd',dx1,dx2,dx3,dy1,dy2,dy3
!        if(dx3 == 0 .and. dy3 == 0)then
!            M1(1,1) = xy_coor(cellset(i).nodes(2),1) - xy_coor(cellset(i).nodes(1),1)
!            M1(2,1) = xy_coor(cellset(i).nodes(3),1) - xy_coor(cellset(i).nodes(2),1)
!            M1(3,1) = xy_coor(cellset(i).nodes(1),1)
!            M1(1,2) = xy_coor(cellset(i).nodes(2),2) - xy_coor(cellset(i).nodes(1),2)
!            M1(2,2) = xy_coor(cellset(i).nodes(3),2) - xy_coor(cellset(i).nodes(2),2)
!            M1(3,2) = xy_coor(cellset(i).nodes(1),2)
!            M1(1,3) = 0.0_prec
!            M1(2,3) = 0.0_prec
!            M1(3,3) = 1.0_prec  
!            !write(*,*)'------+++++'
!        else
!            M1(1,3) = (dx3*dy2-dx2*dy3)/(dx1*dy2-dx2*dy1)
!            M1(2,3) = (dx1*dy3-dx3*dy1)/(dx1*dy2-dx2*dy1)
!            M1(3,3) = 1.0_prec
!            M1(1,1) = xy_coor(cellset(i).nodes(2),1)*(1.0_prec + M1(1,3)) - xy_coor(cellset(i).nodes(1),1)
!            M1(1,2) = xy_coor(cellset(i).nodes(2),2)*(1.0_prec + M1(1,3)) - xy_coor(cellset(i).nodes(1),2)
!            M1(2,1) = xy_coor(cellset(i).nodes(4),1)*(1.0_prec + M1(2,3)) - xy_coor(cellset(i).nodes(1),1)
!            M1(2,2) = xy_coor(cellset(i).nodes(4),2)*(1.0_prec + M1(2,3)) - xy_coor(cellset(i).nodes(1),2)
!            M1(3,1) = xy_coor(cellset(i).nodes(1),1)
!            M1(3,2) = xy_coor(cellset(i).nodes(1),2)
!        end if
!
!        !��M1�棬A = BM1_inv;
!        call Matrix_Inverse(M1,3,M1_inv)                !B = AM1; A = BM1_inv
!
!        !��M3, C = BM3
!        M2 = RESHAPE((/2.0_prec,0.0_prec,-1.0_prec,0.0_prec,2.0_prec,-1.0_prec,0.0_prec,0.0_prec,1.0_prec/),(/3,3/))   !C = AM2
!        M3 = Matmul(M1_inv,M2)                          !C = AM2 = BM1_invM2 = BM3
!
!        !��¼M3_inv��cellset���ñ�׼��Ԫ�Ľ������ֵ����ʵ�������н�������
!        call Matrix_Inverse(M3,3,M3_inv)                !B =CM3_inv
!        cellset(i).MC2P(:,:) = M3_inv(:,:)    !computational to physical
!        
!        !----����M3_inv��⵼��ƫx/ƫkesi�ȣ���Jacobi����---------------------------------------------------------            
!        !�����ڴ��MJacobi
!        allocate(cellset(i).MJacobi(nsp,nsp,4))
!        
!        V_sp_local = 1.0_prec
!        do j = 1,nsp
!            do k = 1,nsp
!                !xkesi
!                V_sp_local(1:2) = SPs_local(j,k,1:2)
!                cellset(i).MJacobi(j,k,1) = (M3_inv(1,1)*(dot_product(V_sp_local(:),M3_inv(:,3)))+&
!                                           M3_inv(1,3)*(dot_product(V_sp_local(:),M3_inv(:,1))))/&
!                                          (dot_product(V_sp_local(:),M3_inv(:,3)))**2
!                !xeta
!                cellset(i).MJacobi(j,k,2) = (M3_inv(2,1)*(dot_product(V_sp_local(:),M3_inv(:,3)))+&
!                                           M3_inv(2,3)*(dot_product(V_sp_local(:),M3_inv(:,1))))/&
!                                          (dot_product(V_sp_local(:),M3_inv(:,3)))**2 
!                !ykesi                         
!                cellset(i).MJacobi(j,k,3) = (M3_inv(1,2)*(dot_product(V_sp_local(:),M3_inv(:,3)))+&
!                                           M3_inv(1,3)*(dot_product(V_sp_local(:),M3_inv(:,2))))/&
!                                          (dot_product(V_sp_local(:),M3_inv(:,3)))**2  
!                !yeta                          
!                cellset(i).MJacobi(j,k,4) = (M3_inv(2,2)*(dot_product(V_sp_local(:),M3_inv(:,3)))+&
!                                           M3_inv(2,3)*(dot_product(V_sp_local(:),M3_inv(:,2))))/&
!                                          (dot_product(V_sp_local(:),M3_inv(:,3)))**2 
!                !write(*,*) cellset(i).MJacobi(j,k,:)                          
!            end do
!        end do
!       
!        !write(*,*) 'cccccccc'
!        !do j = 1,4
!        !    write(*,*) xy_coor(cellset(1).nodes(j),:)
!        !end do
!        !do j =1,3
!        !    write(*,*) M3(j,:)
!        !end do
!    end do 
!end subroutine Perspective_Mappings
!
!subroutine solve_SP_coor
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer i,j,k
!    real(prec),dimension(:,:) :: M_sp_local(3,3),Mtemp(3,3)
!    do i = 1,ncells
!        allocate(cellset(i).sp_coor(nsp,nsp,2))!�ڲ����ȫ������    (nsp,nsp)
!    end do
!      
!    M_sp_local(2:3,:) = 0.0_prec
!    do i = 1,ncells
!        do j = 1,nsp
!            do k =1,nsp
!                M_sp_local(1,1:2) = SPs_local(j,k,1:2)
!                M_sp_local(1,3) = 1.0_prec
!                Mtemp = Matmul(M_sp_local,cellset(i).MC2P)/dot_product(M_sp_local(1,:),cellset(i).MC2P(:,3))
!                cellset(i).sp_coor(j,k,1:2) = Mtemp(1,1:2)          
!                !write(*,*) Mtemp(1,1:2)
!            end do
!        end do        
!    end do
!!!    do j = 1,nsp
!!!        do k =1,nsp
!!!            M_sp_local(1,1) = -1.0_prec
!!!            M_sp_local(1,2) = 1.0_prec
!!!            M_sp_local(1,3) = 1.0_prec
!!!            Mtemp = Matmul(M_sp_local,cellset(1).MC2P)      
!!!            !write(*,*) Mtemp(1,1:2)
!!!        end do
!!!    end do 
!    
!end subroutine solve_SP_coor
!module RHS_module
!    use real_precision
!    use parameter_setting,only:nsp
!    implicit none
!    real(prec),dimension(:,:,:),allocatable :: FG_upw
!    real(prec),dimension(:,:,:),allocatable :: O1,O2,Q1,Q2,F1,F2,G1,G2
!    real(prec),dimension(:,:),allocatable  :: detJ
!end module RHS_module    
!subroutine get_local_fpvalue(index)
!    !(���ͨ��һ��Ԫһ���£����ô�������ӭ��ͨ�����Դ������������ظ�����)
!    use parameter_setting,only:nsp
!     
!    use type_module
!    use RHS_module
!    use global_var
!    implicit none
!    integer index,i,j,k,l,m,n
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(1:2) :: norm
!    real(prec),dimension(1:4) :: u1,u2,fg1,fg2,ori1,ori2
!    integer :: sideIndex,nearCellIndex
!    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4)
!    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2
!    real(prec),dimension(:,:) :: quad_vertex(4,2) 
!    !F G
!    do i =1,4 !side
!        sideIndex = cellset(index).sides(i)
!        !if(sideset(sideindex).fpvalue_upw(1,1) /= 1000.0_prec) cycle    !ÿ��ʱ�䲽���ӭ��ͨ����ʼ��1000������Ϊ1000����ǰ�沽�����󣬲����ظ�����
!        nearCellIndex = cellset(index).nearcells(i)
!        !�������±߽���ͨ��
!        !1��ʾ��Ԫ�ڣ�2��ʾ�ٵ�Ԫ
!        Q1 = cellset(index).spvalue_con_loc
!        F1 = cellset(index).spvalue_fluF_loc 
!        G1 = cellset(index).spvalue_fluG_loc
!        Q2 = cellset(nearCellIndex).spvalue_con_loc
!        F2 = cellset(nearCellIndex).spvalue_fluF_loc 
!        G2 = cellset(nearCellIndex).spvalue_fluG_loc
!        !write(*,*) index,sideIndex,nearCellIndex
!        do l = 1,3 
!            do m =1,4
!                if(i == 1)then                  
!                    u1(m) = LaI_nPs(kesi_l,SPs,Q1(:,l,m),nsp)
!                    u2(m) = LaI_nPs(kesi_r,SPs,Q2(:,l,m),nsp)
!                    fg1(m) = LaI_nPs(kesi_l,SPs,G1(:,l,m),nsp)!��
!                    fg2(m) = LaI_nPs(kesi_r,SPs,G2(:,l,m),nsp)!��
!                    !write(*,*)'g',G1(i,l,m)!,G2,fg1,fg2
!                else if(i == 2 )then
!                    u1(m) = LaI_nPs(kesi_r,SPs,Q1(l,:,m),nsp)
!                    u2(m) = LaI_nPs(kesi_l,SPs,Q2(l,:,m),nsp)
!                    fg1(m) = LaI_nPs(kesi_r,SPs,F1(l,:,m),nsp)!��
!                    fg2(m) = LaI_nPs(kesi_l,SPs,F2(l,:,m),nsp)!��
!                    !write(*,*)'f',fg1(m),fg2(m)
!                else if(i == 3)then                  
!                    u1(m) = LaI_nPs(kesi_r,SPs,Q1(:,l,m),nsp)
!                    u2(m) = LaI_nPs(kesi_l,SPs,Q2(:,l,m),nsp)
!                    fg1(m) = LaI_nPs(kesi_r,SPs,G1(:,l,m),nsp)!��
!                    fg2(m) = LaI_nPs(kesi_l,SPs,G2(:,l,m),nsp)!��
!                else if(i == 4 )then
!                    u1(m) = LaI_nPs(kesi_l,SPs,Q1(l,:,m),nsp)
!                    u2(m) = LaI_nPs(kesi_r,SPs,Q2(l,:,m),nsp)
!                    fg1(m) = LaI_nPs(kesi_l,SPs,F1(l,:,m),nsp)!��
!                    fg2(m) = LaI_nPs(kesi_r,SPs,F2(l,:,m),nsp)!��
!                end if    
!            end do
!            if(i == 1.or.i == 4)then 
!                n = -1
!            else
!                n = 1
!            end if 
!            
!            call Lax_F_flux(u1,u2,fg1,fg2,n,FG_upw(i,l,:))  
!            !write(*,*)'upw',FG_upw(2,l,1)
!        end do            
!        sideset(sideindex).fpvalue_upw(:,:) = FG_upw(i,:,:)!i-th side, l-th fp ,m-th var
!        !write(*,*)sideIndex,sideset(sideindex).fpvalue_upw(:,:)
!        !write(*,*)FG_upw
!    end do
!  
!    
!    
!    
!end subroutine get_local_fpvalue
    
    !subroutine Lax_F_flux(ul,ur,ql,qr,n,q_upw_temp)!array of Un/Un+1
!    !һά
!    use real_precision
!    use global_var,only:gamma
!    implicit none
!    real(prec),dimension(1:4) :: ul,ur,ql,qr,q_upw_temp,ul_conser,ur_conser
!    real(prec) :: p,r,c,cl,cr,lambda
!    integer n
!    ul_conser(1) = ul(1)                        !rho
!    ul_conser(2) = ul(1)*ul(2)   !rho*u
!    ul_conser(3) = ul(1)*ul(3)   !rho*v
!    ul_conser(4) = ul(4)/(gamma-1)+0.5_prec*ul(1)*(ul(2)**2+ul(3)**2)!E
!    
!    ur_conser(1) = ur(1)                        !rho
!    ur_conser(2) = ur(1)*ur(2)   !rho*u
!    ur_conser(3) = ur(1)*ur(3)   !rho*v
!    ur_conser(4) = ur(4)/(gamma-1)+0.5_prec*ur(1)*(ur(2)**2+ur(3)**2)!E
!    
!    p = ul(4)
!    r = ul(1)
!    c = gamma*p/r
!    cl = sqrt(c)!����
!
!    p = ur(4)
!    r = ur(1)
!    c = gamma*p/r
!    cr = sqrt(c)
!    lambda = max(abs(ul(2)+cl),abs(ur(2)+cr))
!    q_upw_temp = 0.5_prec*(ql + qr - lambda*n*(ur_conser-ul_conser))
!
!end subroutine Lax_F_flux
    
    !
!subroutine solve_error_order
!    use global_var
!    implicit none
!    integer i,j
!    real(prec),dimension(10) :: order_L1,order_L2,order_Linf
!    character(len=30) file_name
!    character(len=5) char_T,char_nsd
!!!!!!!!!!��ⷶ����� 
!    write(*,*) " "
!    write(*,*) "   L1 error       ", "   L2 error       ","Linf error"
!    do i = 1,counter-1
!        write(*,*) a_Norm_L1(i),a_Norm_L2(i),a_Norm_Linf(i)
!    end do
!    write(*,*) "   L1 order       ", "   L2 order        ","Linf order "
!    do i = 2,counter-1
!        order_L1(i) = log(a_Norm_L1(i-1)/a_Norm_L1(i))/log(2.0_prec)
!        order_L2(i) = log(a_Norm_L2(i-1)/a_Norm_L2(i))/log(2.0_prec)
!        order_Linf(i) = log(a_Norm_Linf(i-1)/a_Norm_Linf(i))/log(2.0_prec)
!        write(*,*) order_L1(i),order_L2(i),order_Linf(i)
!    end do
!
!    write(*,*) "m      ","   L1 error       ","   L1 order       ","Linf error    ","      Linf order"
!    do i = 1,counter-1
!        write(*,"(1x,I3,E12.3,F10.3,E12.3,F10.3,E12.3,F10.3 )") nsdx/2**(counter-i),a_Norm_L1(i),order_L1(i),a_Norm_L2(i),order_L2(i),a_Norm_Linf(i),order_Linf(i)
!    end do
!!!!!!!!!!Print error result to '*.dat' file    
!    write(char_T,'(F5.3)') T !Tд��Ϊ�ַ���char_T
!    !write(char_nsd,'(I4)') nsd
!    file_name = 'order_'//'T='//char_T// '.txt'
!    open(unit=1001,status='REPLACE',file=file_name)
!    write(1001,*) "Error"
!    write(1001,*) " "
!    write(1001,*) "   L1 error       ","Linf error"
!    do i = 1,counter-1
!        write(1001,*) a_Norm_L1(i),a_Norm_Linf(i)
!    end do
!    write(1001,*) "   L1 order       ","Linf order"
!    do i = 2,counter-1
!        write(1001,*) order_L1(i),order_Linf(i)
!    end do
!    write(1001,*) "m      ","   L1 error       ","   L1 order       ","Linf error    ","      Linf order"
!    !write(1001,*) m,a_Norm_L1(i),a_Norm_Linf(i)
!    do i = 1,counter-1
!        write(1001,"(1x,I3,E20.7,F20.7,E20.7,F20.7 )") nsdx/2**(counter-i),a_Norm_L1(i),order_L1(i),a_Norm_Linf(i),order_Linf(i)
!    end do
!    close(1001,status="keep")
!end subroutine solve_error_order
    !subroutine IC_func(coor,ruvwp)
!    ! Ϊ��ʹT=P/R,ȡgamma*Ma*Ma = 1
!    use real_precision
!    use global_var,only: gamma,gamma1,u00,v00,xc0,yc0,xlong,ylong,pi
!    implicit none
!    real(prec),dimension(:) :: coor(2)       
!    real(prec) :: xx,yy,r,u,v,w,p,ruvwp(1:4)
!    real(prec) :: xcore,ycore,EKC,RC,RFA,RB2,RB,TAOB,TAOB2,temexp,TE,tt
!    
!    xx = coor(1)
!    yy = coor(2)
!    EKC=5.0_prec
!
!    if(xc0>=0.0_prec)then
!        if (xx>=(xc0-0.5_prec*xlong)) then
!            xcore = xc0
!        else
!            xcore = xc0 - xlong
!        end if
!    else
!        if (xx<=(xc0+0.5_prec*xlong)) then
!            xcore = xc0
!        else
!            xcore = xc0 + xlong
!        end if
!    end if
!
!
!    if(yc0>=0.0_prec)then
!        if (yy>=(yc0-0.5_prec*ylong)) then
!            ycore = yc0
!        else
!            ycore = yc0 - ylong
!        end if
!    else
!        if (yy<=(yc0+0.5_prec*ylong)) then
!            ycore = yc0
!        else
!            ycore = yc0 + ylong
!        end if
!    end if
! 
!
!    RB2 =  (xx-xcore)*(xx-xcore) + (yy-ycore)*(yy-ycore)
!    RB  = sqrt(RB2)
!    temexp = exp(0.5_prec*(1.0_prec-RB2))
!
!    TE = -(gamma-1.0_prec)*EKC*EKC*temexp*temexp/(8.0_prec*gamma*pi*pi)
!
!    tt = 1.0_prec + TE
!    r = tt**(1.0_prec/(gamma-1.0_prec))
!    p = r**gamma
!    u = u00 -EKC*temexp*(yy-ycore)/(2.0_prec*pi)
!    v = v00 +EKC*temexp*(xx-xcore)/(2.0_prec*pi)
!    w = 0.0_prec
! 
!    ruvwp(1)=r
!    ruvwp(2)=u
!    ruvwp(3)=v
!    ruvwp(4)=p
!    !ruvwp = 1.0_prec
!end subroutine IC_func
    
                        !!��ʣ��߱��
                    !do k = 1,4
                    !    if(node1 == cellset(cell1).nodes(k))then
                    !        do j = 1,4
                    !            if(j /= k)then
                    !                if(node2 == cellset(cell1).nodes(j))then
                    !                    if((k==1.and.j==4) .or.(k==4.and.j==1))then
                    !                        cellset(cell1).sides(4)=i
                    !                    else
                    !                        cellset(cell1).sides(4)=min(k,j)
                    !                    end if
                    !                end if
                    !            end if
                    !        end do
                    !    end if      
                    !end do

!--------------------------------------------------------------------------------------------------
!module RHS_module
!    use real_precision
!    use parameter_setting,only:nsp
!    implicit none
!    real(prec),dimension(:,:,:),allocatable :: FG_upw
!    real(prec),dimension(:,:,:),allocatable :: O1,O2,Q1,F1,G1
!    real(prec),dimension(:,:),allocatable  :: detJ
!    real(prec),dimension(4) :: fu_sub_kesi,gu_sub_eta
!end module RHS_module
!    
!subroutine RHS_CPR(index,Lsub)
!    !�ⲿ�������һ��cell�����н����Ҷ������Lsub��
!    use real_precision
!    use parameter_setting,only:nsp
!    use RHS_module
!    implicit none
!    integer index,i,j,k,n1,n2    
!    real(prec),dimension(:,:,:) :: Lsub(nsp,nsp,4)
!    
!    allocate(FG_upw(4,nsp,4),O1(nsp,nsp,4),O2(nsp,nsp,4),Q1(nsp,nsp,4),F1(nsp,nsp,4),G1(nsp,nsp,4),detJ(nsp,nsp))
!    
!    call get_local_fpvalue(index)   !���ͨ����ӭ��ͨ��
!    call get_RHS_dif(index,Lsub)    !�Ҷ���
!    
!    deallocate(FG_upw,O1,O2,Q1,F1,G1,detJ)
!    
!end subroutine RHS_CPR
!
!subroutine get_local_fpvalue(index)
!    !(���ͨ��һ��Ԫһ���£����ô�������ӭ��ͨ�����Դ������������ظ�����)
!    use parameter_setting,only:nsp
!     
!    use type_module
!    use RHS_module
!    use global_var
!    implicit none
!    integer index,i,j,k,l,m,n
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(1:2) :: norm
!    real(prec),dimension(1:4) :: u1,u2,fg1,fg2,ori1,ori2
!    integer :: sideIndex,nearCellIndex
!    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4)
!    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2
!    real(prec),dimension(:,:) :: quad_vertex(4,2) 
!   
!    do i =1,4 !side
!        sideIndex = cellset(index).sides(i)
!        !if(sideset(sideindex).fpvalue_upw(1,1) /= 1000.0_prec) cycle    !ÿ��ʱ�䲽���ӭ��ͨ����ʼ��1000������Ϊ1000����ǰ�沽�����󣬲����ظ�����
!        nearCellIndex = cellset(index).nearcells(i)
!        !�������±߽���ͨ��
!        !1��ʾ��Ԫ�ڣ�2��ʾ�ٵ�Ԫ
!        O1 = cellset(index).spvalue_ori
!        O2 = cellset(nearCellIndex).spvalue_ori
!        Q1 = cellset(index).spvalue_con_loc
!        F1 = cellset(index).spvalue_fluF_loc 
!        G1 = cellset(index).spvalue_fluG_loc
!        !ȡ���ı��ε�Ԫ��������ֵ
!        do j = 1,4
!            do k = 1,2
!                quad_vertex(j,k) = xy_coor(cellset(index).nodes(j),k)
!            end do
!        end do
!        !write(*,*) index,sideIndex,nearCellIndex
!        do l = 1,nsp
!            do m =1,4
!                if(i == 1)then                  
!                    ori1(m) = LaI_nPs(kesi_l,SPs,O1(:,l,m),nsp)
!                    !ori2(m) = LaI_nPs(kesi_r,SPs,O2(:,l,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,ori2)                
!                else if(i == 2 )then
!                    ori1(m) = LaI_nPs(kesi_r,SPs,O1(l,:,m),nsp)
!                    !ori2(m) = LaI_nPs(kesi_l,SPs,O2(l,:,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,ori2)
!                    !write(*,*)'f',fg1(m),fg2(m)
!                else if(i == 3)then                  
!                    ori1(m) = LaI_nPs(kesi_r,SPs,O1(:,l,m),nsp)
!                    !ori2(m) = LaI_nPs(kesi_l,SPs,O2(:,l,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,ori2)
!                else if(i == 4 )then
!                    ori1(m) = LaI_nPs(kesi_l,SPs,O1(l,:,m),nsp)
!                    !ori2(m) = LaI_nPs(kesi_r,SPs,O2(l,:,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,ori2)
!                end if    
!            end do
!            
!            !���Ʒ���
!            if(i==1)then
!                n1 = 0.0_prec
!                n2 = -1.0_prec         
!            else if(i==2)then
!                n1 = 1.0_prec
!                n2 = 0.0_prec
!            else if(i==3)then
!                n1 = 0.0_prec
!                n2 = 1.0_prec
!            else if(i==4)then
!                n1 = -1.0_prec
!                n2 = 0.0_prec
!            end if
!
!            if(i==1)then
!                coor_C(1) = SPs(l)
!                coor_C(2) = -1.0_prec
!            else if(i==2)then
!                coor_C(1) = 1.0_prec
!                coor_C(2) = SPs(l)
!            else if(i == 3)then
!                coor_C(1) = SPs(l)
!                coor_C(2) = 1.0_prec
!            else if(i == 4)then
!                coor_C(1) = -1.0_prec
!                coor_C(2) = SPs(l)
!            end if
!                
!            call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)
!
!            kesix = FPdirect(1)
!            kesiy = FPdirect(2)
!            etax  = FPdirect(3)
!            etay  = FPdirect(4)
!            norm(1) = kesix*n1+etax*n2
!            norm(2) = kesiy*n1+etay*n2
!            !norm = norm*FPdetJ
!            !write(*,*)index,i,l
!            call laxf_flux(ori1,ori2,norm,FG_upw(i,l,:))
!            FG_upw(i,l,:)=FG_upw(i,l,:)*FPdetJ
!            !write(*,"(I3,3F10.5)")index,ori1(1),ori2(1),FG_upw(i,l,1)
!        end do            
!        if(i == 1 .or. i == 4)then
!            FG_upw(i,:,:) = -1.0_prec*FG_upw(i,:,:) 
!        end if
!        !sideset(sideindex).fpvalue_upw(:,:) = FG_upw(i,:,:)!i-th side, l-th fp ,m-th var
!        !write(*,*)sideIndex,sideset(sideindex).fpvalue_upw(:,:)
!        !write(*,*)FG_upw
!    end do 
!     
!end subroutine get_local_fpvalue
!
!subroutine get_ori2(index,nearCellIndex,i,l,m,ori2)
!    use global_var
!    use parameter_setting
!    use type_module
!    use RHS_module
!    implicit none
!    
!    integer :: index,nearCellIndex,i,l,m
!    real(prec) :: ori2(4)
!    real(prec),external :: LaI_nPs
!    if(i==1)then
!        if(cellset(nearCellIndex).nearcells(2)==index)then
!            ori2(m) = LaI_nPs(kesi_r,SPs,O2(nsp+1-l,:,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(3)==index)then
!            ori2(m) = LaI_nPs(kesi_r,SPs,O2(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(4)==index)then
!            ori2(m) = LaI_nPs(kesi_l,SPs,O2(l,:,m),nsp)
!        end if
!    elseif(i==2)then
!        if(cellset(nearCellIndex).nearcells(1)==index)then
!            ori2(m) = LaI_nPs(kesi_l,SPs,O2(:,nsp+1-l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(3)==index)then
!            ori2(m) = LaI_nPs(kesi_r,SPs,O2(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(4)==index)then
!            ori2(m) = LaI_nPs(kesi_l,SPs,O2(l,:,m),nsp)
!        end if        
!    elseif(i==3)then
!        if(cellset(nearCellIndex).nearcells(1)==index)then
!            ori2(m) = LaI_nPs(kesi_l,SPs,O2(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(2)==index)then
!            ori2(m) = LaI_nPs(kesi_r,SPs,O2(l,:,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(4)==index)then
!            ori2(m) = LaI_nPs(kesi_l,SPs,O2(nsp+1-l,:,m),nsp)
!        end if          
!    elseif(i==4)then
!        if(cellset(nearCellIndex).nearcells(1)==index)then
!            ori2(m) = LaI_nPs(kesi_l,SPs,O2(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(2)==index)then
!            ori2(m) = LaI_nPs(kesi_r,SPs,O2(l,:,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(3)==index)then
!            ori2(m) = LaI_nPs(kesi_r,SPs,O2(:,nsp+1-l,m),nsp)
!        end if         
!    end if
!    
!end subroutine get_ori2
!
!subroutine get_RHS_dif(index,L_sub_global)
!    use global_var
!    use parameter_setting
!    use type_module
!    use RHS_module
!    implicit none
!    integer k,l,m,index
!    real(prec),external :: LaI_nPs,LaI_nPs_deri,gl_sub_kesi,gr_sub_kesi
!    real(prec),dimension(4) :: f_l,f_r,g_l,g_r
!    real(prec),dimension(:,:,:) :: L_sub_global(nsp,nsp,4),F_sub_kesi(nsp,nsp,4),G_sub_eta(nsp,nsp,4)
!    
!    detJ = cellset(index).det_J
!    do k = 1,nsp
!        do l =1,nsp
!        
!            do m =1,4              
!                f_l(m) = LaI_nPs(kesi_l,SPs,F1(k,:,m),nsp)!ͨ����ֵ
!                f_r(m) = LaI_nPs(kesi_r,SPs,F1(k,:,m),nsp)!ͬһ��Ԫ
!                g_l(m) = LaI_nPs(kesi_l,SPs,G1(:,l,m),nsp)
!                g_r(m) = LaI_nPs(kesi_r,SPs,G1(:,l,m),nsp)
!            end do  
!            
!            call get_fluxDer(k,l)
!                    
!            F_sub_kesi(k,l,:) = fu_sub_kesi(:)+(FG_upw(4,k,:)-f_l(:))*gl_sub_kesi(SPs(l)) + (FG_upw(2,k,:)-f_r(:))*gr_sub_kesi(SPs(l))
!            G_sub_eta(k,l,:)  = gu_sub_eta(:) +(FG_upw(1,l,:)-g_l(:))*gl_sub_kesi(SPs(k)) + (FG_upw(3,l,:)-g_r(:))*gr_sub_kesi(SPs(k))                
!        end do
!    end do
!    do m =1,4
!        L_sub_global(:,:,m) = -(F_sub_kesi(:,:,m)+G_sub_eta(:,:,m))/detJ(:,:)
!    end do
!end subroutine
!
!subroutine get_fluxDer(k,l)
!    !���������
!    use RHS_module
!    use real_precision
!    use parameter_setting
!    use global_var
!    implicit none
!    real(prec),external :: LaI_nPs,LaI_nPs_deri
!    integer m,k,l
!    real(prec) :: r,u,v,e
!    real(prec),dimension(4,4) :: Fdq,Gdq
!    real(prec),dimension(4)   :: Qder1,Qder2 
!    if(fluxD_type == LP)then
!        do m = 1,4
!            fu_sub_kesi(m) = LaI_nPs_deri(SPs(l),SPs,F1(k,:,m),nsp)
!            gu_sub_eta(m)  = LaI_nPs_deri(SPs(k),SPs,G1(:,l,m),nsp)
!        end do
!    elseif(fluxD_type == CR)then
!        r = Q1(k,l,1)
!        u = Q1(k,l,2)/r
!        v = Q1(k,l,3)/r
!        e = Q1(k,l,4)
!        !A(U)
!        Fdq(1,1) = 0.0_prec
!        Fdq(1,2) = 1.0_prec
!        Fdq(1,3) = 0.0_prec
!        Fdq(1,4) = 0.0_prec
!        Fdq(2,1) = 0.5_prec*((gamma-3.0_prec)*u**2+(gamma-1.0_prec)*v**2)
!        Fdq(2,2) = (3.0_prec-gamma)*u
!        Fdq(2,3) = (1.0_prec-gamma)*v
!        Fdq(2,4) = gamma-1.0_prec
!        Fdq(3,1) = -u*v
!        Fdq(3,2) = v
!        Fdq(3,3) = u
!        Fdq(3,4) = 0.0_prec
!        Fdq(4,1) = -gamma*e*u/r+(gamma-1.0_prec)*(u**3+u*v**2)
!        Fdq(4,2) = gamma*e/r-0.5_prec*(gamma-1.0_prec)*(3.0_prec*u**2+v**2)
!        Fdq(4,3) = (1.0_prec-gamma)*u*v
!        Fdq(4,4) = gamma*u
!        !B(U)
!        Gdq(1,1) = 0.0_prec
!        Gdq(1,2) = 0.0_prec
!        Gdq(1,3) = 1.0_prec
!        Gdq(1,4) = 0.0_prec
!        Gdq(2,1) = -u*v
!        Gdq(2,2) = v
!        Gdq(2,3) = u
!        Gdq(2,4) = 0.0_prec
!        Gdq(3,1) = 0.5_prec*((gamma-3.0_prec)*v**2+(gamma-1.0_prec)*u**2)
!        Gdq(3,2) = (1.0_prec-gamma)*u
!        Gdq(3,3) = (3.0_prec-gamma)*v
!        Gdq(3,4) = gamma-1.0_prec
!        Gdq(4,1) = -gamma*e*v/r+(gamma-1.0_prec)*(v**3+v*u**2)
!        Gdq(4,2) = (1-gamma)*u*v
!        Gdq(4,3) = gamma*e/r-0.5_prec*(gamma-1.0_prec)*(3.0_prec*v**2+u**2)
!        Gdq(4,4) = gamma*v
!        do m = 1,4
!            Qder1(m) = LaI_nPs_deri(SPs(l),SPs,Q1(k,:,m),nsp)
!            Qder2(m) = LaI_nPs_deri(SPs(k),SPs,Q1(:,l,m),nsp)
!        end do
!        do m = 1,4
!            fu_sub_kesi(m) = dot_product(Fdq(m,:),Qder1(:))
!            gu_sub_eta(m)  = dot_product(Gdq(m,:),Qder2(:))
!        end do
!    end if
!
!end subroutine get_fluxDer
!
!function gl_sub_kesi(kesi)
!    !!gl = (-1)**K/2*(P_k-P_k-1),Pk is Legendre polynomial,��˽�������
!    !!K=3,gl = -1/4*(5*x**3-3*x**2-3*x+1),gl_sub_kesi = -1/4*(15*x**2-6*x-3)
!    !!K=4,gl = 1/16*(35*x**4-20*x**3-30*x**2+12*x+3),gl_sub_kesi = 1/16*(140*x**3-60*x**2-60*x+12)
!    !!K=5,gl = -1/16*(63*x**5-35*x**4-70*x**3+30*x**2+15*x-3),gl_sub_kesi = -1/16*(315*x**4-140*x**3-210*x**2+60*x+15)
!    use real_precision
!    use parameter_setting ,only:nsp
!    implicit none
!    real(prec) :: gl_sub_kesi,kesi
!    if (nsp .EQ. 3)then
!        gl_sub_kesi = -0.25_prec*(15.0_prec*kesi**2.0_prec-6.0_prec*kesi-3.0_prec)
!    elseif(nsp .EQ. 4)then
!        gl_sub_kesi = (140.0_prec*kesi**3-60.0_prec*kesi**2-60.0_prec*kesi+12.0_prec)/16.0_prec
!    elseif(nsp .EQ. 5)then
!        gl_sub_kesi = (315.0_prec*kesi**4-140.0_prec*kesi**3-210.0_prec*kesi**2+60.0_prec*kesi+15.0_prec)*(-1.0_prec)/16.0_prec
!    endif
!    return
!end function gl_sub_kesi
!
!function gr_sub_kesi(kesi)
!    !!gr = 1/2*(P_k+P_k-1),�Ҷ˽�������
!    !!K=3,gr = 1/4*(5*x**3+3*x**2-3*x-1),gr_sub_kesi = 1/4*(15*x**2+6*x-3)
!    !!K=4,gl = 1/16*(35*x**4+20*x**3-30*x**2-12*x+3),gl_sub_kesi = 1/16*(140*x**3+60*x**2-60*x-12)
!    !!K=5,gl = 1/16*(63*x**5+35*x**4-70*x**3-30*x**2+15*x+3),gl_sub_kesi = 1/16*(315*x**4+140*x**3-210*x**2-60*x+15)
!    use real_precision
!    use parameter_setting ,only:nsp
!    implicit none
!    real(prec) :: gr_sub_kesi,kesi
!
!    if (nsp .EQ. 3)then
!        gr_sub_kesi = 0.25_prec*(15.0_prec*kesi**2.0_prec+6.0_prec*kesi-3.0_prec)
!    elseif(nsp .EQ.4)then
!        gr_sub_kesi = (140.0_prec*kesi**3+60.0_prec*kesi**2-60.0_prec*kesi-12.0_prec)/16.0_prec
!    elseif(nsp .EQ. 5)then
!        gr_sub_kesi = (315.0_prec*kesi**4+140.0_prec*kesi**3-210.0_prec*kesi**2-60.0_prec*kesi+15.0_prec)/16.0_prec
!    endif
!    return
!end function gr_sub_kesi
    
!subroutine get_local_fpvalue(index)
!    !(���ͨ��һ��Ԫһ���£����ô�������ӭ��ͨ�����Դ������������ظ�����)
!    use parameter_setting,only:nsp
!    use type_module
!    use RHS_module
!    use global_var
!    implicit none
!    integer index,i,j,k,l,m,n
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(1:2) :: norm
!    real(prec),dimension(1:4) :: u1,u2,fg1,fg2,ori1,ori2
!    integer :: sideIndex,nearCellIndex
!    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4)
!    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2
!    real(prec),dimension(:,:) :: quad_vertex(4,2) 
!   
!    do i =1,4 !side
!        sideIndex = cellset(index).sides(i)
!        !if(sideset(sideindex).fpvalue_upw(1,1) /= 1000.0_prec) cycle    !ÿ��ʱ�䲽���ӭ��ͨ����ʼ��1000������Ϊ1000����ǰ�沽�����󣬲����ظ�����
!        nearCellIndex = cellset(index).nearcells(i)
!        !�������±߽���ͨ��
!        !1��ʾ��Ԫ�ڣ�2��ʾ�ٵ�Ԫ
!        O1 = cellset(index).spvalue_ori
!        O2 = cellset(nearCellIndex).spvalue_ori
!        Q1 = cellset(index).spvalue_con_loc
!        F1 = cellset(index).spvalue_fluF_loc 
!        G1 = cellset(index).spvalue_fluG_loc
!        !ȡ���ı��ε�Ԫ��������ֵ
!        do j = 1,4
!            do k = 1,2
!                quad_vertex(j,k) = xy_coor(cellset(index).nodes(j),k)
!            end do
!        end do
!        !write(*,*) index,sideIndex,nearCellIndex
!        do l = 1,nsp
!            do m =1,4
!                if(i == 1)then                  
!                    ori1(m) = LaI_nPs(kesi_l,SPs,O1(:,l,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,O2,ori2)                
!                else if(i == 2 )then
!                    ori1(m) = LaI_nPs(kesi_r,SPs,O1(l,:,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,O2,ori2)
!                    !write(*,*)'f',fg1(m),fg2(m)
!                else if(i == 3)then                  
!                    ori1(m) = LaI_nPs(kesi_r,SPs,O1(:,l,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,O2,ori2)
!                else if(i == 4 )then
!                    ori1(m) = LaI_nPs(kesi_l,SPs,O1(l,:,m),nsp)
!                    call  get_ori2(index,nearCellIndex,i,l,m,O2,ori2)
!                end if    
!            end do
!            
!            !���Ʒ���
!            if(i==1)then
!                n1 = 0.0_prec
!                n2 = -1.0_prec         
!            else if(i==2)then
!                n1 = 1.0_prec
!                n2 = 0.0_prec
!            else if(i==3)then
!                n1 = 0.0_prec
!                n2 = 1.0_prec
!            else if(i==4)then
!                n1 = -1.0_prec
!                n2 = 0.0_prec
!            end if
!
!            if(i==1)then
!                coor_C(1) = SPs(l)
!                coor_C(2) = -1.0_prec
!            else if(i==2)then
!                coor_C(1) = 1.0_prec
!                coor_C(2) = SPs(l)
!            else if(i == 3)then
!                coor_C(1) = SPs(l)
!                coor_C(2) = 1.0_prec
!            else if(i == 4)then
!                coor_C(1) = -1.0_prec
!                coor_C(2) = SPs(l)
!            end if
!                
!            call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)
!
!            kesix = FPdirect(1)
!            kesiy = FPdirect(2)
!            etax  = FPdirect(3)
!            etay  = FPdirect(4)
!            norm(1) = kesix*n1+etax*n2
!            norm(2) = kesiy*n1+etay*n2
!            !norm = norm*FPdetJ
!            !write(*,*)index,i,l
!            call laxf_flux(ori1,ori2,norm,FG_upw(i,l,:))
!            FG_upw(i,l,:)=FG_upw(i,l,:)*FPdetJ
!            !write(*,"(I3,3F10.5)")index,ori1(1),ori2(1),FG_upw(i,l,1)
!        end do            
!        if(i == 1 .or. i == 4)then
!            FG_upw(i,:,:) = -1.0_prec*FG_upw(i,:,:) 
!        end if
!        !sideset(sideindex).fpvalue_upw(:,:) = FG_upw(i,:,:)!i-th side, l-th fp ,m-th var
!        !write(*,*)sideIndex,sideset(sideindex).fpvalue_upw(:,:)
!        !write(*,*)FG_upw
!    end do 
!     
!end subroutine get_local_fpvalue
    
!    subroutine Face_Flux_upw
!!���㵥Ԫ����ϵ�ӭ��ͨ��
!    use parameter_setting
!    use global_var
!    use type_module
!    use real_precision
!    implicit none
!    integer :: i,j,m,k,l,LC_sideth
!    integer :: indexCellL,indexCellR
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(1:2) :: norm
!    real(prec),dimension(1:4) :: u1,u2,fg1,fg2,ori1,ori2
!    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4)
!    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2
!    real(prec),dimension(:,:) :: quad_vertex(4,2)
!    real(prec),dimension(:,:,:) :: varOri(nsp,nsp,4),varOri2(nsp,nsp,4)
!    do i = 1,nsides
!        indexCellL = sideset(i).nearcells(1)
!        indexCellR = sideset(i).nearcells(2)
!        if(indexCellL==0)then
!            !Ϊ0˵���Ǳ߽��ϵıߣ�����ݱ߽�����ȡ��Ԫ
!            do j = 1,4
!                if(cellset(indexCellR).sides(j)==i)then
!                    indexCellL = cellset(indexCellR).nearcells(j)
!                    exit
!                end if
!            end do           
!        else if(indexCellR==0)then
!            !Ϊ0˵���Ǳ߽��ϵıߣ�����ݱ߽�����ȡ�ҵ�Ԫ
!            do j = 1,4
!                if(cellset(indexCellL).sides(j)==i)then
!                    indexCellR = cellset(indexCellL).nearcells(j)
!                    exit
!                end if
!            end do 
!        end if
!        !�ж�����������Ԫ�ĵڼ������
!        do j = 1,4
!            if(cellset(indexCellL).sides(j) == i)then
!                LC_sideth = j
!                exit
!            end if
!        end do
!        !������������ֵ��ԭʼ����      
!        varOri(:,:,:) = cellset(indexCellL).spvalue_ori(:,:,:)
!        varOri2(:,:,:)= cellset(indexCellR).spvalue_ori(:,:,:)
!        
!        do l = 1,nsp
!            do m =1,4
!                if(LC_sideth == 1)then                                      
!                    ori1(m) = LaI_nPs(kesi_l,SPs,varOri(:,l,m),nsp)
!                    call get_ori2(indexCellL,indexCellR,LC_sideth,l,m,varOri2,ori2)                
!                else if(LC_sideth == 2 )then
!                    ori1(m) = LaI_nPs(kesi_r,SPs,varOri(l,:,m),nsp)
!                    call get_ori2(indexCellL,indexCellR,LC_sideth,l,m,varOri2,ori2)
!                    !write(*,*)'f',fg1(m),fg2(m)
!                else if(LC_sideth == 3)then                  
!                    ori1(m) = LaI_nPs(kesi_r,SPs,varOri(:,l,m),nsp)
!                    call get_ori2(indexCellL,indexCellR,LC_sideth,l,m,varOri2,ori2)
!                else if(LC_sideth == 4 )then
!                    ori1(m) = LaI_nPs(kesi_l,SPs,varOri(l,:,m),nsp)
!                    call get_ori2(indexCellL,indexCellR,LC_sideth,l,m,varOri2,ori2)
!                end if    
!            end do
!            
!            !���Ʒ���
!            if(LC_sideth == 1)then
!                n1 = 0.0_prec
!                n2 = -1.0_prec         
!            else if(LC_sideth == 2)then
!                n1 = 1.0_prec
!                n2 = 0.0_prec
!            else if(LC_sideth == 3)then
!                n1 = 0.0_prec
!                n2 = 1.0_prec
!            else if(LC_sideth == 4)then
!                n1 = -1.0_prec
!                n2 = 0.0_prec
!            end if
!            
!            if(LC_sideth==1)then
!                coor_C(1) = SPs(l)
!                coor_C(2) = -1.0_prec
!            else if(LC_sideth == 2)then
!                coor_C(1) = 1.0_prec
!                coor_C(2) = SPs(l)
!            else if(LC_sideth == 3)then
!                coor_C(1) = SPs(l)
!                coor_C(2) = 1.0_prec
!            else if(LC_sideth == 4)then
!                coor_C(1) = -1.0_prec
!                coor_C(2) = SPs(l)
!            end if
!            
!            !ȡ�����ı��ε�Ԫ��������ֵ
!            do j = 1,4
!                do k = 1,2
!                    quad_vertex(j,k) = xy_coor(cellset(indexCellL).nodes(j),k)
!                end do
!            end do    
!            call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)
!            
!            kesix = FPdirect(1)
!            kesiy = FPdirect(2)
!            etax  = FPdirect(3)
!            etay  = FPdirect(4)
!            norm(1) = kesix*n1+etax*n2
!            norm(2) = kesiy*n1+etay*n2
!            !norm = norm*FPdetJ
!            !write(*,*)index,i,l
!            call laxf_flux(ori1,ori2,norm,sideset(i).fpvalue_upw(l,:))
!            sideset(i).fpvalue_upw(l,:) = sideset(i).fpvalue_upw(l,:)*FPdetJ
!            !write(*,"(I3,3F10.5)")index,ori1(1),ori2(1),FG_upw(i,l,1)   
!        end do
!                          
!        if(LC_sideth == 1 .or. LC_sideth == 4)then
!            sideset(i).fpvalue_upw(:,:) = -1.0_prec*sideset(i).fpvalue_upw(:,:)
!        end if            
!    end do
!
!end subroutine Face_Flux_upw
    
function lim2(Var,VarC,varMax,varMin)
    use real_precision
    implicit none
    
    real(prec) :: lim2,Var,VarC,varMax,varMin,temp
    
    if(VarC<Var)then
        temp = (Var-VarC)/(varMax-VarC)
        lim2 = min(1.0_prec,temp)
    elseif(VarC>Var)then
        temp = (Var-VarC)/(varMin-VarC)
        lim2 = min(1.0_prec,temp)
    elseif(Var==VarC)then
        lim2 = 1.0_prec
    endif

    return
    end function lim2

!subroutine C2NNW2_LR(xu,xf,u,u_cell_L,u_cell_R)
!    use real_precision
!    implicit none
!    real(prec),dimension(:) :: xu(1:3),xf(1:2),u(1:3)
!    real(prec) :: u1,u2,u3,u_cell_L,u_cell_R
!    real(prec) :: d1,d2,d3,d4,dd1,dd2,dd3,dd4
!    real(prec) :: cc,w1,w2,w3,w4,w5,w6
!    real(prec) :: uA1,uB1,uA2,uB2,du1,du2,du
!    real(prec) :: limA,limB,Fai,m1,m2 
!    real(prec),external :: lim,lim2 
!    !ͨ��������֮��ľ���
!    d1 = xf(1)-xu(1)
!    d2 = xu(2)-xf(1)
!    d3 = xf(2)-xu(2)
!    d4 = xu(3)-xf(2)
!    
!    !����Ȩ
!    cc = 1.0_prec!����1������д���������滻
!    dd1 = cc/d1
!    dd2 = cc/d2
!    dd3 = cc/d3
!    dd4 = cc/d4
!    
!    w1 = dd1/(dd1+dd2)
!    w2 = dd2/(dd1+dd2)
!    w3 = dd3/(dd3+dd4)
!    w4 = dd4/(dd3+dd4)
!    
!    !ͨ���㴦ֵ
!    uA1 = w1*u(1)+w2*u(2)
!    uB1 = w3*u(2)+w4*u(3)
!    
!    !��Ԫ����
!    w5 = dd2/(dd2+dd3)
!    w6 = dd3/(dd2+dd3)
!
!    du1 = (u(2)-uA1)/d2
!    du2 = (uB1-u(2))/d3
!    du = w5*du1 + w6*du2
!    !write(*,*)w1,w2,w3,w4,w5,w6
!    !���¼���uA,uB
!    uA2 = u(2)-du*d2
!    uB2 = u(2)+du*d3
!    m1 = maxVal(u)
!    m2 = minVal(u)
!
!    limA = lim2(u(2),uA2,m1,m2)
!    limB = lim2(u(2),uB2,m1,m2)
!    !write(*,*)'lim',du,limA,limB
!    Fai = min(limA,limB)
!    
!    !��Ԫ��ֵ�ǽ�����ֵ����Ԫ��ֵ�ǽ�����ֵ
!    u_cell_L = u(2)-Fai*du1*d2
!    u_cell_R = u(2)+Fai*du2*d3
!    
!end subroutine C2NNW2_LR
  !if(grid_set == self)then
    !    select case(bound_set)
    !    case(bound_cyc)     !������������ֻ��nsdyΪż��ʱ������ѭ���߽��ͨ�����غ�,����Ϊż������Ҫ������������߽���λ�ã������²�ֵ�õ�ֵ           
    !            do i = 1,ncells
    !                if(cellset(i).nearcells(1)==0)then
    !                    cellset(i).nearcells(1) = cellset(i+(nsdpx-1)*(nsdpy-2)).index 
    !                    cellset(i+(nsdpx-1)*(nsdpy-2)).nearcells(3) = cellset(i).index
    !                end if
    !                if(cellset(i).nearcells(2)==0)then
    !                    cellset(i).nearcells(2) = cellset(i-nsdpx+2).index
    !                    cellset(i-nsdpx+2).nearcells(4) = cellset(i).index
    !                end if
    !            end do          
    !    case(bound_non)
    !    
    !        do i = 1,ncells
    !            if(cellset(i).nearcells(1)==0)then
    !                cellset(i).nearcells(1) = cellset(i).index
    !            end if
    !            if(cellset(i).nearcells(2)==0)then
    !                cellset(i).nearcells(2) = cellset(i).index
    !            end if
    !            if(cellset(i).nearcells(3)==0)then
    !                cellset(i).nearcells(3) = cellset(i).index
    !            end if
    !            if(cellset(i).nearcells(4)==0)then
    !                cellset(i).nearcells(4) = cellset(i).index
    !            end if
    !        end do
    !    case(99)
    !        do i = 1,ncells
    !            if(cellset(i).nearcells(1)==0)then
    !                cellset(i).nearcells(1) = cellset(i+(nsdpx-1)*(nsdpy-2)).index 
    !                cellset(i+(nsdpx-1)*(nsdpy-2)).nearcells(3) = cellset(i).index
    !            end if
    !            if(cellset(i).nearcells(2)==0)then
    !                cellset(i).nearcells(2) = cellset(i).index
    !            end if
    !            if(cellset(i).nearcells(4)==0)then
    !                cellset(i).nearcells(4) = cellset(i).index
    !            end if
    !        end do
    !    case(100)
    !        do i = 1,ncells
    !            if(cellset(i).nearcells(2)==0)then
    !                cellset(i).nearcells(2) = cellset(i-nsdpx+2).index
    !                cellset(i-nsdpx+2).nearcells(4) = cellset(i).index
    !            end if
    !            if(cellset(i).nearcells(1)==0)then
    !                cellset(i).nearcells(1) = cellset(i).index
    !            end if                
    !            if(cellset(i).nearcells(3)==0)then
    !                cellset(i).nearcells(3) = cellset(i).index
    !            end if
    !        end do
    !    case default
    !        write(*,*) 'Error! Occur in subroutine set_boundary'
    !    end select
    !else if(grid_set == fluent_cas)then!cas����




        !do i = 1,ncells
        !    do j = 1,nsp
        !        do k = 1,nsp
        !            write(60,'(6F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori_exa(j,k,1),cellset(i).spvalue_ori_exa(j,k,2),cellset(i).spvalue_ori_exa(j,k,3),cellset(i).spvalue_ori_exa(j,k,4)
        !        end do
        !    end do
        !end do
        !do i = 1,ncells
        !    do j = 1,nsp
        !        do k = 1,nsp
        !            write(60,*) (i-1)*nsp*nsp+(j-1)*nsp+k,(i-1)*nsp*nsp+(j-1)*nsp+k,(i-1)*nsp*nsp+(j-1)*nsp+k,(i-1)*nsp*nsp+(j-1)*nsp+k
        !        end do
        !    end do
        !end do
!�����������ȷ���߽�����equEntropy_case = 0, SodShockTube_case = 1, LaxShockTube_case = 2,ShuOsher_case = 3,Riemann2D_case = 4
    !if(case_comp == equEntropy_case)then
    !    bound_set = bound_cyc
    !elseif(case_comp == SodShockTube_case)then
    !    select case(dire_shock)
    !    case(direct_x)
    !        bound_set = 99
    !    case(direct_y)
    !        bound_set = 100
    !    end select
    !elseif(case_comp == LaxShockTube_case)then
    !    select case(dire_shock)
    !    case(direct_x)
    !        bound_set = 99
    !    case(direct_y)
    !        bound_set = 100
    !    end select
    !elseif(case_comp == ShuOsher_case)then
    !    select case(dire_shock)
    !    case(direct_x)
    !        bound_set = 99
    !    case(direct_y)
    !        bound_set = 100
    !    end select
    !elseif(case_comp == Riemann2D_case)then
    !    bound_set = bound_non
    !else
    !    write(*,*)'other'!��˫��շ���ĳ�ʼ�߽�����
    !end if
    
 !do m = 1,4
 !       !����Ԫ����ƽ��ֵ��ԭʼ����
 !       aveCell = Gauss_double_integral(varOriCell(:,:,m))/4.0_prec
 !       !�ڵ�Ԫ����ƽ��ֵ
 !       do i = 1,4
 !           aveNearCell(i) = Gauss_double_integral(varOriCellNearCell(i,:,:,m))/4.0_prec
 !       end do
 !   
 !       !TVB_M = 0.01_prec
 !   
 !       do i = 1,nsp
 !           !dh2 = 1+SPs(i)**2   !4M*dh^2
 !           dh2 = 2.0_prec*cellset(index).MJacobi(i,1,1)
 !           !write(*,*)dh2
 !           varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(i,:,m),nsp)
 !           varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(i,:,m),nsp)        
 !           varmodL = TVB_minmod(aveCell-varFluxL,aveNearCell(2)-aveCell,aveCell-aveNearCell(4),dh2)
 !           varmodR = TVB_minmod(varFluxR-aveCell,aveNearCell(2)-aveCell,aveCell-aveNearCell(4),dh2)
 !           if(varmodL /= aveCell-varFluxL .or. varmodR /= varFluxR-aveCell)then
 !               cellset(index).Beta = 1
 !               exit
 !           end if
 !           dh2 = 2.0_prec*cellset(index).MJacobi(1,i,4) 
 !           varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(:,i,m),nsp)
 !           varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(:,i,m),nsp)  
 !           varmodL = TVB_minmod(aveCell-varFluxL,aveNearCell(3)-aveCell,aveCell-aveNearCell(1),dh2)
 !           varmodR = TVB_minmod(varFluxR-aveCell,aveNearCell(3)-aveCell,aveCell-aveNearCell(1),dh2)
 !           if(varmodL /= aveCell-varFluxL .or. varmodR /= varFluxR-aveCell)then
 !               cellset(index).Beta = 1
 !               write(*,*)'Trouble Cell',index
 !               exit
 !           end if     
 !       end do
 !   end do
 !do i = 1,ncells
        !���ⵥԪ����� -- ncells+1 : ncells+n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
    !    if(cellset(i).nearcells(1)==0)then
    !        cellset(i).nearcells(1) = ncells+index1
    !        index1 = index1+1
    !    end if
    !    if(cellset(i).nearcells(2)==0)then
    !        cellset(i).nearcells(2) = ncells+n_BCells_D+index2
    !        index2 = index2+1
    !    end if
    !    if(cellset(i).nearcells(3)==0)then
    !        cellset(i).nearcells(3) = ncells+n_BCells_D+n_BCells_R+index3
    !        index3 = index3+1
    !    end if
    !    if(cellset(i).nearcells(4)==0)then
    !        cellset(i).nearcells(4) = ncells+n_BCells_D+n_BCells_R+n_BCells_U+index4
    !        index4 = index4+1
    !    end if
    !end do
    !write(*,*)index1,index2,index3,index4
    !equEntropy_case = 0, SodShockTube_case = 1, LaxShockTube_case = 2,ShuOsher_case = 3,Riemann2D_case = 4
    !select case(case_comp)
    !    case(equEntropy_case)  !�����в��� ��  x y ѭ���߽�����             
    !        do i = 1,ncells
    !            if(cellset(i).nearcells(1)==0)then
    !                do j = 1,ncells                 
    !                    if(cellset(j).nearcells(3)==0)then
    !                        if(abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(4),1))<1.0e-10&
    !                            .and.abs(xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(3),1))<1.0e-10 )then
    !                            cellset(i).nearcells(1) = j
    !                            cellset(j).nearcells(3) = i
    !                            !write(*,*)'00',i,j
    !                        end if
    !                    endif
    !                end do
    !            end if
    !            if(cellset(i).nearcells(2)==0)then
    !                do j = 1,ncells
    !                    if(cellset(j).nearcells(4)==0)then
    !                        
    !                        if(abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))<1.0e-10&
    !                            .and.abs(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))<1.0e-10 )then
    !                            cellset(i).nearcells(2) = j
    !                            cellset(j).nearcells(4) = i
    !                            !sideset(cellset(i).sides(2)).nearcells(2) = j
    !                            !sideset(cellset(j).sides(4)).nearcells(2) = i
    !                            !write(*,*)'01',i,j
    !                        end if
    !                    endif
    !                end do
    !            end if
    !        end do
    !        !do i = 1,ncells
    !        !    write(*,*) cellset(i).nearcells(:)
    !        !end do
    !    case(bound_non)
    !        do i = 1,ncells
    !            do k =1,4
    !                if(cellset(i).nearcells(k)==0)then                       
    !                    cellset(i).nearcells(k) = i    
    !                    
    !                end if    
    !            end do                    
    !        end do         
    !    case(99)
    !        do i = 1,ncells
    !            if(cellset(i).nearcells(1)==0)then                   
    !                do j = 1,ncells
    !                    
    !                    if(cellset(j).nearcells(3)==0 .and.abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(4),1))<1.0e-10&
    !                        .and.abs(xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(3),1))<1.0e-10 )then
    !                        !write(*,*) i,j
    !                        cellset(i).nearcells(1) = j
    !                        cellset(j).nearcells(3) = i
    !                    endif
    !                end do
    !            end if
    !            if(cellset(i).nearcells(2)==0)then
    !                !write(*,*) i
    !                cellset(i).nearcells(2) = i
    !                !write(*,*) i,cellset(i).nearcells(2)
    !                do j = 1,ncells!�����ڲ��߽�㣬����������С���������򽻽��бߣ��ҽڵ��غϣ�������������������������x���򴫲��ļ��������⣩
    !                    dis1 = (xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(1),1))**2+(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))**2
    !                    dis2 = (xy_coor(cellset(i).nodes(3),1)-xy_coor(cellset(j).nodes(4),1))**2+(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))**2
    !                    
    !                    if(dis1+dis2 < 1.0e-9)then                                                         
    !                        !write(*,*) i,j                   
    !                        cellset(i).nearcells(2) = j
    !                        cellset(j).nearcells(4) = i
    !                        exit                            
    !                    endif
    !                end do
    !                !cellset(i).nearcells(2) = i
    !            end if
    !            if(cellset(i).nearcells(4)==0)then
    !                cellset(i).nearcells(4) = i
    !                do j = 1,ncells
    !                    dis1 = (xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(2),1))**2+(xy_coor(cellset(i).nodes(1),2)-xy_coor(cellset(j).nodes(2),2))**2
    !                    dis2 = (xy_coor(cellset(i).nodes(4),1)-xy_coor(cellset(j).nodes(3),1))**2+(xy_coor(cellset(i).nodes(4),2)-xy_coor(cellset(j).nodes(3),2))**2
    !                    if(dis1+dis2 < 1.0e-9)then
    !                        !write(*,*) i,j
    !                        cellset(i).nearcells(4) = j
    !                        cellset(j).nearcells(2) = i
    !                        exit
    !                    endif
    !                end do
    !                !cellset(i).nearcells(4) = i
    !            end if
    !        end do
    !        
    !    case(100)
    !        do i = 1,ncells
    !            if(cellset(i).nearcells(2)==0)then
    !                do j = 1,ncells
    !                    if(cellset(j).nearcells(4)==0 .and.abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))<1.0e-10&
    !                        .and.abs(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))<1.0e-10 )then
    !                        cellset(i).nearcells(2) = j
    !                        cellset(j).nearcells(4) = i
    !                        sideset(cellset(i).sides(2)).nearcells(2) = j
    !                        sideset(cellset(j).sides(4)).nearcells(2) = i
    !                    endif
    !                end do
    !            end if
    !            if(cellset(i).nearcells(1)==0)then
    !                cellset(i).nearcells(1) = i
    !            end if                
    !            if(cellset(i).nearcells(3)==0)then
    !                cellset(i).nearcells(3) = i
    !            end if
    !        end do
    !    case default
    !        write(*,*) 'Error! Occur in subroutine set_boundary'
    !    end select
     !        if(cellset(i).nearcells(1)==0)then
    !            do j = 1,ncells                 
    !                if(cellset(j).nearcells(3)==0)then
    !                    if(abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(4),1))<1.0e-10&
    !                        .and.abs(xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(3),1))<1.0e-10 )then
    !                        cellset(i).nearcells(1) = j
    !                        cellset(j).nearcells(3) = i
    !                        !write(*,*)'00',i,j
    !                    end if
    !                endif
    !            end do
    !        end if
    !        if(cellset(i).nearcells(2)==0)then
    !            do j = 1,ncells
    !                if(cellset(j).nearcells(4)==0)then
    !                        
    !                    if(abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))<1.0e-10&
    !                        .and.abs(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))<1.0e-10 )then
    !                        cellset(i).nearcells(2) = j
    !                        cellset(j).nearcells(4) = i
    !                        !sideset(cellset(i).sides(2)).nearcells(2) = j
    !                        !sideset(cellset(j).sides(4)).nearcells(2) = i
    !                        !write(*,*)'01',i,j
    !                    end if
    !                endif
    !            end do
    !        end if
    !    end do
    !    !do i = 1,ncells
    !    !    write(*,*) cellset(i).nearcells(:)
    !    !end do
    !else if(case_comp==Riemann2D_case)then
    !    do i = 1,ncells
    !        do k =1,4
    !            if(cellset(i).nearcells(k)==0)then                       
    !                cellset(i).nearcells(k) = i    
    !                    
    !            end if    
    !        end do                    
    !    end do         
    !elseif((case_comp==SodShockTube_case.OR.case_comp==LaxShockTube_case .OR. case_comp==ShuOsher_case).AND.dire_shock	== 0)then
    !    do i = 1,ncells
    !        if(cellset(i).nearcells(1)==0)then                   
    !            do j = 1,ncells
    !                    
    !                if(cellset(j).nearcells(3)==0 .and.abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(4),1))<1.0e-10&
    !                    .and.abs(xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(3),1))<1.0e-10 )then
    !                    !write(*,*) i,j
    !                    cellset(i).nearcells(1) = j
    !                    cellset(j).nearcells(3) = i
    !                endif
    !            end do
    !        end if
    !        if(cellset(i).nearcells(2)==0)then
    !            !write(*,*) i
    !            cellset(i).nearcells(2) = i
    !            !write(*,*) i,cellset(i).nearcells(2)
    !            do j = 1,ncells!�����ڲ��߽�㣬����������С���������򽻽��бߣ��ҽڵ��غϣ�������������������������x���򴫲��ļ��������⣩
    !                dis1 = (xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(1),1))**2+(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))**2
    !                dis2 = (xy_coor(cellset(i).nodes(3),1)-xy_coor(cellset(j).nodes(4),1))**2+(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))**2
    !                    
    !                if(dis1+dis2 < 1.0e-9)then                                                         
    !                    !write(*,*) i,j                   
    !                    cellset(i).nearcells(2) = j
    !                    cellset(j).nearcells(4) = i
    !                    exit                            
    !                endif
    !            end do
    !            !cellset(i).nearcells(2) = i
    !        end if
    !        if(cellset(i).nearcells(4)==0)then
    !            cellset(i).nearcells(4) = i
    !            do j = 1,ncells
    !                dis1 = (xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(2),1))**2+(xy_coor(cellset(i).nodes(1),2)-xy_coor(cellset(j).nodes(2),2))**2
    !                dis2 = (xy_coor(cellset(i).nodes(4),1)-xy_coor(cellset(j).nodes(3),1))**2+(xy_coor(cellset(i).nodes(4),2)-xy_coor(cellset(j).nodes(3),2))**2
    !                if(dis1+dis2 < 1.0e-9)then
    !                    !write(*,*) i,j
    !                    cellset(i).nearcells(4) = j
    !                    cellset(j).nearcells(2) = i
    !                    exit
    !                endif
    !            end do
    !            !cellset(i).nearcells(4) = i
    !        end if
    !    end do
    !        
    !elseif((case_comp==SodShockTube_case.OR.case_comp==LaxShockTube_case .OR. case_comp==ShuOsher_case).AND.dire_shock	== 1)then
    !    do i = 1,ncells
    !        if(cellset(i).nearcells(2)==0)then
    !            do j = 1,ncells
    !                if(cellset(j).nearcells(4)==0 .and.abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))<1.0e-10&
    !                    .and.abs(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))<1.0e-10 )then
    !                    cellset(i).nearcells(2) = j
    !                    cellset(j).nearcells(4) = i
    !                    sideset(cellset(i).sides(2)).nearcells(2) = j
    !                    sideset(cellset(j).sides(4)).nearcells(2) = i
    !                endif
    !            end do
    !        end if
    !        if(cellset(i).nearcells(1)==0)then
    !            cellset(i).nearcells(1) = i
    !        end if                
    !        if(cellset(i).nearcells(3)==0)then
    !            cellset(i).nearcells(3) = i
    !        end if
    !    end do
    !end if
    
      ! !ͨ���㴦��Jacobi����ֱ�Ӷ�����
            !if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then
            !    RC_l = nsp+1-l
            !else
            !    RC_l = l
            !end if
            !if(RC_sideth==1)then
            !    coor_C(1) = SPs(RC_l)
            !    coor_C(2) = -1.0_prec
            !else if(RC_sideth == 2)then
            !    coor_C(1) = 1.0_prec
            !    coor_C(2) = SPs(RC_l)
            !else if(RC_sideth == 3)then
            !    coor_C(1) = SPs(RC_l)
            !    coor_C(2) = 1.0_prec
            !else if(RC_sideth == 4)then
            !    coor_C(1) = -1.0_prec
            !    coor_C(2) = SPs(RC_l)
            !end if
            !
            !!ȡ�����ı��ε�Ԫ��������ֵ
            !do j = 1,4
            !    do k = 1,2
            !        quad_vertex(j,k) = xy_coor(cellset(indexCellR).nodes(j),k)
            !    end do
            !end do    
            !call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)
            !if((LC_sideth*RC_sideth.ne.3).and.(LC_sideth*RC_sideth.ne.8))then
            !    write(*,*) 'R',FPJacobi(:),'--',FPdetJ
            !end if

!subroutine revolveArray(nearCell,nearCellTemp)
!    use real_precision
!    use parameter_setting
!    implicit none 
!    integer :: l,m,j
!    real(prec),dimension(:,:,:) :: nearCell(nsp,nsp,4),nearCellTemp(nsp,nsp,4)
!
!    if(LC_sideth==1)then
!        if(cellset(nearCellIndex).nearcells(2)==index)then
!            oriR(l,m) = LaI_nPs(kesi_r,SPs,varOriR(nsp+1-l,:,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(3)==index)then
!            oriR(l,m) = LaI_nPs(kesi_r,SPs,varOriR(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(4)==index)then
!            oriR(l,m) = LaI_nPs(kesi_l,SPs,varOriR(l,:,m),nsp)
!        end if
!    elseif(LC_sideth==2)then
!        if(cellset(nearCellIndex).nearcells(1)==index)then
!            oriR(l,m) = LaI_nPs(kesi_l,SPs,varOriR(:,nsp+1-l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(3)==index)then
!            oriR(l,m) = LaI_nPs(kesi_r,SPs,varOriR(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(4)==index)then
!            oriR(l,m) = LaI_nPs(kesi_l,SPs,varOriR(l,:,m),nsp)
!        end if        
!    elseif(LC_sideth==3)then
!        if(cellset(nearCellIndex).nearcells(1)==index)then
!            oriR(l,m) = LaI_nPs(kesi_l,SPs,varOriR(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(2)==index)then
!            oriR(l,m) = LaI_nPs(kesi_r,SPs,varOriR(l,:,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(4)==index)then
!            oriR(l,m) = LaI_nPs(kesi_l,SPs,varOriR(nsp+1-l,:,m),nsp)
!        end if          
!    elseif(LC_sideth==4)then
!        if(cellset(nearCellIndex).nearcells(1)==index)then
!            oriR(l,m) = LaI_nPs(kesi_l,SPs,varOriR(:,l,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(2)==index)then
!            oriR(l,m) = LaI_nPs(kesi_r,SPs,varOriR(l,:,m),nsp)
!        else if(cellset(nearCellIndex).nearcells(3)==index)then
!            oriR(l,m) = LaI_nPs(kesi_r,SPs,varOriR(:,nsp+1-l,m),nsp)
!        end if         
!    end if
!end subroutine revolveArray

!subroutine set_BC 
!    !�����������ñ߽磬��δ�����������������ã�Ŀǰ����ض�����
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer i,j,k,l
!    integer :: index1 = 1,index2 = 1,index3 = 1,index4 = 1,index0 = 1,index_bd1,index_bd2
!    real(prec) :: dis1,dis2,d_D,d_R,d_U,d_L
!    !integer :: n_BCells_D = 0,n_BCells_R = 0,n_BCells_U = 0,n_BCells_L = 0
!    integer :: i_D,i_R,i_U,i_L,O_sideth
!    
!    !do i = 1,nnodes
!    !    d_D = abs(xy_coor(i,2)-yl)!Down
!    !    d_R = abs(xy_coor(i,1)-xr)!Right
!    !    d_U = abs(xy_coor(i,2)-yr)!Up
!    !    d_L = abs(xy_coor(i,1)-xl)!Left
!    !    if(d_D < 1.0e-13)then
!    !        n_BCells_D = n_BCells_D + 1 !�²�߽ڵ������
!    !    end if
!    !    if(d_R < 1.0e-13)then
!    !        n_BCells_R = n_BCells_R + 1 !�Ҳ�߽ڵ������
!    !    end if
!    !    if(d_U < 1.0e-13)then
!    !        n_BCells_U = n_BCells_U + 1 !�ϲ�߽ڵ������
!    !    end if
!    !    if(d_L < 1.0e-13)then
!    !        n_BCells_L = n_BCells_L + 1 !���߽ڵ������
!    !    end if
!    !end do
!    !
!    !n_BCells_D = n_BCells_D - 1 !�ߵ�����
!    !n_BCells_R = n_BCells_R - 1 
!    !n_BCells_U = n_BCells_U - 1 
!    !n_BCells_L = n_BCells_L - 1 
!    !!write(*,*) n_BCells_D, n_BCells_R, n_BCells_U, n_BCells_L
!    !nbdsides = n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
!    !allocate(BoundCells_index_set(nbdsides))
!    !do i = 1,ncells    
!    !    i_D = 0
!    !    do j =1,4
!    !        d_D = abs(xy_coor(cellset(i).nodes(j),2)-yl)!Down
!    !        if(d_D < 1.0e-13)then
!    !            i_D = i_D+1
!    !        end if
!    !    end do
!    !    if(i_D==2)then   !�ױ��ϵĵ�Ԫ
!    !        O_sideth = 1 !�߽絥Ԫ�ڵױ��ϵĲ�߱��һ��Ϊ1
!    !        if(cellset(i).nearcells(1)/=0)then  !�Ų��Ų�Ϊ1�ĵ�Ԫ
!    !            do j = 1,4
!    !                if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ1
!    !                    O_sideth = j 
!    !                end if
!    !            end do
!    !        end if
!    !        cellset(i).nearcells(O_sideth) = ncells+index1!�����Ϊncells+index1�Ľṹ�帳Ϊi���ڵ�Ԫ
!    !        cellset(ncells+index1).nearcells(3) = i
!    !        BoundCells_index_set(index1)=i      !��¼�߽絥Ԫ��ŵ�����
!    !        index1 = index1+1
!    !    end if
!    !    i_R = 0
!    !    do j =1,4
!    !        d_R = abs(xy_coor(cellset(i).nodes(j),1)-xr)!Right
!    !        if(d_R < 1.0e-13)then
!    !            i_R = i_R+1
!    !        end if
!    !    end do
!    !    if(i_R==2)then   !�ױ��ϵĵ�Ԫ
!    !        O_sideth = 2 !�߽絥Ԫ���ұ��ϵĲ�߱��һ��Ϊ2
!    !        if(cellset(i).nearcells(2)/=0)then  !�Ų��Ų�Ϊ2�ĵ�Ԫ
!    !            do j = 1,4
!    !                if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ2
!    !                    O_sideth = j 
!    !                end if
!    !            end do
!    !        end if
!    !        cellset(i).nearcells(O_sideth) = ncells+n_BCells_D+index2
!    !        cellset(ncells+n_BCells_D+index2).nearcells(4) = i
!    !        BoundCells_index_set(n_BCells_D+index2)=i
!    !        index2 = index2+1
!    !    end if
!    !    i_U = 0
!    !    do j =1,4
!    !        d_U = abs(xy_coor(cellset(i).nodes(j),2)-yr)!Up
!    !        if(d_U < 1.0e-13)then
!    !            i_U = i_U+1
!    !        end if
!    !    end do
!    !    if(i_U==2)then   !�ױ��ϵĵ�Ԫ
!    !        O_sideth = 3 !�߽絥Ԫ���ϱ��ϵĲ�߱��һ��Ϊ3
!    !        if(cellset(i).nearcells(2)/=0)then  !�Ų��Ų�Ϊ3�ĵ�Ԫ
!    !            do j = 1,4
!    !                if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ3
!    !                    O_sideth = j 
!    !                end if
!    !            end do
!    !        end if
!    !        cellset(i).nearcells(O_sideth) = ncells+n_BCells_D+n_BCells_R+index3
!    !        cellset(ncells+n_BCells_D+n_BCells_R+index3).nearcells(1) = i
!    !        BoundCells_index_set(n_BCells_D+n_BCells_R+index3)=i
!    !        index3 = index3+1
!    !    end if
!    !    i_L = 0
!    !    do j =1,4
!    !        d_L = abs(xy_coor(cellset(i).nodes(j),1)-xl)!Up
!    !        if(d_L < 1.0e-13)then
!    !            i_L = i_L+1
!    !        end if
!    !    end do
!    !    if(i_L==2)then   !�ױ��ϵĵ�Ԫ
!    !        O_sideth = 4 !�߽絥Ԫ������ϵĲ�߱��һ��Ϊ4
!    !        if(cellset(i).nearcells(4)/=0)then  !�Ų��Ų�Ϊ4�ĵ�Ԫ
!    !            do j = 1,4
!    !                if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ4
!    !                    O_sideth = j 
!    !                end if
!    !            end do
!    !        end if
!    !        cellset(i).nearcells(O_sideth) = ncells+n_BCells_D+n_BCells_R+n_BCells_U+index4
!    !        cellset(ncells+n_BCells_D+n_BCells_R+n_BCells_U+index4).nearcells(2) = i
!    !        BoundCells_index_set(n_BCells_D+n_BCells_R+n_BCells_U+index4)=i
!    !        index4 = index4+1
!    !    end if
!    !   
!    !end do
!    !index0 = index0-1       !�߽����
!    !!write(*,*)index1,index2,index3,index4,index0
!    !!write(*,*)BoundCells_index_set
!    !!do i = ncells+1,ncells+nbdsides
!    !!    write(*,*)cellset(i).nearcells(:)    
!    !!end do
!    !
!    !stop
!
!    
!    !�����������ȷ���߽�����equEntropy_case = 0, SodShockTube_case = 1, LaxShockTube_case = 2,ShuOsher_case = 3,Riemann2D_case = 4
!    if(case_comp==equEntropy_case)then            
!        do i = 1,ncells
!            if(cellset(i).nearcells(1)==0)then
!                do j = 1,ncells                 
!                    if(cellset(j).nearcells(3)==0)then
!                        if(abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(4),1))<1.0e-10&
!                            .and.abs(xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(3),1))<1.0e-10 )then
!                            cellset(i).nearcells(1) = j
!                            cellset(j).nearcells(3) = i
!                            !write(*,*)'00',i,j
!                        end if
!                    endif
!                end do
!            end if
!            if(cellset(i).nearcells(2)==0)then
!                do j = 1,ncells
!                    if(cellset(j).nearcells(4)==0)then
!                            
!                        if(abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))<1.0e-10&
!                            .and.abs(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))<1.0e-10 )then
!                            cellset(i).nearcells(2) = j
!                            cellset(j).nearcells(4) = i
!                            !sideset(cellset(i).sides(2)).nearcells(2) = j
!                            !sideset(cellset(j).sides(4)).nearcells(2) = i
!                            !write(*,*)'01',i,j
!                        end if
!                    endif
!                end do
!            end if
!        end do
!        !do i = 1,ncells
!        !    write(*,*) cellset(i).nearcells(:)
!        !end do
!    else if(case_comp==Riemann2D_case .OR. case_comp==DoubleMach_case)then
!        do i = 1,ncells
!            do k =1,4
!                if(cellset(i).nearcells(k)==0)then                       
!                    cellset(i).nearcells(k) = i                           
!                end if    
!            end do                    
!        end do         
!    elseif((case_comp==SodShockTube_case.OR.case_comp==LaxShockTube_case .OR. case_comp==ShuOsher_case).AND.dire_shock	== 0)then
!        do i = 1,ncells
!            if(cellset(i).nearcells(1)==0)then                   
!                do j = 1,ncells
!                        
!                    if(cellset(j).nearcells(3)==0 .and.abs(xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(4),1))<1.0e-10&
!                        .and.abs(xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(3),1))<1.0e-10 )then
!                        !write(*,*) i,j
!                        cellset(i).nearcells(1) = j
!                        cellset(j).nearcells(3) = i
!                    endif
!                end do
!            end if
!            if(cellset(i).nearcells(2)==0)then
!                !write(*,*) i
!                cellset(i).nearcells(2) = i
!                !write(*,*) i,cellset(i).nearcells(2)
!                do j = 1,ncells!�����ڲ��߽�㣬����������С���������򽻽��бߣ��ҽڵ��غϣ�������������������������x���򴫲��ļ��������⣩
!                    dis1 = (xy_coor(cellset(i).nodes(2),1)-xy_coor(cellset(j).nodes(1),1))**2+(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))**2
!                    dis2 = (xy_coor(cellset(i).nodes(3),1)-xy_coor(cellset(j).nodes(4),1))**2+(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))**2
!                        
!                    if(dis1+dis2 < 1.0e-9)then                                                         
!                        !write(*,*) i,j                   
!                        cellset(i).nearcells(2) = j
!                        cellset(j).nearcells(4) = i
!                        exit                            
!                    endif
!                end do
!                !cellset(i).nearcells(2) = i
!            end if
!            if(cellset(i).nearcells(4)==0)then
!                cellset(i).nearcells(4) = i
!                do j = 1,ncells
!                    dis1 = (xy_coor(cellset(i).nodes(1),1)-xy_coor(cellset(j).nodes(2),1))**2+(xy_coor(cellset(i).nodes(1),2)-xy_coor(cellset(j).nodes(2),2))**2
!                    dis2 = (xy_coor(cellset(i).nodes(4),1)-xy_coor(cellset(j).nodes(3),1))**2+(xy_coor(cellset(i).nodes(4),2)-xy_coor(cellset(j).nodes(3),2))**2
!                    if(dis1+dis2 < 1.0e-9)then
!                        !write(*,*) i,j
!                        cellset(i).nearcells(4) = j
!                        cellset(j).nearcells(2) = i
!                        exit
!                    endif
!                end do
!                !cellset(i).nearcells(4) = i
!            end if
!        end do
!            
!    elseif((case_comp==SodShockTube_case.OR.case_comp==LaxShockTube_case .OR. case_comp==ShuOsher_case).AND.dire_shock	== 1)then
!        do i = 1,ncells
!            if(cellset(i).nearcells(2)==0)then
!                do j = 1,ncells
!                    if(cellset(j).nearcells(4)==0 .and.abs(xy_coor(cellset(i).nodes(2),2)-xy_coor(cellset(j).nodes(1),2))<1.0e-10&
!                        .and.abs(xy_coor(cellset(i).nodes(3),2)-xy_coor(cellset(j).nodes(4),2))<1.0e-10 )then
!                        cellset(i).nearcells(2) = j
!                        cellset(j).nearcells(4) = i
!                        sideset(cellset(i).sides(2)).nearcells(2) = j
!                        sideset(cellset(j).sides(4)).nearcells(2) = i
!                    endif
!                end do
!            end if
!            if(cellset(i).nearcells(1)==0)then
!                cellset(i).nearcells(1) = i
!            end if                
!            if(cellset(i).nearcells(3)==0)then
!                cellset(i).nearcells(3) = i
!            end if
!        end do
!    end if
!           
!   !do i  = 1,ncells
!   ! write(*,*) i,cellset(i).nearcells(:)
!   !end do
!end subroutine set_BC 
                    !if(k==1)nextNearSideth = 4
                    !indexNearCell = cellset(j).nearcells(nextNearSideth)  
                    !call search_sameVertexNearCells(j,i,indexNearCell,cells_contain_vertex)      
!RECURSIVE subroutine search_sameVertexNearCells(index_start,indexVertex,indexNearCell,array)
!    !�ݹ��ѯ
!    use real_precision
!    use type_module
!    use global_var
!    implicit none
!    integer :: indexVertex,index_start  !������
!    integer :: index,i,j,k,l,m,indexNearCell,nextNearSideth
!    integer :: array(10)
!    integer :: count
!    !write(*,*)count,indexNearCell
!    count = 2
!    do j = 1,4
!        if(cellset(indexNearCell).nodes(j)==indexVertex)then
!            array(count) = indexNearCell 
!            count = count+1
!            nextNearSideth = j-1
!            if(j==1)nextNearSideth = 4
!            index = indexNearCell
!            indexNearCell = cellset(index).nearcells(nextNearSideth)       
!        end if
!    end do 
!   
!    if(indexNearCell == index_start .OR. indexNearCell > ncells)then
!        array(1) = index_start
!    else       
!        call search_sameVertexNearCells(index_start,indexVertex,indexNearCell,array)
!    end if
!    
!    
!end subroutine search_sameVertexNearCells 
        !do l = 1,4
        !    if(cellset(startCell).nearcells(l) > ncells)then    !ͳ���Ƿ��Ǳ߽��ϵĵ�Ԫ��count1>0���ǡ�����߽絥Ԫ�ı����ncell֮��
        !        count1 = count1 + 1 
        !    end if
        !end do

    !do i = 1,ncells
    !    do j = 1,nsp
    !        do k = 1,nsp
    !            write(20,'(2F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2)
    !        end do
    !    end do
    !end do     
    !cycle1:do i = 1,ncells
    !    cycle2:do j = 1,nsp-1
    !        !if(j==nsp)cycle
    !        cycle3:do k = 1,nsp-1  
    !            P1 = (i-1)*nsp**2+(j-1)*nsp+k
    !            P2 = P1+1
    !            P3 = P2+nsp
    !            P4 = P3-1
    !            write(20,*) P1,P2,P3,P4
    !        end do cycle3
    !    end do cycle2
    !end do cycle1
    !!do i = 1,nsides
    !!    write(*,*)sideset(i).nearcells
    !!end do
    !
    !!�����Ǹ��ݵ�Ԫ��ߣ�������ԵĽ�㣬�����ı���
    !allocate(SC_BP_index(4,nsp))
    !do j = 1,nsp!�����ң����µ���
    !    SC_BP_index(1,j) = j
    !    SC_BP_index(2,j) = j*nsp
    !    SC_BP_index(3,j) = (nsp-1)*nsp+j
    !    SC_BP_index(4,j) = (j-1)*nsp+1
    !end do
    !!�����ǿ�����Ԫ����Ľ����
    !Vertex_P_index(1) = 1
    !Vertex_P_index(2) = nsp
    !Vertex_P_index(3) = nsp*nsp
    !Vertex_P_index(4) = (nsp-1)*nsp+1
    !do i = 1,nsides
    !    iL = sideset(i).nearcells(1)
    !    iR = sideset(i).nearcells(2)    
    !    if(iL==0 .OR. iR == 0)then
    !        cycle
    !    else
    !        do j = 1,4
    !            if(cellset(iL).nearcells(j)==iR)then
    !                LC_sideth = j !Left cell Side -th.�˱��Ǹ���Ԫ��j�����
    !            end if
    !            if(cellset(iR).nearcells(j)==iL)then
    !                RC_sideth = j
    !            end if
    !        end do
    !        do k = 1,nsp-1
    !            P1 = SC_BP_index(LC_sideth,k)+(iL-1)*nsp*nsp
    !            P2 = SC_BP_index(LC_sideth,k+1)+(iL-1)*nsp*nsp
    !            if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then
    !                rk = nsp+1-k
    !                P3 = SC_BP_index(RC_sideth,rk-1)+(iR-1)*nsp*nsp
    !                P4 = SC_BP_index(RC_sideth,rk)+(iR-1)*nsp*nsp
    !            else
    !                P3 = SC_BP_index(RC_sideth,k+1)+(iR-1)*nsp*nsp
    !                P4 = SC_BP_index(RC_sideth,k)+(iR-1)*nsp*nsp
    !            end if
    !            write(20,*) P1,P2,P3,P4
    !        end do     
    !    end if
    !end do
    !!�����ǵ�Ԫ���㴦�����ӡ����˴��е㸴�ӣ����ܴ��������Σ��ı��Σ�����Σ������Σ�����������Ρ������ʱ��Ҫ����հף�ʹ�����м䲻���ס�
    !!   ���ݶ�����ң�Ȼ��ȡ��������Χ��Ԫ������Ľ�㡣���ӳ��ı���+�����ε���ʽ
    !!write(*,*) ncells,nbdsides  
    !
    !sum_vertexCell = 0
    !Loop1:do i = 1,nnodes    
    !    cells_contain_vertex = 0
    !    Loop2:do j = 1,ncells
    !        Loop3:do k = 1,4
    !            if(cellset(j).nodes(k)==i)then!��ѯ�������˽ڵ�ĵ�һ����Ԫ�����Ϊj
    !                startCell = j
    !                exit Loop2
    !            end if
    !        end do Loop3
    !    end do Loop2
    !    !write(*,*)'k',k
    !    !֮�����ڵ�Ԫ��
    !    d_D = abs(xy_coor(i,2)-yl)!Down
    !    d_R = abs(xy_coor(i,1)-xr)!Right
    !    d_U = abs(xy_coor(i,2)-yr)!Up
    !    d_L = abs(xy_coor(i,1)-xl)!Left
    !    
    !    if(d_D < 1.0e-13 .OR. d_R < 1.0e-13 .OR. d_U < 1.0e-13 .OR. d_L < 1.0e-13)then
    !        count1 = 1
    !    else
    !        count1 = 0
    !    end if
    !
    !    if(count1 == 1)then     !count1 = 1��˵���˶�������ڼ�����߽���
    !        !��������ѯ�������˶�������е�Ԫ��������֪�ĵ�Ԫ���Դ������ڵ�Ԫ���֣���С��ѯ��Χ  
    !        !����߶���Ҫ��ѯ����Ϊ���߽�ֹͣ����ȷ����һ����Ԫ����λ��
    !        !��һ���
    !        nextNearSideth = k-1
    !        if(nextNearSideth == 0)nextNearSideth=4
    !        nextCell = cellset(startCell).nearcells(nextNearSideth)
    !        count2 = 1                       
    !        cells_contain_vertex(count2,1) = startCell!��¼��һ����Ԫ
    !        cells_contain_vertex(count2,2) = k        !��¼�ڼ�������
    !        do while((nextCell < ncells+1) .and. (nextCell > 0))   
    !  
    !            do m = 1,4
    !                if(cellset(nextCell).nodes(m)==i) vertex_th = m!mΪ�ö�����nextCell�еı��
    !            end do
    !
    !            nextNearSideth = vertex_th-1
    !            if(nextNearSideth == 0)nextNearSideth=4
    !            if(cellset(nextCell).nearcells(nextNearSideth)==startCell)then!�ų�startCell���nextCell���ڵ�Ԫ
    !                nextNearSideth = vertex_th
    !                startCell = nextCell
    !                nextCell = cellset(startCell).nearcells(nextNearSideth)                    
    !            elseif(cellset(nextCell).nearcells(vertex_th)==startCell)then
    !                startCell = nextCell
    !                nextCell = cellset(startCell).nearcells(nextNearSideth)
    !            else
    !                nextCell = 0
    !            end if                          
    !            count2 = count2+1
    !            cells_contain_vertex(count2,1) = startCell!��¼��Ԫ
    !            cells_contain_vertex(count2,2) = vertex_th
    !        end do
    !        !�ڶ����
    !        nextNearSideth = k
    !        startCell = j
    !        nextCell = cellset(startCell).nearcells(nextNearSideth)
    !        do while((nextCell < ncells+1) .and. (nextCell > 0))   
    !            !write(*,*)nextCell
    !            
    !            do m = 1,4
    !                if(cellset(nextCell).nodes(m)==i)vertex_th = m!mΪ�ö�����nextCell�еı��
    !            end do
    !            nextNearSideth = vertex_th-1
    !            if(nextNearSideth == 0)nextNearSideth=4
    !            if(cellset(nextCell).nearcells(nextNearSideth)==startCell)then!�ų�startCell���nextCell���ڵ�Ԫ
    !                nextNearSideth = vertex_th
    !                startCell = nextCell
    !                nextCell = cellset(startCell ).nearcells(nextNearSideth)
    !            elseif(cellset(nextCell).nearcells(vertex_th)==startCell)then
    !                startCell = nextCell
    !                nextCell = cellset(startCell).nearcells(nextNearSideth)
    !            else
    !                nextCell = 0
    !            end if                        
    !            count2 = count2+1
    !            cells_contain_vertex(count2,1) = startCell!��¼��Ԫ
    !            cells_contain_vertex(count2,2) = vertex_th
    !        end do
    !    elseif(count1 == 0)then     !�������ڲ��Ķ��㣬���ܰ��������������ϵ�Ԫ��
    !        !k�Ǹö����ڵ�һ����Ԫ�ı�š��ɴ˿���֪����k-1,k��������Ե����ڵ�ԪҲ�����˶��㡣k=1ʱ��k-1��Ϊ4
    !        !����whileѭ��һֱ��ѯ�ڵ�Ԫ��ֱ��ָ���һ����Ԫ       
    !        nextNearSideth = k-1
    !        if(nextNearSideth == 0)nextNearSideth=4
    !        nextCell = cellset(startCell).nearcells(nextNearSideth)
    !        count2 = 1                      
    !        cells_contain_vertex(count2,1) = startCell!��¼��һ����Ԫ
    !        cells_contain_vertex(count2,2) = k
    !        do while(nextCell /= j)                           
    !
    !            do m = 1,4
    !                if(cellset(nextCell).nodes(m)==i)vertex_th = m!mΪ�ö�����nextCell�еı��
    !            end do
    !
    !            nextNearSideth = vertex_th-1
    !            if(nextNearSideth == 0)nextNearSideth=4
    !            if(cellset(nextCell).nearcells(nextNearSideth)==startCell)then!�ų�startCell���nextCell���ڵ�Ԫ
    !                nextNearSideth = vertex_th
    !                startCell = nextCell
    !                nextCell = cellset(startCell).nearcells(nextNearSideth)
    !            elseif(cellset(nextCell).nearcells(vertex_th)==startCell)then
    !                startCell = nextCell
    !                nextCell = cellset(startCell).nearcells(nextNearSideth)
    !            end if
    !            count2 = count2+1
    !            cells_contain_vertex(count2,1) = startCell
    !            cells_contain_vertex(count2,2) = vertex_th
    !        end do
    !    end if                            
    !    !write(*,*)i,'-----------------------------------------------'
    !    !write(*,*)cells_contain_vertex(1:5,1)
    !    
    !    count3 = 0
    !    do j = 1,10
    !        if(cells_contain_vertex(j,1)/=0)then
    !            count3 = count3+1   !ͳ�ư����ö��㵥Ԫ��
    !        end if
    !    end do
    !    if(count3 == 3)then  !count3 <= 2�Ǳ߽��ϵĵ�Ԫ��ֻ��������Ԫ�����û����ı���
    !        P1 = (cells_contain_vertex(1,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1,2))
    !        P2 = (cells_contain_vertex(2,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(2,2))
    !        P3 = (cells_contain_vertex(3,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(3,2))
    !        P4 = (cells_contain_vertex(1,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1,2))
    !        write(20,*) P1,P2,P3,P4
    !        sum_vertexCell = sum_vertexCell + 1
    !    elseif(count3 == 4)then !�ĸ���Ԫ
    !        P1 = (cells_contain_vertex(1,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1,2))
    !        P2 = (cells_contain_vertex(2,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(2,2))
    !        P3 = (cells_contain_vertex(3,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(3,2))
    !        P4 = (cells_contain_vertex(4,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(4,2))
    !        write(20,*) P1,P2,P3,P4
    !        sum_vertexCell = sum_vertexCell + 1
    !    else                    
    !        !4�����ϵ�Ԫ.�����������:  p = count3 |  �ı���    ������
    !        !                           pΪ����    |  p/2 - 1      1                               
    !        !                           pΪż��    | (p-1)/2 - 1   0    
    !        !ż������������������д��һ������ʽ
    !        count4 = (count3+1)/2 - 1   !��ȡ��
    !        do j = 1,count4
    !            P1 = (cells_contain_vertex(1+(j-1),1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(1+(j-1),2))
    !            P2 = (cells_contain_vertex(2+(j-1),1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(2+(j-1),2))
    !            P3 = (cells_contain_vertex(count3-j,1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(count3-j,2))
    !            P4 = (cells_contain_vertex(count3-(j-1),1)-1)*nsp*nsp + Vertex_P_index(cells_contain_vertex(count3-(j-1),2))
    !            write(20,*) P1,P2,P3,P4
    !            sum_vertexCell = sum_vertexCell + 1
    !        end do     
    !    end if
    !end do Loop1
    !!write(*,*)sum_vertexCell
    !!do i =1,ncells
    !!    write(*,*)i,cellset(i).nodes
    !!end do
    !
    !!write(*,*)((SC_BP_index(i,j),j=1,nsp),i=1,4)

 !    cycle1:do i = 1,ncells
    !        cycle2:do j = 1,nsp-1
    !            !if(j==nsp)cycle
    !            cycle3:do k = 1,nsp-1  
    !                P1 = (i-1)*nsp**2+(j-1)*nsp+k
    !                P2 = P1+1
    !                P3 = P2+nsp
    !                P4 = P3-1
    !                write(40,*) P1,P2,P3,P4
    !            end do cycle3
    !        end do cycle2
    !    end do cycle1
!
    !allocate(SC_BP_index(4,nsp))
    !do j = 1,nsp!�����ң����µ���
    !    SC_BP_index(1,j) = j
    !    SC_BP_index(2,j) = j*nsp
    !    SC_BP_index(3,j) = (nsp-1)*nsp+j
    !    SC_BP_index(4,j) = (j-1)*nsp+1
    !end do
    !do i = 1,nsides
    !    iL = sideset(i).nearcells(1)
    !    iR = sideset(i).nearcells(2)    
    !    if(iL==0 .OR. iR == 0)then
    !        cycle
    !    else
    !        do j = 1,4
    !            if(cellset(iL).nearcells(j)==iR)then
    !                LC_sideth = j !Left cell Side -th.�˱��Ǹ���Ԫ��j�����
    !            end if
    !            if(cellset(iR).nearcells(j)==iL)then
    !                RC_sideth = j
    !            end if
    !        end do
    !        do k = 1,nsp-1
    !            P1 = SC_BP_index(LC_sideth,k)+(iL-1)*nsp*nsp
    !            P2 = SC_BP_index(LC_sideth,k+1)+(iL-1)*nsp*nsp
    !            if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then
    !                rk = nsp+1-k
    !                P3 = SC_BP_index(RC_sideth,rk-1)+(iR-1)*nsp*nsp
    !                P4 = SC_BP_index(RC_sideth,rk)+(iR-1)*nsp*nsp
    !            else
    !                P3 = SC_BP_index(RC_sideth,k+1)+(iR-1)*nsp*nsp
    !                P4 = SC_BP_index(RC_sideth,k)+(iR-1)*nsp*nsp
    !            end if
    !            write(40,*) P1,P2,P3,P4
    !        end do
    !        
    !    end if
    !end do

        !write(40,*) 'ZONE, N=',ncells*nsp*nsp,', E=',ncells*(nsp-1)*(nsp-1)+(nsides-nbdsides)*(nsp-1),', F=FEPOINT, ET=QUADRILATERAL'
        !do i = 1,ncells
        !    do j = 1,nsp
        !        do k = 1,nsp
        !            write(40,'(6F20.10)')cellset(i).sp_coor(j,k,1),cellset(i).sp_coor(j,k,2),cellset(i).spvalue_ori(j,k,1),cellset(i).spvalue_ori(j,k,2),cellset(i).spvalue_ori(j,k,3),cellset(i).spvalue_ori(j,k,4)
        !        end do
        !    end do
        !end do     

 !do j = 1,4
    !    if(cellset(nearCellIndex).nearcells(j)==index)then
    !        RC_sideth = j
    !    end if
    !end do
    !
    !
    !if(index == sideset(sideIndex).nearCells(1))then
    !    FP_upw_temp(:) = sideset(sideIndex).fpvalue_upw(l,:) 
    !    !cellset(index).fluxF_innerfp(:,1,:) = sideset(cellset(index).sides(4)).fpvalue_upw(:,:)       !���ǵ�Ԫ�ڱߵ������ͨ����˳��һ�� 
    !elseif(index == sideset(sideIndex).nearCells(2))then
    !
    !    if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12))then!1-2 3-4 ���ڣ�ͨ����˳���෴���߽絥Ԫ���Ǳ߽�ߵ���Ԫ�����ÿ���
    !        FP_upw_temp(:)=sideset(sideIndex).fpvalue_upw(nsp+1-l,:)      !ͨ����˳���෴         
    !    else
    !        FP_upw_temp(:)=sideset(sideIndex).fpvalue_upw(l,:)      !ͨ����˳��һ�� 
    !    end if
    !end if
    !
!subroutine get_oriL(index,nearCellIndex,LC_sideth,RC_sideth,l,oriL)
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none 
!    integer ::index,nearCellIndex,LC_sideth,RC_sideth,l,m,k
!    real(prec),dimension(:,:) :: u(3,4),conL(4),ll(4,4),rr(4,4),chara_con(3,4)
!    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4),varConL(nsp,nsp,4),varConR(nsp,nsp,4)
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),oriL(4),chara_oriL(4),chara_conL(4)
!    real(prec) :: realDis1,realDis2
!    
!    !ȡ�����ҵ�Ԫ�ϵ�ԭʼ����   
!    varOriL(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!    varOriR(:,:,:) = cellset(nearCellIndex).spvalue_ori(:,:,:)
!    !ȡ�����ҵ�Ԫ�ϵ��غ���� 
!    varConL(:,:,:) = cellset(index).spvalue_con(:,:,:)
!    varConR(:,:,:) = cellset(nearCellIndex).spvalue_con(:,:,:)
!    if(cellset(index).Beta == 0)then
!        do m = 1,4
!            if(LC_sideth == 1)then                                      
!                oriL(m) = LaI_nPs(kesi_l,SPs,varOriL(:,l,m),nsp)
!            else if(LC_sideth == 2)then
!                oriL(m) = LaI_nPs(kesi_r,SPs,varOriL(l,:,m),nsp)
!            else if(LC_sideth == 3)then
!                oriL(m) = LaI_nPs(kesi_r,SPs,varOriL(:,l,m),nsp)
!            else if(LC_sideth == 4)then
!                oriL(m) = LaI_nPs(kesi_l,SPs,varOriL(l,:,m),nsp)
!            end if
!        end do
!    else if(cellset(index).Beta == 1)then 
!        !write(*,*)'23',LC_sideth
!        !if(index == 160)then
!        !    write(*,*)'001',index,nearCellIndex,LC_sideth,RC_sideth,l
!        !    write(*,*)cellset(160).spvalue_ori(:,:,:)
!        !    pause
!        !end if
!        select case(var_type)
!        case(ori_type)
!            !!----ԭʼ����---------------------
!            if(LC_sideth == 1)then    
!                u(2,:) = varOriL(1,l,:)
!                u(3,:) = varOriL(2,l,:)
!            else if(LC_sideth == 2)then
!                u(2,:) = varOriL(l,nsp,:)
!                u(3,:) = varOriL(l,nsp-1,:)
!            else if(LC_sideth == 3)then
!                u(2,:) = varOriL(nsp,l,:)
!                u(3,:) = varOriL(nsp-1,l,:)
!            else if(LC_sideth == 4)then
!                u(2,:) = varOriL(l,1,:)
!                u(3,:) = varOriL(l,2,:)
!            end if     
!            call get_nearCellInfo(index,nearCellIndex,LC_sideth,RC_sideth,l,varOriR,u(:,:),realDis1)
!            call FaceFluxC2NNW2(index,nearCellIndex,LC_sideth,l,u,realDis1,oriL(:),innerDisOri)
!            !write(*,*)'oriLR',oriL(1),innerDisOri(1)
!            !!----ԭʼ����---------------------
!        case(con_type)
!            !!----�غ����---------------------
!            if(LC_sideth == 1)then    
!                u(2,:) = varConL(1,l,:)
!                u(3,:) = varConL(2,l,:)
!            else if(LC_sideth == 2)then
!                u(2,:) = varConL(l,nsp,:)
!                u(3,:) = varConL(l,nsp-1,:)
!            else if(LC_sideth == 3)then
!                u(2,:) = varConL(nsp,l,:)
!                u(3,:) = varConL(nsp-1,l,:)
!            else if(LC_sideth == 4)then
!                u(2,:) = varConL(l,1,:)
!                u(3,:) = varConL(l,2,:)
!            end if        
!            call get_nearCellInfo(index,nearCellIndex,LC_sideth,RC_sideth,l,varConR,u(:,:),realDis1)
!            call FaceFluxC2NNW2(index,nearCellIndex,LC_sideth,l,u,realDis1,conL(:),innerDisCon)
! 
!            call Func_con_to_ori(conL(:),oriL(:))
!            call Func_con_to_ori(innerDisCon(:),innerDisOri(:))
!            !write(*,*)'oriLR',oriL(1),innerDisOri(1)
!            !!----�غ����---------------------
!        case(character_type)
!            !!----��������---------------------
!            if(LC_sideth == 1)then    
!                u(2,:) = varOriL(1,l,:)
!                u(3,:) = varOriL(2,l,:)
!            else if(LC_sideth == 2)then
!                u(2,:) = varOriL(l,nsp,:)
!                u(3,:) = varOriL(l,nsp-1,:)
!            else if(LC_sideth == 3)then
!                u(2,:) = varOriL(nsp,l,:)
!                u(3,:) = varOriL(nsp-1,l,:)
!                !if(index==625)then
!                !    write(*,*)nsp-1,l,u(3,:)
!                !end if
!            else if(LC_sideth == 4)then
!                u(2,:) = varOriL(l,1,:)
!                u(3,:) = varOriL(l,2,:)
!            end if        
!            !write(*,*)'LL',l
!            call get_nearCellInfo(index,nearCellIndex,LC_sideth,RC_sideth,l,varOriR,u(:,:),realDis1)            
!            call proj_matrix(u(2,:),ll,rr,LC_sideth)  !(uԭʼ����)
!            !write(*,*)'u',u(2,1)
!            do k = 1,3
!                call Characteristic_projection(u(k,:),ll,chara_con(k,:)) !(chara_con�����غ����)
!                !write(*,"(4F10.5)")chara_u(k,:)
!            end do
!            call FaceFluxC2NNW2(index,nearCellIndex,LC_sideth,l,chara_con,realDis1,chara_conL,chara_innerDisCon)
!            call Inverse_Characteristic_projection(chara_conL(:),rr,oriL(:))
!            call Inverse_Characteristic_projection(chara_innerDisCon(:),rr,innerDisOri(:))
!            !write(*,*)'ch',oriL(1),innerDisOri(1)
!        end select
!        !write(*,*)'24'
!        if(method_subcell == method_NNW)then
!            if(LC_sideth==1)then    !�ӵ�Ԫ�Ȱѽ���ԭʼ������������֮��絥Ԫ��������һ��
!                cellset(index).fluxG_innerfp(l,2,:) = innerDisOri(:)         
!            elseif(LC_sideth==2)then
!                cellset(index).fluxF_innerfp(l,nsp,:) = innerDisOri(:)
!            elseif(LC_sideth==3)then
!                cellset(index).fluxG_innerfp(l,nsp,:) = innerDisOri(:)
!            elseif(LC_sideth==4)then
!                cellset(index).fluxF_innerfp(l,2,:) = innerDisOri(:)
!            end if
!        elseif(method_subcell == method_Godunov)then
!            oriL(:) = u(2,:)
!        end if
!        
!        !if(index==48 )then
!        !    write(*,*)nt_temp
!        !    write(*,*)'1111111',LC_sideth,l,index,nearCellIndex,realDis1
!        !    !write(*,*)'u2',u(2,1)
!        !    write(*,*)'u(:,1)',u(:,1)
!        !    !write(*,*)'ll',ll
!        !    !write(*,*)'rr',rr
!        !    !write(*,*)'oril',oriL(:)
!        !    !write(*,*)'chara_conL(:)',chara_conL(:)
!        !    !write(*,*)'chara_innerDisCon(:)',chara_innerDisCon(:)
!        !    !write(*,*)'innerDisOri',innerDisOri(:)
!        !    !write(*,*)cellset(index).Mdirect
!        !    !if(LC_sideth==1)then
!        !    !    write(*,*)'gg',index,LC_sideth,innerDisOri(1),cellset(index).fluxG_innerfp(l,2,1)
!        !    !elseif(LC_sideth==2)then
!        !    !    write(*,*)'ff',index,LC_sideth,innerDisOri(1),cellset(index).fluxF_innerfp(l,nsp,1)
!        !    !elseif(LC_sideth==3)then
!        !    !    write(*,*)'gg',index,LC_sideth,innerDisOri(1),cellset(index).fluxG_innerfp(l,nsp,1)
!        !    !elseif(LC_sideth==4)then
!        !    !    write(*,*)'ff',index,LC_sideth,innerDisOri(1),cellset(index).fluxF_innerfp(l,2,1)
!        !    !end if
!        !    
!        !end if
!    !    if(index==26527 .and. nearCellIndex==26528.and.l==1)then
!    !        write(*,*)'26527',l,u(:,1)
!    !    end if
!        !do m =1,4
!        !    if(isnan(innerDisOri(m)))then
!        !        write(*,*)'innerDisOri(m)',nt_temp,index,m,innerDisOri
!        !        stop
!        !    end if
!        !end do
!    end if
!end subroutine get_oriL
!
!subroutine get_oriR(index,nearCellIndex,LC_sideth,RC_sideth,l,oriR)
!    !�ǽṹ�ı������в����򣬿��Ǽ������Ӧ�Ĳ�߱�Ź�ϵ������ͨ����������
!    !����һ�ַ�������һ���м�������ڵ�Ԫ���С���ת��ʹ֮��Ӧ������1-3��2-4��Ӧ������
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none 
!    integer :: index,nearCellIndex,LC_sideth,l,m,j,k
!    integer :: RC_sideth,RC_l
!    real(prec),dimension(:,:) :: u(3,4),ll(4,4),rr(4,4),chara_u(3,4),chara_con(3,4)
!    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4),varConL(nsp,nsp,4),varConR(nsp,nsp,4)
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),conR(4),oriR(4),chara_conR(4)
!    real(prec) :: realDis
!    
!    !ȡ�����ҵ�Ԫ�ϵ�ԭʼ����      
!    varOriL(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!    varOriR(:,:,:) = cellset(nearCellIndex).spvalue_ori(:,:,:)   
!    !ȡ�����ҵ�Ԫ�ϵ��غ���� 
!    varConL(:,:,:) = cellset(index).spvalue_con(:,:,:)
!    varConR(:,:,:) = cellset(nearCellIndex).spvalue_con(:,:,:)
!    if(cellset(nearCellIndex).Beta == 0)then
!        do m = 1,4
!            if(LC_sideth==1)then
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(:,RC_l,m),nsp)
!                elseif(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(RC_l,:,m),nsp)
!                else if(RC_sideth == 3)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(:,RC_l,m),nsp)
!                else if(RC_sideth == 4)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(RC_l,:,m),nsp)
!                end if
!            elseif(LC_sideth==2)then
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(:,RC_l,m),nsp)
!                elseif(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(RC_l,:,m),nsp)
!                else if(RC_sideth == 3)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(:,RC_l,m),nsp)
!                else if(RC_sideth == 4)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(RC_l,:,m),nsp)
!                end if        
!            elseif(LC_sideth==3)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(:,RC_l,m),nsp)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(RC_l,:,m),nsp)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(:,RC_l,m),nsp)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(RC_l,:,m),nsp)
!                end if          
!            elseif(LC_sideth==4)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(:,RC_l,m),nsp)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(RC_l,:,m),nsp)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_r,SPs,varOriR(:,RC_l,m),nsp)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    oriR(m) = LaI_nPs(kesi_l,SPs,varOriR(RC_l,:,m),nsp)
!                end if         
!            end if
!        end do
!    else if(cellset(nearCellIndex).Beta == 1)then
!        select case(var_type)
!        case(ori_type)
!            !----ԭʼ����-------------------------------------
!            if(LC_sideth==1)then
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:)   
!                else if(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)               
!                else if(RC_sideth == 3)then
!                    RC_l = l    
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = l  
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if
!            elseif(LC_sideth==2)then     
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:) 
!                else if(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)     
!                else if(RC_sideth == 3)then
!                    RC_l = l
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if        
!            elseif(LC_sideth==3)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if          
!            elseif(LC_sideth==4)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if         
!            end if
!            !write(*,*)'RR'
!            call get_nearCellInfo(nearCellIndex,index,RC_sideth,LC_sideth,RC_l,varOriL,u(:,:),realDis)
!            call FaceFluxC2NNW2(nearCellIndex,index,RC_sideth,RC_l,u,realDis,oriR(:),innerDisOri)
!            !----ԭʼ����-------------------------------------
!        case(con_type)
!            !!----�غ����-------------------------------------
!             if(LC_sideth==1)then
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(1,RC_l,:)
!                    u(3,:) = varConR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(RC_l,nsp,:)
!                    u(3,:) = varConR(RC_l,nsp-1,:)               
!                else if(RC_sideth == 3)then
!                    RC_l = l    
!                    u(2,:) = varConR(nsp,RC_l,:)
!                    u(3,:) = varConR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = l  
!                    u(2,:) = varConR(RC_l,1,:)
!                    u(3,:) = varConR(RC_l,2,:)
!                end if
!            elseif(LC_sideth==2)then     
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(RC_l,1,:)
!                    u(3,:) = varConR(RC_l,2,:)      
!                else if(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(RC_l,nsp,:)
!                    u(3,:) = varConR(RC_l,nsp-1,:) 
!                else if(RC_sideth == 3)then
!                    RC_l = l
!                    u(2,:) = varConR(nsp,RC_l,:)
!                    u(3,:) = varConR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = l
!                    u(2,:) = varConR(RC_l,1,:)
!                    u(3,:) = varConR(RC_l,2,:)
!                end if        
!            elseif(LC_sideth==3)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    u(2,:) = varConR(1,RC_l,:)
!                    u(3,:) = varConR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    u(2,:) = varConR(RC_l,nsp,:)
!                    u(3,:) = varConR(RC_l,nsp-1,:)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(nsp,RC_l,:)
!                    u(3,:) = varConR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(RC_l,1,:)
!                    u(3,:) = varConR(RC_l,2,:)
!                end if          
!            elseif(LC_sideth==4)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    u(2,:) = varConR(1,RC_l,:)
!                    u(3,:) = varConR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    u(2,:) = varConR(RC_l,nsp,:)
!                    u(3,:) = varConR(RC_l,nsp-1,:)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(nsp,RC_l,:)
!                    u(3,:) = varConR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varConR(RC_l,1,:)
!                    u(3,:) = varConR(RC_l,2,:)
!                end if         
!            end if        
!            call get_nearCellInfo(nearCellIndex,index,RC_sideth,LC_sideth,RC_l,varConL,u(:,:),realDis)
!            call FaceFluxC2NNW2(nearCellIndex,index,RC_sideth,RC_l,u,realDis,conR(:),innerDisCon)
!
!            call Func_con_to_ori(conR(:),oriR(:))
!            call Func_con_to_ori(innerDisCon(:),innerDisOri(:))
!            !!----�غ����-------------------------------------
!        case(character_type)
!            !��������
!             if(LC_sideth == 1)then
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:)    
!                else if(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)               
!                else if(RC_sideth == 3)then
!                    RC_l = l    
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = l  
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if
!            elseif(LC_sideth==2)then     
!                if(RC_sideth == 1)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:) 
!                else if(RC_sideth == 2)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)  
!                else if(RC_sideth == 3)then
!                    RC_l = l
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if        
!            elseif(LC_sideth==3)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)
!                else if(RC_sideth == 3)then
!                    RC_l = l
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if          
!            elseif(LC_sideth==4)then
!                if(RC_sideth == 1)then
!                    RC_l = l
!                    u(2,:) = varOriR(1,RC_l,:)
!                    u(3,:) = varOriR(2,RC_l,:)
!                else if(RC_sideth == 2)then
!                    RC_l = l
!                    u(2,:) = varOriR(RC_l,nsp,:)
!                    u(3,:) = varOriR(RC_l,nsp-1,:)
!                else if(RC_sideth == 3)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(nsp,RC_l,:)
!                    u(3,:) = varOriR(nsp-1,RC_l,:)
!                else if(RC_sideth == 4)then
!                    RC_l = nsp+1-l
!                    u(2,:) = varOriR(RC_l,1,:)
!                    u(3,:) = varOriR(RC_l,2,:)
!                end if         
!            end if        
!            call get_nearCellInfo(nearCellIndex,index,RC_sideth,LC_sideth,RC_l,varOriL,u(:,:),realDis)
!            
!            call proj_matrix(u(2,:),ll,rr,RC_sideth)  
!            do k = 1,3
!                call Characteristic_projection(u(k,:),ll,chara_con(k,:)) !��k���� 
!            end do            
!            call FaceFluxC2NNW2(nearCellIndex,index,RC_sideth,RC_l,chara_con,realDis,chara_conR(:),chara_innerDisCon)
!            call Inverse_Characteristic_projection(chara_conR(:),rr,oriR(:))
!            call Inverse_Characteristic_projection(chara_innerDisCon(:),rr,innerDisOri(:))
!            
!        end select        
!        
!        if(method_subcell == method_NNW)then
!            if(RC_sideth==1)then    !�ӵ�Ԫ�ȰѼ��ͨ����������֮��絥Ԫ��������һ��
!                cellset(nearCellIndex).fluxG_innerfp(RC_l,2,:) = innerDisOri
!            elseif(RC_sideth==2)then
!                cellset(nearCellIndex).fluxF_innerfp(RC_l,nsp,:) = innerDisOri
!            elseif(RC_sideth==3)then
!                cellset(nearCellIndex).fluxG_innerfp(RC_l,nsp,:) = innerDisOri
!            elseif(RC_sideth==4)then
!                cellset(nearCellIndex).fluxF_innerfp(RC_l,2,:) = innerDisOri
!            end if
!        elseif(method_subcell == method_Godunov)then
!            oriR(:) = u(2,:)
!        end if
!        !write(*,*)'gg',nearCellIndex,RC_sideth,innerDisOri(:)
!        !do m =1,4
!            !if(innerDisOri(1)<0.0_prec)then
!            !    write(*,*)'22',nearCellIndex,oriR(l,1),innerDisOri(1)
!            !end if
!        !end do
!        !if(index==26527 .and. nearCellIndex==26528.and.l==5)then
!        !    write(*,*)'26528',u(:,1)
!        !end if
!        !if(index==26527 .and. nearCellIndex==26528.and.l==1)then
!        !    write(*,*)'26528',l,u(:,1)
!        !end if
!        !if(index==9)then
!        !    write(*,*)nt_temp,l,u(:,1)
!        !end if
!        !if(nearCellIndex==160)then
!        !    write(*,*)'222222222',l
!        !    write(*,*)chara_con
!        !    !if(RC_sideth==1)then
!        !    !    write(*,*)'gg',nearCellIndex,RC_sideth,innerDisOri(1),cellset(nearCellIndex).fluxG_innerfp(l,2,1)
!        !    !elseif(RC_sideth==2)then
!        !    !    write(*,*)'ff',nearCellIndex,RC_sideth,innerDisOri(1),cellset(nearCellIndex).fluxF_innerfp(l,nsp,1)
!        !    !elseif(RC_sideth==3)then
!        !    !    write(*,*)'gg',nearCellIndex,RC_sideth,innerDisOri(1),cellset(nearCellIndex).fluxG_innerfp(l,nsp,1)
!        !    !elseif(RC_sideth==4)then
!        !    !    write(*,*)'ff',nearCellIndex,RC_sideth,innerDisOri(1),cellset(nearCellIndex).fluxF_innerfp(l,2,1)
!        !    !end if
!        !end if
!    end if
!end subroutine get_oriR
 !if(method_subcell == method_NNW)then
        !    if(LC_sideth==1)then    !�ӵ�Ԫ�Ȱѽ���ԭʼ������������֮��絥Ԫ��������һ��
        !        oriL(:) = oriL(:)/cellset(index).fpdet_J_G(l,1)
        !        innerDisOri(:)   = innerDisOri(:)/cellset(index).fpdet_J_G(l,1)
        !        cellset(index).fluxG_innerfp(l,2,:) = innerDisOri(:)         
        !    elseif(LC_sideth==2)then
        !        oriL(:) = oriL(:)/cellset(index).fpdet_J_F(l,nsp+1)
        !        innerDisOri(:)   = innerDisOri(:)/cellset(index).fpdet_J_F(l,nsp+1)
        !        cellset(index).fluxF_innerfp(l,nsp,:) = innerDisOri(:)
        !    elseif(LC_sideth==3)then
        !        oriL(:) = oriL(:)/cellset(index).fpdet_J_G(l,nsp+1)
        !        innerDisOri(:)   = innerDisOri(:)/cellset(index).fpdet_J_G(l,nsp+1)
        !        cellset(index).fluxG_innerfp(l,nsp,:) = innerDisOri(:)
        !    elseif(LC_sideth==4)then
        !        oriL(:) = oriL(:)/cellset(index).fpdet_J_F(l,1)
        !        innerDisOri(:)   = innerDisOri(:)/cellset(index).fpdet_J_F(l,1)
        !        cellset(index).fluxF_innerfp(l,2,:) = innerDisOri(:)
        !    end if
        !elseif(method_subcell == method_Godunov)then
        !    if(LC_sideth==1)then    
        !        oriL(:) = u(2,:)/cellset(index).det_J(1,l)       
        !    elseif(LC_sideth==2)then
        !        oriL(:) = u(2,:)/cellset(index).det_J(l,nsp)
        !    elseif(LC_sideth==3)then
        !        oriL(:) = u(2,:)/cellset(index).det_J(nsp,l)
        !    elseif(LC_sideth==4)then
        !        oriL(:) = u(2,:)/cellset(index).det_J(l,1)
        !    end if
        !end if
!call get_nearCellInfo(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,varOriR,u(:,:),realDis1)  !u1��������ԭʼ������u2,u3�Ǽ�����          
!            call proj_matrix(u(2,:),ll,rr,LC_sideth)  !(uԭʼ����)
!            do k = 1,3
!                call Characteristic_projection(u(k,:),ll,chara_con(k,:)) !(chara_con�����غ����)
!                !write(*,"(4F10.5)")chara_u(k,:)
!            end do
!            call FaceFluxC2NNW2(index,nearCellIndex,LC_sideth,l,chara_con,realDis1,chara_conL,chara_innerDisCon)
!            call Inverse_Characteristic_projection(chara_conL(:),rr,oriL(:))
!            call Inverse_Characteristic_projection(chara_innerDisCon(:),rr,innerDisOri(:))        
    !do k = 1,nsp
    !    do l =1,nsp           
    !        kesix = cellset(index).Mdirect(k,l,1)
    !        kesiy = cellset(index).Mdirect(k,l,2)
    !        etax  = cellset(index).Mdirect(k,l,3)
    !        etay  = cellset(index).Mdirect(k,l,4)
    !        Q1(k,l,:) = detJ(k,l)*cellset(index).spvalue_con(k,l,:)
    !        F1(k,l,:) = detJ(k,l)*(kesix*cellset(index).spvalue_fluF(k,l,:) + kesiy*cellset(index).spvalue_fluG(k,l,:))
    !        G1(k,l,:) = detJ(k,l)*(etax *cellset(index).spvalue_fluF(k,l,:) + etay *cellset(index).spvalue_fluG(k,l,:))  
    !    end do
    !end do 
   !if(method_subcell == method_NNW)then
   !         if(RC_sideth==1)then    !�ӵ�Ԫ�ȰѼ��ͨ����������֮��絥Ԫ��������һ��
   !             oriR(:) = oriR(:)
   !             innerDisOri(:)  = innerDisOri(:)
   !             cellset(nearCellIndex).fluxG_innerfp(RC_l,2,:) = innerDisOri
   !         elseif(RC_sideth==2)then
   !             oriR(:) = oriR(:)
   !             innerDisOri(:)  = innerDisOri(:)
   !             cellset(nearCellIndex).fluxF_innerfp(RC_l,nsp,:) = innerDisOri
   !         elseif(RC_sideth==3)then
   !             oriR(:) = oriR(:)
   !             innerDisOri(:)  = innerDisOri(:)
   !             cellset(nearCellIndex).fluxG_innerfp(RC_l,nsp,:) = innerDisOri
   !         elseif(RC_sideth==4)then
   !             oriR(:) = oriR(:)
   !             innerDisOri(:)  = innerDisOri(:)
   !             cellset(nearCellIndex).fluxF_innerfp(RC_l,2,:) = innerDisOri
   !         end if
   !     elseif(method_subcell == method_Godunov)then
   !         !oriR = u2 Godunov 1st oeder
   !         if(RC_sideth==1)then    
   !             call Func_ComTranPhy(nearCellIndex,cellset(nearCellIndex).det_J(1,l),u(2,:),oriR)
   !         elseif(RC_sideth==2)then
   !             call Func_ComTranPhy(nearCellIndex,cellset(nearCellIndex).det_J(l,nsp),u(2,:),oriR)
   !         elseif(RC_sideth==3)then
   !             call Func_ComTranPhy(nearCellIndex,cellset(nearCellIndex).det_J(nsp,l),u(2,:),oriR)
   !         elseif(RC_sideth==4)then
   !             call Func_ComTranPhy(nearCellIndex,cellset(nearCellIndex).det_J(l,1),u(2,:),oriR)
   !         end if
   !     end if
                !WCNS ��߽���C2NNW2��ֵ�õ�����ֵ����Ԫ�ڲ����ӵ�Ԫʹ��WCNS��ֵ
                !do m = 1,4
                !    x = FPs(l)
                !    xk = SPs(l-1:l+1)
                !    oriDis(l,2,m) = NonLinear_WCNS(x,xk,u(:,m))
                !    x = FPs(l+1)
                !    oriDis(l+1,1,m) = NonLinear_WCNS(x,xk,u(:,m))
                !end do
    !do k = 1,nsp
    !    do l =1,nsp           
    !        kesix = cellset(index).Mdirect(k,l,1)
    !        kesiy = cellset(index).Mdirect(k,l,2)
    !        etax  = cellset(index).Mdirect(k,l,3)
    !        etay  = cellset(index).Mdirect(k,l,4)
    !        Q1(k,l,:) = detJJ(k,l)*cellset(index).spvalue_con(k,l,:)
    !        F1(k,l,:) = detJJ(k,l)*(kesix*cellset(index).spvalue_fluF(k,l,:) + kesiy*cellset(index).spvalue_fluG(k,l,:))
    !        G1(k,l,:) = detJJ(k,l)*(etax *cellset(index).spvalue_fluF(k,l,:) + etay *cellset(index).spvalue_fluG(k,l,:))  
    !    end do
    !end do 


!subroutine get_nearCellInfo(index,nearCellIndex,Cell_sideth,nearCell_sideth,l,near_l,CellOri,nearCellOri,ui,uA,realDis)
!    !C2NNW2�絥Ԫ����ȡ�ڵ�Ԫ��ģ��㣬�����ڵ�Ԫ�Ķ�Ӧ��ϵȡ��Ӧ��
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none 
!    integer :: index,nearCellIndex,Cell_sideth,nearCell_sideth,l,near_l,m,i
!    real(prec),dimension(:)   :: coor_C(2),coor_P(2),coor_Pfp(2),uA(4)!��׼��Ԫ��ʵ�ʵ�Ԫ
!    real(prec),dimension(:,:) :: cell_vertex(4,2)   !ʵ�ʵ�Ԫ���㣬��׼��Ԫ��֪ 
!    real(prec),dimension(:,:) :: oriR(nsp,4),ui(3,4)
!    real(prec),dimension(:,:,:) :: nearCellOri(nsp,nsp,4),CellOri(nsp,nsp,4)
!    real(prec) :: realDis,kesifp,etafp,det_u1,det_u2
!    integer :: ikesi,ieta
!    !write(*,*)'l',l
!    if(nearCell_sideth == 1)then
!        ikesi = 1;            ieta  = near_l
!        kesifp = SPs_local(ikesi,ieta,1); etafp = -1.0_prec
!    elseif(nearCell_sideth == 2)then
!        ikesi = near_l;            ieta  = nsp
!        kesifp = 1.0_prec; etafp = SPs_local(ikesi,ieta,2)
!    else if(nearCell_sideth == 3)then
!        ikesi = nsp;            ieta  = near_l
!        kesifp = SPs_local(ikesi,ieta,1); etafp = 1.0_prec
!    else if(nearCell_sideth == 4)then
!        ikesi = near_l;            ieta  = 1
!        kesifp = -1.0_prec; etafp = SPs_local(ikesi,ieta,2) 
!    end if
!    det_u1 = cellset(nearCellIndex).det_J(ikesi,ieta)!u1����Jacobi
!    !ȡ�����ڵ�Ԫ��ģ���
!    ui(1,:) = nearCellOri(ikesi,ieta,:)    
!    !write(*,"('0',F20.9)") ui(1,1)
!    call Func_ComTranPhy(nearCellIndex,det_u1,ui(1,:),ui(1,:))  !�Ӽ���ռ�任������ռ䣬����������Ҫ�����򷴾����Ȩ�õ�ͨ����ֵ
!    !write(*,*) nearCellIndex,det_u1,cellset(nearCellIndex).det_J
!    !write(*,"('1',4F20.9)") ui(1,1)
!    !write(*,"('2',4F20.9)") cellset(nearCellIndex).spvalue_ori(ikesi,ieta,:)
!    !write(*,"(4F20.9)") 
!    if(case_comp == Riemann2D_case .and. nearCellIndex > ncells)then
!        if(Cell_sideth==1)then
!            det_u2 = cellset(index).det_J(1,l)
!        elseif(Cell_sideth==2)then
!            det_u2 = cellset(index).det_J(l,nsp)
!        elseif(Cell_sideth==3)then
!            det_u2 = cellset(index).det_J(nsp,l)
!        elseif(Cell_sideth==4)then
!            det_u2 = cellset(index).det_J(l,1)
!        end if
!        call Func_ComTranPhy(nearCellIndex,det_u2,ui(2,:),ui(1,:))   !����ռ���u1 = u2
!    end if
!    !�����1����������
!    do i = 1,4
!        cell_vertex(i,:) = xy_coor(cellset(nearCellIndex).nodes(i),:)
!    end do
!    coor_C = SPs_local(ikesi,ieta,:)
!    call quad_C2Phy(cell_vertex,coor_C,coor_P)
!    coor_C(1) = kesifp;    coor_C(2) = etafp
!    call quad_C2Phy(cell_vertex,coor_C,coor_Pfp)
!    !write(*,*)Cell_sideth,l,ikesi,ieta
!    realDis = sqrt((coor_P(1)-coor_Pfp(1))**2+(coor_P(2)-coor_Pfp(2))**2)
!end subroutine get_nearCellInfo
!
!subroutine FaceFluxC2NNW2(index,nearCellIndex,Cell_sideth,l,u,uA,realDis,u_cell_L,u_cell_R)
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    real(prec),dimension(:,:) :: u(3,4)
!    real(prec),dimension(2) :: coor_C,coor_P1,coor_P2,coor_Pfp
!    real(prec),dimension(4) :: uA1,uB1,uA2,uB2,du1,du2,du,u_cell_L(4),u_cell_R(4),u2_phy(4),uA(4)
!    real(prec) :: d1,d2,d3,d4,dd1,dd2,dd3,dd4
!    real(prec) :: cc,w1,w2,w3,w4,w5,w6,Fai,varMax,varMin,limA,limB    
!    real(prec),external :: lim
!    real(prec),dimension(:,:) :: cell_vertex(4,2)   !ʵ�ʵ�Ԫ����
!    integer :: index,nearCellIndex,Cell_sideth,m,l,j,ikesi2,ieta2,k
!    real(prec) :: kesifp,etafp,realDis,det_A,detJ_u2
!    real(prec),dimension(:,:) :: ll(4,4),rr(4,4)
!    !ͨ��������֮��ľ���
!    !d1,d2�漰����Ԫ���磬��ʵ�ʾ��� real distance
!    !d3,d4�ڵ�Ԫ�ڲ�����local distance
!    
!    !ȡ����Ԫ����
!    do j = 1,4
!        cell_vertex(j,:) = xy_coor(cellset(index).nodes(j),:)
!    end do
!    !(ikesi2��ieta2) �ǵ�2�ڵ�Ԫ�ı�ţ�(kesifp,etafp)ͨ���������
!    if(Cell_sideth==1)then      !��1�ϵ�ͨ����
!        ikesi2  = 1;     ieta2 = l;
!        kesifp = SPs_local(ikesi2,ieta2,1);    etafp = -1.0_prec;
!        det_A = cellset(index).fpdet_J_G(l,1)
!    elseif(Cell_sideth==2)then  !��2�ϵ�ͨ����
!        ikesi2 = l;     ieta2 = nsp;
!        kesifp = 1.0_prec;     etafp = SPs_local(ikesi2,ieta2,2);
!        det_A = cellset(index).fpdet_J_F(l,nsp+1)
!    elseif(Cell_sideth==3)then  !��3�ϵ�ͨ����
!        ikesi2 = nsp;   ieta2 = l;
!        kesifp = SPs_local(ikesi2,ieta2,1);    etafp = 1.0_prec;
!        det_A = cellset(index).fpdet_J_G(l,nsp+1)
!    elseif(Cell_sideth==4)then  !��4�ϵ�ͨ����
!        ikesi2 = l;     ieta2 = 1;
!        kesifp = -1.0_prec;     etafp = SPs_local(ikesi2,ieta2,2);
!        det_A = cellset(index).fpdet_J_F(l,1)
!    end if
!
!    !��2
!    coor_C = SPs_local(ikesi2,ieta2,:)
!    call quad_C2Phy(cell_vertex,coor_C,coor_P2)
!    !��fp_1/fp_nsp+1
!    coor_C(1) = kesifp; coor_C(2) = etafp
!    call quad_C2Phy(cell_vertex,coor_C,coor_Pfp)
!    
!    d1 = realDis
!    d2 = sqrt((coor_Pfp(1)-coor_P2(1))**2 + (coor_Pfp(2)-coor_P2(2))**2 ) 
!    !write(*,*) 'd',Cell_sideth,d1/d2
!    d3 = dis_sp_fp(2)
!    d4 = dis_sp_fp(3)
!    
!    !������Ȩ
!    cc = 1.0_prec!����1������д���������滻
!    dd1 = cc/d1
!    dd2 = cc/d2    
!    dd3 = cc/d3
!    dd4 = cc/d4
!    
!    w1 = dd1/(dd1+dd2)
!    w2 = dd2/(dd1+dd2)
!    w3 = dd3/(dd3+dd4)
!    w4 = dd4/(dd3+dd4)
!    
!    !ͨ���㴦ֵ
!    detJ_u2 = cellset(index).det_J(ikesi2,ieta2)
!    call Func_ComTranPhy(index,detJ_u2,u(2,:),u2_phy)
!    !write(*,"('1',2F20.9)") u(1,1),u2_phy(1)
!    uA1 = w1*u(1,:)+w2*u2_phy!����ռ䣬u(1,:)�洢��������ռ�ı���ֵ
!    !write(*,"(F20.9)")uA1(:)
!    write(*,*)uA1(:)
!    call Func_PhyTranCom(index,det_A,uA1,uA1)
!
!    !write(*,"('2',4F20.9)") cellset(index).spvalue_ori(ikesi2,ieta2,:)
!    if(var_type == character_type)then  !���ѡ������������ֵ������Ҫ���������任
!        call proj_matrix(u(2,:),ll,rr,Cell_sideth)  
!        do k = 2,3
!            call Characteristic_projection(u(k,:),ll,u(k,:)) !��k���� ����u������,�������������
!        end do   
!        call Characteristic_projection(uA1,ll,uA1)
!    end if
!    uB1 = w3*u(2,:)+w4*u(3,:)   !������ 
!    !��Ԫ����
!    d2 = dis_sp_fp(1)   !���������ڵ�Ԫ�ڲ�������ʹ��local distance
!    dd2 = cc/d2
!    w5 = dd2/(dd2+dd3)
!    w6 = dd3/(dd2+dd3)
!    du1 = (u(2,:)-uA1)/d2
!    du2 = (uB1-u(2,:))/d3
!    du = w5*du1 + w6*du2
!    !write(*,"('3',F20.9)")uA2(1)
!    !���¼���uA,uB
!    uA2 = u(2,:)-du*d2
!    uB2 = u(2,:)+du*d3
!    !u(1,:) = uA1            !�絥Ԫ��Ϊ��varMax,varMinĿǰ������������
!    !u(1,:)�任�����������������ò�����������Ҫ
!    call Func_PhyTranCom(index,det_A,u(1,:),u(1,:))     
!    if(var_type == character_type)then
!        call Characteristic_projection(u(1,:),ll,u(1,:))
!    end if
!    !!---------------------------------------------------------------------------
!    do m = 1,4
!        varMax = maxVal(u(:,m))
!        varMin = minVal(u(:,m))
!        limA = lim(u(2,m),uA2(m),varMax,varMin)
!        limB = lim(u(2,m),uB2(m),varMax,varMin)
!        Fai = min(limA,limB)    
!        !��Ԫ��ֵ�ǽ�����ֵ����Ԫ��ֵ�ǽ�����ֵ
!        u_cell_L(m) = u(2,m)-Fai*du(m)*d2
!        u_cell_R(m) = u(2,m)+Fai*du(m)*d3        
!    end do
!    if(var_type == character_type)then  !���ѡ������������ֵ������Ҫ���������任
!        call Inverse_Characteristic_projection(u_cell_L(:),rr,u_cell_L(:))
!        call Inverse_Characteristic_projection(u_cell_R(:),rr,u_cell_R(:))
!    end if
!end subroutine FaceFluxC2NNW2
!subroutine indicator_TVB2(index,Beta)
!    !TVB��⣬�����غ����
!    ! ÿһά�����бȽϣ����ֲ�����minmod�����ı���Ϊ���ⵥԪ��ֹͣ����
!    !Refs:
!    !   
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m,sideth,near_i 
!    integer :: index,Beta,indexNearCell
!    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh2
!    real(prec),dimension(4) :: aveNearCell
!    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
!    real(prec),dimension(nsp,nsp,4) :: varOriCell
!    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
!    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa
!   
!    !!!����Ԫ -- ������ԭʼ����
!    varOriCell(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!    !�ڵ�Ԫ
!    do i = 1,4
!        indexNearCell  = cellset(index).nearcells(i)
!        varOriCellNearCell(i,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
!    end do 
!    !��ά���
!    !x����
!    loop1:do m = 1,4,3!ÿ����������Ҫ���
!        !write(*,*)m
!        loop2:do i = 1,nsp!ÿ�зֱ����
!            !dh2 = 1+SPs(i)**2   !4M*dh^2
!            dh2 = (2.0_prec*cellset(index).MJacobi(i,1,1))**2
!            !write(*,*)TVB_M*dh2
!            !����Ԫ�����άƽ��ֵ��ԭʼ����//�����㲻ΪGauss�㣬����Ҫ�Ȳ�ֵ��Gauss�����ٻ��֣������Ժ����ӵ�Ԫ���������
!            aveCell = Gauss_integral_SPs(varOriCell(i,:,m))/2.0_prec
!
!            varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(i,:,m),nsp)
!            varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(i,:,m),nsp)
!            
!            !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!            indexNearCell  = cellset(index).nearcells(2)
!            do j = 1,4
!                if(cellset(indexNearCell).nearcells(j)==index)then
!                    sideth = j
!                end if
!            end do
!            if((sideth==1).or.(sideth==2))then
!                near_i = nsp+1-i
!            else
!                near_i = i
!            end if
!            if(sideth==2.OR.sideth==4)then
!                aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,near_i,:,m))/2.0_prec
!            else
!                aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,:,near_i,m))/2.0_prec
!            end if
!            
!            !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!            indexNearCell  = cellset(index).nearcells(4)
!            !write(*,*) cellset(index).fpdet_J_F(i,1),cellset(indexNearCell).fpdet_J_F(i,nsp+1)
!            do j = 1,4
!                if(cellset(indexNearCell).nearcells(j)==index)then
!                    sideth = j
!                end if
!            end do
!            if((sideth==3).or.(sideth==4))then
!                near_i = nsp+1-i
!            else
!                near_i = i
!            end if
!            if(sideth==2.OR.sideth==4)then
!                aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,near_i,:,m))/2.0_prec
!            else
!                aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,:,near_i,m))/2.0_prec
!            end if
!  
!            deta_L =aveCell-aveNearCell(4)
!            deta_R = aveNearCell(2)-aveCell           
!            deta_CL = aveCell-varFluxL
!            deta_CR = varFluxR-aveCell
!
!            varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh2)      
!            varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh2)
!            !if(m==1)then
!            !    !write(*,*)'1', index,varmodR
!            !    !write(*,*)'2',TVB_M*dh2,
!            !    write(*,*)varmodL,deta_CL,aveCell,varFluxL
!            !end if
!            !write(*,*)aveCell-varFluxL,aveNearCell(2)-aveCell,aveCell-aveNearCell(4)
!            if(varmodL /= deta_CL .or. varmodR /= deta_CR)then
!                Beta = 1
!                !cellset(index).Beta_line(i,1)=1
!                if(detect_type == ByDim )then
!                    cellset(index).Beta_line(i,1)=1
!                elseif(detect_type == ByCell)then
!                    cellset(index).Beta_line = 1
!                    write(*,*)TVB_M*dh2
!                    exit loop1
!                end if
!                
!            end if           
!            
!            !y����
!            dh2 = (2.0_prec*cellset(index).MJacobi(1,i,4))**2
!            aveCell = Gauss_integral_SPs(varOriCell(:,i,m))/2.0_prec
!
!            varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(:,i,m),nsp)
!            varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(:,i,m),nsp)        
!            
!            !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!            indexNearCell = cellset(index).nearcells(1)
!            do j = 1,4
!                if(cellset(indexNearCell).nearcells(j)==index)then
!                    sideth = j
!                end if
!            end do
!            if((sideth==2).or.(sideth==1))then
!                near_i = nsp+1-i
!            else
!                near_i = i
!            end if
!            if(sideth==2.OR.sideth==4)then
!                aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,near_i,:,m))/2.0_prec
!            else
!                aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,:,near_i,m))/2.0_prec
!            end if
!            !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!            indexNearCell = cellset(index).nearcells(3)
!            do j = 1,4
!                if(cellset(indexNearCell).nearcells(j)==index)then
!                    sideth = j
!                end if
!            end do
!            if((sideth==3).or.(sideth==4))then
!                near_i = nsp+1-i
!            else
!                near_i = i
!            end if
!            if(sideth==2.OR.sideth==4)then
!                aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,near_i,:,m))/2.0_prec
!            else
!                aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,:,near_i,m))/2.0_prec
!            end if
!            !varmodL = TVB_minmod(aveCell-varFluxL,aveNearCell(3)-aveCell,aveCell-aveNearCell(1),dh2)
!            !varmodR = TVB_minmod(varFluxR-aveCell,aveNearCell(3)-aveCell,aveCell-aveNearCell(1),dh2)
! 
!            deta_L = aveCell-aveNearCell(1)
!            deta_R = aveNearCell(3)-aveCell
!            deta_CL = aveCell-varFluxL
!            deta_CR = varFluxR-aveCell
!            
!            varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh2)           
!            varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh2)
!            
!            !if(index==165)then
!            !    write(*,*)'y---------------------'
!            !    write(*,*) TVB_M*dh2
!            !    write(*,*)varmodL,aveCell-varFluxL
!            !    write(*,*)varmodR,varFluxR-aveCell                
!            !    write(*,*)'---------------------'
!            !end if
!            if(varmodL /= deta_CL .or. varmodR /= deta_CR)then
!                Beta = 1
!                !cellset(index).Beta_line(i,2)=1
!                if(detect_type == ByDim )then
!                    cellset(index).Beta_line(i,2)=1
!                elseif(detect_type == ByCell)then
!                    cellset(index).Beta_line = 1
!                    exit loop1
!                end if
!            end if        
!        end do loop2
!    end do loop1
!   
!end subroutine indicator_TVB2

  !!������������߽絥Ԫָ���Ӧλ�ã�֮�󣬲���Ҫ������ֵ
  !      do i = 1,nbdsides
  !          do j = 1,4
  !              indexCell = BoundCells_index_set(i)
  !              indexNearCell = cellset(indexCell).nearcells(j)
  !              if(indexNearCell>ncells)then
  !                  !index_temp = cellset(indexNearCell).index       !����߽絥ԪӦ�õ�λ��ָ��
  !                  !cellset(indexNearCell) = cellset(index_temp)    !��ֵ
  !                  !cellset(indexNearCell).index = index_temp       !����λ��ָ��
  !          
  !                  !���߲�����������Ϊ����index����ͬ
  !                  !write(*,*)indexNearCell,cellset(indexNearCell).spvalue_ori(1,1,1)
  !                  !cellset(indexNearCell).spvalue_ori = cellset(cellset(indexNearCell).index).spvalue_ori
  !                  !write(*,*)indexNearCell,cellset(indexNearCell).spvalue_ori(1,1,1)!ǰ��Ա���֤
  !                  !stop
  !                  !cellset(indexNearCell) = cellset(cellset(indexNearCell).index)    !��ֵ
  !                  cellset(indexNearCell).det_J = cellset(index).det_J
  !                  cellset(indexNearCell).MJacobi = cellset(index).MJacobi
  !                  cellset(indexNearCell).Mdirect = cellset(index).Mdirect 
  !                  cellset(indexNearCell).spvalue_ori= cellset(index).spvalue_ori 
  !                  !write(*,*)indexNearCell,cellset(indexNearCell).spvalue_con(1,1,1) !ǰ��Ա���֤һ��û�и�����ı߽絥Ԫ�����ڴ棬��������𣿻����Զ��������ڴ�
  !              end if            
  !          end do        
  !      end do
  !  end if
!subroutine print_TrouCell2
!    !��ӡ���ⵥԪ
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer :: i,j,k,l,sumTC 
!    character(len=100)filename,filename2
!    character(len=5)char_nsdpx,char_nsdpy,char_nsdx,char_nsdy
!      
!    !----���.plt�鿴 Mesh------------------------------------------------------------------------------   
!    select case(grid_set)
!    case(self)
!        write(char_nsdpx,"(TL1,I4)") nsdpx-1
!        write(char_nsdpy,"(TL1,I4)") nsdpy-1
!        filename = '.\Result\TrouCell_x_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'.plt'   
!        filename2 = '.\Result\TrouCell_y_nx_'//trim(adjustl(char_nsdpx))//'_ny_'//trim(adjustl(char_nsdpy))//'.plt' 
!    case(fluent_cas)
!        write(char_nsdx,"(TL1,I4)") nsdx
!        write(char_nsdy,"(TL1,I4)") nsdy
!        filename = '.\Result\TrouCell_x_'//trim(mesh_file)//'.plt'
!        filename2 = '.\Result\TrouCell_y_'//trim(mesh_file)//'.plt'
!    end select
!    
!    open(81,status='REPLACE',file = filename)
!    write(81,*)'title="Trouble Cells"'
!    write(81,*) 'variables =x,y'
!    write(81,*) 'ZONE, N=',nnodes,', E=',ncells,', F=FEPOINT, ET=QUADRILATERAL'
!    
!    !FEPOINT��ʽ
!    do i = 1,nnodes
!        write(81,*) xy_coor(i,:)
!    end do
! 
!    !�ı��ε�Ԫ������
!    do i = 1,ncells
!        write(81,*) cellset(i).nodes(1:4)
!    end do
!    sumTC = 0
!    do i = 1,ncells
!        do l = 1,nsp
!            sumTC = sumTC + cellset(i).Beta_line(l,1)
!        end do
!    end do
!    if(sumTC>0)then
!        !�ڲ��������
!        write(81,*) 'ZONE, i=',sumTC*nsp ,', j=',1
!        !��ά���
!        do i = 1,ncells
!            do l = 1,nsp
!                if(cellset(i).Beta_line(l,1) == 1)then
!                    do k =1,nsp         
!                        write(81,*)cellset(i).sp_coor(l,k,1),cellset(i).sp_coor(l,k,2)
!                    end do   
!                end if
!            end do    
!        end do
!    end if
!        
!    open(82,status='REPLACE',file = filename2)
!    write(82,*)'title="Trouble Cells"'
!    write(82,*) 'variables =x,y'
!    write(82,*) 'ZONE, N=',nnodes,', E=',ncells,', F=FEPOINT, ET=QUADRILATERAL'
!    
!    !FEPOINT��ʽ
!    do i = 1,nnodes
!        write(82,*) xy_coor(i,:)
!    end do
! 
!    !�ı��ε�Ԫ������
!    do i = 1,ncells
!        write(82,*) cellset(i).nodes(1:4)
!    end do
!    sumTC = 0
!    do i = 1,ncells
!        do l = 1,nsp
!            sumTC = sumTC + cellset(i).Beta_line(l,2)
!        end do
!    end do
!    if(sumTC>0)then
!        !�ڲ��������
!        write(82,*) 'ZONE, i=',sumTC*nsp ,', j=',1
!        !��ά���
!        do i = 1,ncells
!            do l = 1,nsp
!                if(cellset(i).Beta_line(l,2) == 1)then
!                    do k =1,nsp         
!                        write(82,*)cellset(i).sp_coor(k,l,1),cellset(i).sp_coor(k,l,2)
!                    end do   
!                end if
!            end do    
!        end do
!    end if
!    
!    close(81)
!    close(82)
!end subroutine print_TrouCell2
!!������ⵥԪ
!subroutine mark_TroCell
!    !��Ⲣ������ⵥԪ��
!    use parameter_setting
!    use global_var
!    use type_module
!    implicit none
!    integer :: i,j,k,l
!    do i = 1,ncells+nbdsides
!        !��ʼ��
!        cellset(i).Beta = 0
!        cellset(i).Beta_line = 0
!    end do
!    !ATV��⡣��������TVB֮ǰ�������������ⷽ����������ATV�����TroubleCell����С��Χ��Ȼ������Щ��Ԫ��������TVB��⣬�����ų�TVB�ڼ�ֵ�㴦����ⲻ����
!
!    !call indicator_ATV
!    
!    !buffer_cell_switch = buffer_cell_status            
!    !TVB���
!    !call indicator_TVB
!        
!    !KXRCF���
!    !��д       
!
!    !modal 
!    call indicator_MDH
!    
!    !������ⵥԪ
!    call add_buffer_cell
!end subroutine mark_TroCell
!
!subroutine indicator_TVB
!    !TVB��⣬�����غ����
!    ! ÿһά�����бȽϣ����ֲ�����minmod�����ı���Ϊ���ⵥԪ��ֹͣ����
!    !Refs:
!    !   [1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29. 
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m
!    integer :: sideth1,sideth2,sideth3,sideth4,near_i1,near_i2,near_i3,near_i4,indexNearCell,indexNearCell1, indexNearCell2, indexNearCell3, indexNearCell4 
!    integer :: index,Beta,Beta_line_x,Beta_line_y
!    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh2
!    real(prec),dimension(4) :: aveNearCell
!    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
!    real(prec),dimension(nsp,nsp,4) :: varOriCell
!    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
!    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa
!    do index = 1,ncells        
!        !!����Ԫ
!        !do k = 1,nsp
!        !    do l = 1,nsp
!        !        call Func_con_to_ori(cellset(index).spvalue_con_loc(k,l,:),varOriCell(k,l,:))
!        !    end do 
!        !end do    
!        !!�ڵ�Ԫ -- ������ԭʼ����
!        !do i = 1,4
!        !    indexNearCell  = cellset(index).nearcells(i)
!        !    do k = 1,nsp
!        !        do l = 1,nsp
!        !            call Func_con_to_ori(cellset(indexNearCell).spvalue_con_loc(k,l,:),varOriCellNearCell(i,k,l,:))
!        !        end do 
!        !    end do
!        !end do    
!    
!        !!����Ԫ -- ������ԭʼ����
!        varOriCell(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!
!        !varOriCell(:,:,1) = varOriCell(:,:,1)*varOriCell(:,:,4) !rho*p
!        !varOriCell(:,:,1) = varOriCell(:,:,4)/varOriCell(:,:,1) !p/rho
!        !�ڵ�Ԫ
!        do i = 1,4
!            indexNearCell  = cellset(index).nearcells(i)
!            varOriCellNearCell(i,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,1)*varOriCellNearCell(i,:,:,4) !rho*p
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,4)/varOriCellNearCell(i,:,:,1) !p/rho
!        end do 
!        !��ά���
!        !x����
!        Beta = 0
!
!        loop1:do m = 1,1!������
!            !write(*,*)m
!            loop2:do i = 1,nsp!ÿ�зֱ����
!                Beta_line_x = 0
!                Beta_line_y = 0
!                dh2 = (2.0_prec*cellset(index).MJacobi(i,1,1))**2
!                !write(*,*)2.0_prec*cellset(index).MJacobi(i,1,1)
!                !����Ԫ�����άƽ��ֵ��ԭʼ����//�����㲻ΪGauss�㣬����Ҫ�Ȳ�ֵ��Gauss�����ٻ��֣������Ժ����ӵ�Ԫ���������
!                aveCell = Gauss_integral_SPs(varOriCell(i,:,m))/2.0_prec
!                varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(i,:,m),nsp)
!                varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(i,:,m),nsp)
!                !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell1  = cellset(index).nearcells(4)
!                !write(*,*) cellset(index).fpdet_J_F(i,1),cellset(indexNearCell).fpdet_J_F(i,nsp+1)
!                do j = 1,4
!                    if(cellset(indexNearCell1).nearcells(j)==index)then
!                        sideth1 = j
!                    end if
!                end do
!                if((sideth1==3).or.(sideth1==4))then
!                    near_i1 = nsp+1-i
!                else
!                    near_i1 = i
!                end if
!                if(sideth1==2.OR.sideth1==4)then
!                    aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,near_i1,:,m))/2.0_prec
!                else
!                    aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,:,near_i1,m))/2.0_prec
!                end if
!                !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell2  = cellset(index).nearcells(2)
!                do j = 1,4
!                    if(cellset(indexNearCell2).nearcells(j)==index)then
!                        sideth2 = j
!                    end if
!                end do
!                if((sideth2==1).or.(sideth2==2))then
!                    near_i2 = nsp+1-i
!                else
!                    near_i2 = i
!                end if
!                if(sideth2==2.OR.sideth2==4)then
!                    aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,near_i2,:,m))/2.0_prec
!                else
!                    aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,:,near_i2,m))/2.0_prec
!                end if
!    
!                deta_L =aveCell-aveNearCell(4)
!                deta_R = aveNearCell(2)-aveCell           
!                deta_CL = aveCell-varFluxL
!                deta_CR = varFluxR-aveCell
!
!                varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh2)      
!                varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh2)
!            
!                !write(*,*)TVB_M*dh2,a1
!                !if(m==1)then
!                !    !write(*,*)'1', index,varmodR
!                    !if(index==10)then
!                        !write(*,*)varmodL,deta_CL!,deta_R,deta_L
!                    !end if
!                    !pause
!                    !write(*,*)varmodL,deta_CL,aveCell,varFluxL
!                !end if
!                !write(*,"(5F15.10)")varFluxL,varFluxR,aveCell,aveNearCell(2),aveNearCell(4)
!                if(varmodL /= deta_CL .or. varmodR /= deta_CR)then
!                    Beta = 1
!                    if(detect_type == ByDim )then
!                        cellset(index).Beta_line(i,1)=1
!                        cellset(index).Beta = 1
!                        Beta_line_x = 1
!                    elseif(detect_type == ByCell)then
!                        cellset(index).Beta_line = 1
!                        cellset(index).Beta = 1
!                        if(buffer_cell_switch == buffer_cell_no)then
!                            exit loop1    
!                        end if
!                    end if       
!                end if     
!            
!                !y����
!                dh2 = (2.0_prec*cellset(index).MJacobi(1,i,4))**2
!                !write(*,*)2.0_prec*cellset(index).MJacobi(1,i,4)
!                aveCell = Gauss_integral_SPs(varOriCell(:,i,m))/2.0_prec
!                varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(:,i,m),nsp)
!                varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(:,i,m),nsp)        
!            
!                !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell3 = cellset(index).nearcells(1)
!                do j = 1,4
!                    if(cellset(indexNearCell3).nearcells(j)==index)then
!                        sideth3 = j
!                    end if
!                end do
!                if((sideth3==2).or.(sideth3==1))then
!                    near_i3 = nsp+1-i
!                else
!                    near_i3 = i
!                end if
!                if(sideth3==2.OR.sideth3==4)then
!                    aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,near_i3,:,m))/2.0_prec
!                else
!                    aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,:,near_i3,m))/2.0_prec
!                end if
!                !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell4 = cellset(index).nearcells(3)
!                do j = 1,4
!                    if(cellset(indexNearCell4).nearcells(j)==index)then
!                        sideth4 = j
!                    end if
!                end do
!                if((sideth4==3).or.(sideth4==4))then
!                    near_i4 = nsp+1-i
!                else
!                    near_i4 = i
!                end if
!                if(sideth4==2.OR.sideth4==4)then
!                    aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,near_i4,:,m))/2.0_prec
!                else
!                    aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,:,near_i4,m))/2.0_prec
!                end if
! 
!                deta_L = aveCell-aveNearCell(1)
!                deta_R = aveNearCell(3)-aveCell
!                deta_CL = aveCell-varFluxL
!                deta_CR = varFluxR-aveCell
!            
!                varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh2)           
!                varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh2)
!            
!                if(varmodL /= deta_CL .or. varmodR /= deta_CR)then
!                    Beta = 1               
!                    !cellset(index).Beta_line(i,2)=1
!                    if(detect_type == ByDim )then
!                        cellset(index).Beta_line(i,2)=1
!                        Beta_line_y = 1
!                    elseif(detect_type == ByCell)then
!                        cellset(index).Beta_line = 1
!                        cellset(index).Beta = 1
!                        if(buffer_cell_switch == buffer_cell_no)then                        
!                            exit loop1         
!                        end if
!                    end if
!                end if  
!                !if(cellset(index).Beta == 1.or.Beta==1)then
!                !    cellset(index).Beta = 1
!                !end if
!                !!---���ɵ�Ԫ---------------------------------------------------------------------------------
!                !if(buffer_cell_switch == buffer_cell_yes .and. detect_type == ByDim.and.Beta==1)then
!                !
!                !    !x����
!                !    if(Beta_line_x == 1)then
!                !        if(sideth1==2.or.sideth1==4)then
!                !            cellset(indexNearCell1).Beta_line(near_i1,1) = 1
!                !        elseif(sideth1==1.or.sideth1==3)then
!                !            cellset(indexNearCell1).Beta_line(near_i1,2) = 1
!                !        end if
!                !        if(sideth2==2.or.sideth2==4)then
!                !            cellset(indexNearCell2).Beta_line(near_i2,1) = 1
!                !        elseif(sideth2==1.or.sideth2==3)then
!                !            cellset(indexNearCell2).Beta_line(near_i2,2) = 1
!                !        end if
!                !    end if
!                !    !y����
!                !    if(Beta_line_y == 1)then
!                !        if(sideth3==2.or.sideth3==4)then
!                !            cellset(indexNearCell3).Beta_line(near_i3,1) = 1
!                !        elseif(sideth3==1.or.sideth3==3)then
!                !            cellset(indexNearCell3).Beta_line(near_i3,2) = 1
!                !        end if
!                !        if(sideth4==2.or.sideth4==4)then
!                !            cellset(indexNearCell4).Beta_line(near_i4,1) = 1
!                !        elseif(sideth4==1.or.sideth4==3)then
!                !            cellset(indexNearCell4).Beta_line(near_i4,2) = 1
!                !        end if
!                !    end if
!                !    if(Beta==1)then
!                !        cellset(indexNearCell1).Beta = 1
!                !        cellset(indexNearCell2).Beta = 1
!                !        cellset(indexNearCell3).Beta = 1
!                !        cellset(indexNearCell4).Beta = 1
!                !    end if
!                !elseif(buffer_cell_switch == buffer_cell_yes .and. detect_type == ByCell)then
!                !    if(Beta==1)then
!                !        cellset(indexNearCell1).Beta_line = 1
!                !        cellset(indexNearCell2).Beta_line = 1
!                !        cellset(indexNearCell3).Beta_line = 1
!                !        cellset(indexNearCell4).Beta_line = 1
!                !
!                !        cellset(indexNearCell1).Beta = 1
!                !        cellset(indexNearCell2).Beta = 1
!                !        cellset(indexNearCell3).Beta = 1
!                !        cellset(indexNearCell4).Beta = 1
!                !        exit loop1
!                !    end if
!                !end if
!            end do loop2
!        end do loop1
!    end do
!end subroutine indicator_TVB
!
!function TVB_minmod(a1,a2,a3,dh2)
!    use real_precision
!    use global_var,only : TVB_M
!    implicit none
!    real(prec),external :: m_func
!    real(prec) :: TVB_minmod,a1,a2,a3,dh2
!    
!    if(abs(a1) <= TVB_M*dh2)then
!        TVB_minmod = a1
!    else
!        TVB_minmod = m_func(a1,a2,a3)
!    endif
!
!end function TVB_minmod
!
!function m_func(a1,a2,a3)
!    !m(a1,a2,a3,...,ak)
!    use real_precision 
!    implicit none
!    real(prec),external::min_abs
!    real(prec)::a1,a2,a3,m_func
!    
!    !��ͬ����
!    if((a1*a2>0).and.(a1*a3>0).and.(a1*a2*a3>0))then
!        m_func = min_abs(a1,a2,a3)
!    elseif((a1*a2>0).and.(a1*a3>0).and.(a1*a2*a3<0))then
!        m_func = -1.0_prec*min_abs(a1,a2,a3)
!    else
!        m_func = 0.0_prec
!    end if
!    return
!end function
!
!function min_abs(a,b,c)
!    !ȡ����ֵ��Сֵ
!    use real_precision
!    implicit none
!    real(prec)::a,b,c,min_abs
!    if(abs(a)<=abs(b).and.abs(a)<=abs(c))then
!        min_abs=abs(a)
!    elseif(abs(b)<=abs(a).and.abs(b)<=abs(c))then
!        min_abs=abs(b)
!    else
!        min_abs=abs(c)
!    end if
!    return
!end function
!!!!--------------------------------------------------------------------------------------------------------------------
!
!!!!---KXRCF-----------------------------------------------------------------------------------------------------------------
!subroutine indicator_KXRCF()
!    !Refs:
!    !   [2]	Krivodonova L, Xin J, Remacle J F, et al. Shock detection and limiting with discontinuous Galerkin methods for hyperbolic conservation laws [J]. Applied Numerical Mathematics, 2004, 48(3-4): 323-38.
!
!
!
!end subroutine indicator_KXRCF
!!!!--------------------------------------------------------------------------------------------------------------------
!
!!!!---!modal energy-----------------------------------------------------------------------------------------------------------------
!subroutine indicator_MDH!modal energy
!    !Refs:
!    !   [3]	Hennemann S, Rueda-Ram��rez A M, Hindenlang F J, et al. A provably entropy stable subcell shock capturing approach for high order split form DG for the compressible Euler equations [J]. Journal of Computational Physics, 2021, 426(109935).    use parameter_setting
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m,index,Beta,Beta_line_x,Beta_line_y
!    real(prec) :: sum_M1,sum_M2,E_Le ,T_Le
!    real(prec),dimension(:) :: M_Le(nsp)  
!
!    real(prec),dimension(nsp,nsp,4) :: varOriCell
!                
!    !K_La_Le=reshape( (/  -0.44542134632044532872_prec,  0.44856404312571887327_prec,    1.4079281687625479597_prec,     0.44856404312571887327_prec,    -0.44542134632044532872_prec,&
!    !                     -0.26295072534780090115_prec,  -0.31564963775813707471_prec,   0.0_prec,                       0.31564963775813707471_prec,    0.26295072534780090115_prec, &
!    !                     0.27412134137104513456_prec,   -0.04924826331462704873_prec,   -0.44974615611283617166_prec,   -0.04924826331462704873_prec,   0.27412134137104513456_prec, &
!    !                     -0.222081873567294817185_prec, 0.37373739635316142603_prec,    0.0_prec,                       -0.37373739635316142603_prec,   0.222081873567294817185_prec,&
!    !                     0.123506106334137037799_prec,  -0.349780276313832245607_prec,  0.452548339959390415617_prec,   -0.349780276313832245607_prec,  0.123506106334137037799_prec/),   &
!    !                    shape(K_La_Le), order=(/2,1/) )  
!    
!    T_Le = 0.5_prec*10**(-1.8*(nsp+1)**0.25)
!    do i  = 1,ncells
!        varOriCell(:,:,1) = cellset(i).spvalue_ori(:,:,1)*cellset(i).spvalue_ori(:,:,4)
!        do m = 1,1
!            
!            do k = 1,nsp
!                !x direction :varOriCell(k,:,m)
!                do l = 1,nsp
!                    M_Le(l) = dot_product(K_La_Le(l,:),varOriCell(k,:,m))                      
!                end do
!                !write(*,"(2I4 ,5F20.10)")i,k,M_Le
!                sum_M1 = 0.0_prec
!                sum_M2 = 0.0_prec
!                do j = 1,nsp
!                    sum_M1 = sum_M1+M_Le(j)**2
!                end do
!                do j = 1,nsp-1
!                    sum_M2 = sum_M2+M_Le(j)**2
!                end do
!                E_Le = max(M_Le(nsp)**2/sum_M1,M_Le(nsp-1)**2/sum_M2)
!                
!                if(E_Le > T_Le)then                
!                    !write(*,*) i,k,E_Le,T_Le
!                    cellset(i).Beta_line(k,1) = 1
!                    cellset(i).Beta = 1     
!                end if                     
!                
!                !y direction :varOriCell(:,k,m)
!                do l = 1,nsp
!                    M_Le(l) = dot_product(K_La_Le(l,:),varOriCell(:,k,m))                      
!                end do
!                !write(*,"(2I4 ,5F20.10)")i,k,M_Le
!                sum_M1 = 0.0_prec
!                sum_M2 = 0.0_prec
!                do j = 1,nsp
!                    sum_M1 = sum_M1+M_Le(j)**2
!                end do
!                do j = 1,nsp-1
!                    sum_M2 = sum_M2+M_Le(j)**2
!                end do
!                E_Le = max(M_Le(nsp)**2/sum_M1,M_Le(nsp-1)**2/sum_M2)
!                if(E_Le > T_Le)then                
!                    !write(*,*) i,k,E_Le,T_Le
!                    cellset(i).Beta_line(k,2) = 1
!                    cellset(i).Beta = 1    
!                end if           
!            end do              
!        end do
!    end do                                                                                
!
!end subroutine indicator_MDH
!
!
!
!subroutine indicator_ATV
!    !Refs:[1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29.
!    !   ʵ���ϣ���ATV����С�ĸĶ���ԭ�����ǲ�ȡAverage Total Variation,�������������ء��ڴ˽�ƽ���ܱ���Ϊ����
!    !   ������ˣ�Ҳ�������⣬�Ǿ��ǹ⻬����������ⵥԪ������ǿ�������ڵ�����£����ܶ��������ⲻ������
!    !   ������֤������1����ȡȫ��������ֵ*theta��Ϊ�ж����ݣ�2����ȡ��ȫ��������ֵ-��Сֵ��*theta��Ϊ�ж����ݡ�
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m,near_k,ikesi1,ikesi2,ieta1,ieta2,ikesi,ieta
!    integer :: sideth(4),near_kk(4),sidethL,sidethR,indexNearCell,indexNearCell1, indexNearCell2, indexNearCell3, indexNearCell4 
!    integer :: index,indexCellL,indexCellR,Beta,Beta_line_x,Beta_line_y
!    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh2
!    real(prec),dimension(4) :: aveNearCell,TV,ATV,Max_TV,Max_v,Min_v
!    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
!    real(prec),dimension(nsp,nsp,4) :: varOriCell,varOriCellL,varOriCellR
!    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
!    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa,temp,theta,Vtemp_x(4),Vtemp_y(4)
!    theta = 0.1_prec
!    ! �ȼ���TV��the total variation of the solution��.��������֮�䶼��Ҫ���㡣
!    TV = 0.0_prec
!    Max_TV = 0.0_prec
!    Max_v = 0.0_prec
!    Min_v = 0.0_prec
!    !���ȼ���һ����Ԫ�ڲ������TV
!    do i = 1,ncells               
!        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)!!����ȡ����Ԫ -- ������ԭʼ����             
!        do m = 1,1  !������
!            Max_v = MaxVal(varOriCell(:,:,m))
!            Min_v = MinVal(varOriCell(:,:,m))
!            do k = 1,nsp
!                do l = 1,nsp-1
!                    temp = abs(varOriCell(k,l+1,m) - varOriCell(k,l,m)) !x direction
!                    TV(m) = TV(m) + temp
!                    Max_TV(m) = max(Max_TV(m),temp)
!                    temp = abs(varOriCell(l+1,k,m) - varOriCell(l,k,m)) !y
!                    TV(m) = TV(m) + temp
!                    Max_TV(m) = max(Max_TV(m),temp)
!                end do
!            end do      
!        end do
!    end do        
!    !write(*,*)'0000'
!    !Ȼ����ݱ߼���ߵ����ڵ�Ԫ�ڽ��ӱߵ��������
!    do i = 1,nsides        
!        indexCellL = sideset(i).nearcells(1)
!        indexCellR = sideset(i).nearcells(2)
!        if(indexCellL == 0 .OR. indexCellR == 0)cycle
!        varOriCellL(:,:,:) = cellset(indexCellL).spvalue_ori(:,:,:)
!        varOriCellR(:,:,:) = cellset(indexCellR).spvalue_ori(:,:,:)
!        do j = 1,4
!            if(i == cellset(indexCellL).sides(j)) sidethL = j            !���������ڵ�Ԫ�ĵڼ����
!            if(i == cellset(indexCellR).sides(j)) sidethR = j
!        end do
!                
!        do m = 1,1
!            do k = 1,nsp
!                if((sidethL*sidethR==2).or.(sidethL*sidethR==12).or.(sidethL==sidethR))then!������������������������෴
!                    near_k = nsp+1-k
!                else
!                    near_k = k
!                end if
!                if(sidethL==1)then
!                    ikesi1 = 1;     ieta1 = k;
!                elseif(sidethL==2)then
!                    ikesi1 = k;     ieta1 = nsp
!                elseif(sidethL==3)then
!                    ikesi1 = nsp;   ieta1 = k;
!                elseif(sidethL==4)then
!                    ikesi1 = k;     ieta1 = 1
!                end if
!                if(sidethR==1)then
!                    ikesi2 = 1;     ieta2 = near_k
!                elseif(sidethR==2)then
!                    ikesi2 = near_k;ieta2 = nsp
!                elseif(sidethR==3)then
!                    ikesi2 = nsp;   ieta2 = near_k
!                elseif(sidethR==4)then
!                    ikesi2 = near_k;ieta2 = 1
!                end if
!                temp = abs(varOriCellL(ikesi1,ieta1,m)-varOriCellR(ikesi2,ieta2,m))
!                TV(m) = TV(m) + temp
!                Max_TV(m) = max(Max_TV(m),temp)
!            end do
!        end do
!    end do
!    ATV = TV/(ncells*nsp*(nsp-1)*2+(nsides-nbdsides)*nsp)
!    !write(*,*)ATV(1),Max_TV(1)
!    !-------------------------------------------------------
!    !
!    !-------------------------------------------------------
!    !write(*,*)'0001'
!    !!����Ԫ -- ������ԭʼ����
!    do i = 1,ncells
!        Vtemp_x = 0.0_prec
!        Vtemp_y = 0.0_prec
!        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)    
!        !varOriCell(:,:,1) = varOriCell(:,:,1)*varOriCell(:,:,4) !rho*p
!        !varOriCell(:,:,1) = varOriCell(:,:,4)/varOriCell(:,:,1) !p/rho
!        !�ڵ�Ԫ
!        do j = 1,4  !4�����ڵ�Ԫ
!            indexNearCell  = cellset(i).nearcells(j)!���ڵ�Ԫ������
!            !write(*,*)j,indexNearCell,cellset(indexNearCell).nearcells
!            do k = 1,4
!                if(i==cellset(indexNearCell).nearcells(k))then
!                    sideth(j) = k   !��¼����Ԫ�Ĳ�������ڵ�Ԫ �ĵڼ����
!                end if
!            end do
!        
!            varOriCellNearCell(j,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,1)*varOriCellNearCell(i,:,:,4) !rho*p
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,4)/varOriCellNearCell(i,:,:,1) !p/rho
!        end do 
!        !write(*,*)i,sideth
!        !----x direction --------------------------------------------------------------------
!        beta = 0
!        do m = 1,1
!            !write(*,*)'12'
!            do k =1,nsp
!                Beta_line_x = 0
!                Beta_line_y = 0
!                if(sideth(4)==1)then    !x���򣬿������4�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                    ikesi = 1;          ieta = k;
!                elseif(sideth(4)==2)then
!                    ikesi = k;    ieta = nsp
!                elseif(sideth(4)==3)then
!                    ikesi = nsp;   ieta = nsp+1-k;
!                elseif(sideth(4)==4)then
!                    ikesi = nsp+1-k;     ieta = 1
!                end if  
!                !write(*,*)'13',ikesi,ieta
!                
!                temp =  abs(varOriCell(k,1,m)-varOriCellNearCell(4,ikesi,ieta,m))
!                Vtemp_x(m) = max(Vtemp_x(m),temp)                
!                
!                do l = 2,nsp
!                    temp =  abs(varOriCell(k,l,m)-varOriCell(k,l-1,m))
!                    Vtemp_x(m) = max(Vtemp_x(m),temp)
!                    if(temp>theta*Max_TV(m))then    !MV : theta*Max_TV(m);  ATV : theta*ATV(m)
!                        cellset(i).Beta_line(k,1) = 1
!                        cellset(i).Beta = 1
!                    end if
!                    !write(*,*)'2',temp,theta*Max_TV(m)
!                end do 
!                !x���򣬿������2�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                if(sideth(2)==1)then    
!                    ikesi = 1;          ieta = nsp+1-k;
!                elseif(sideth(2)==2)then
!                    ikesi = nsp+1-k;    ieta = nsp
!                elseif(sideth(2)==3)then
!                    ikesi = nsp;   ieta = k;
!                elseif(sideth(2)==4)then
!                    ikesi = k;     ieta = 1
!                end if  
!                temp =  abs(varOriCell(k,nsp,m)-varOriCellNearCell(2,ikesi,ieta,m))
!                Vtemp_x(m) = max(Vtemp_x(m),temp)
!
!            !--------------------------------------------------------------------------------------------------       
!    
!            !----y direction ----------------------------------------------------------------------------------
!                if(sideth(1)==1)then    !x���򣬿������1�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                    ikesi = 1;          ieta = nsp+1-k;
!                elseif(sideth(1)==2)then
!                    ikesi = nsp+1-k;    ieta = nsp
!                elseif(sideth(1)==3)then
!                    ikesi = nsp;   ieta = k;
!                elseif(sideth(1)==4)then
!                    ikesi = k;     ieta = 1
!                end if  
!                temp =  abs(varOriCell(1,k,m)-varOriCellNearCell(1,ikesi,ieta,m))
!                Vtemp_y(m) = max(Vtemp_y(m),temp)
!        
!                do l = 2,nsp
!                    temp =  abs(varOriCell(l,k,m)-varOriCell(l-1,k,m))
!                    Vtemp_y(m) = max(Vtemp_y(m),temp)
!                end do 
!                !x���򣬿������3�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                if(sideth(3)==1)then    
!                    ikesi = 1;          ieta = k;
!                elseif(sideth(3)==2)then
!                    ikesi = k;    ieta = nsp
!                elseif(sideth(3)==3)then
!                    ikesi = nsp;   ieta = nsp+1-k;
!                elseif(sideth(3)==4)then
!                    ikesi = nsp+1-k;     ieta = 1
!                end if  
!                temp =  abs(varOriCell(nsp,k,m)-varOriCellNearCell(3,ikesi,ieta,m))
!                Vtemp_y(m) = max(Vtemp_y(m),temp)
!                
!                if(Vtemp_x(m)>theta*Max_TV(m))then
!                    cellset(i).Beta_line(k,1) = 1
!                    cellset(i).Beta = 1
!                    Beta_line_x = 1
!                    beta = 1
!                end if
!                if(Vtemp_y(m)>theta*Max_TV(m))then
!                    cellset(i).Beta_line(k,2) = 1
!                    cellset(i).Beta = 1
!                    Beta_line_y = 1
!                    beta = 1
!                end if
!                
!                !���ɵ�Ԫ
!                do j = 1,4
!                    if(sideth(j)*j==2.or.sideth(j)*j==12.or.sideth(j)==j)then
!                        near_kk(j) = nsp+1-k
!                    else 
!                        near_kk(j) = k
!                    end if
!                end do
!                !write(*,*)near_kk
!                !if(buffer_cell_switch == buffer_cell_yes .and. detect_type == ByDim.and.Beta==1)then
!                !    !x����
!                !    if(Beta_line_x == 1)then
!                !        if(sideth(4)==2.or.sideth(4)==4)then
!                !            cellset(cellset(i).nearcells(4)).Beta_line(near_kk(4),1) = 1
!                !        elseif(sideth(4)==1.or.sideth(4)==3)then
!                !            cellset(cellset(i).nearcells(4)).Beta_line(near_kk(4),2) = 1
!                !        end if
!                !        if(sideth(2)==2.or.sideth(2)==4)then
!                !            cellset(cellset(i).nearcells(2)).Beta_line(near_kk(2),1) = 1
!                !        elseif(sideth(2)==1.or.sideth(2)==3)then
!                !            cellset(cellset(i).nearcells(2)).Beta_line(near_kk(2),2) = 1
!                !        end if
!                !    end if
!                !    !y����
!                !    if(Beta_line_y == 1)then
!                !        if(sideth(1)==2.or.sideth(1)==4)then
!                !            cellset(cellset(i).nearcells(1)).Beta_line(near_kk(1),1) = 1
!                !        elseif(sideth(1)==1.or.sideth(1)==3)then
!                !            cellset(cellset(i).nearcells(1)).Beta_line(near_kk(1),2) = 1
!                !        end if
!                !        if(sideth(3)==2.or.sideth(3)==4)then
!                !            cellset(cellset(i).nearcells(3)).Beta_line(near_kk(3),1) = 1
!                !        elseif(sideth(3)==1.or.sideth(3)==3)then
!                !            cellset(cellset(i).nearcells(3)).Beta_line(near_kk(3),2) = 1
!                !        end if              
!                !    end if
!                !    if(Beta==1)then
!                !        do j = 1,4
!                !            cellset(cellset(i).nearcells(j)).Beta = 1
!                !        end do
!                !    end if
!                !elseif(buffer_cell_switch == buffer_cell_yes .and. detect_type == ByCell)then
!                !    if(Beta==1)then
!                !        do j = 1,4
!                !            cellset(cellset(i).nearcells(j)).Beta_line = 1
!                !            cellset(cellset(i).nearcells(j)).Beta = 1
!                !        end do
!                !    end if
!                !end if
!            end do     
!        end do        
!    end do
!    
!
!end subroutine indicator_ATV
!
!subroutine add_buffer_cell
!    !���Ǵ����ظ��ԣ��ڴ˽�������ⵥԪ����д���ӳ���
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer,dimension(:),allocatable :: Cell_Beta
!    integer,dimension(:,:,:),allocatable :: Cell_Beta_line
!    integer :: i,j,k,l,m,indexNearCell
!    integer :: sideth(4),near_kk(4)
!    
!    allocate(Cell_Beta(ncells),Cell_Beta_line(ncells,nsp,2))
!    do i = 1,ncells
!        Cell_Beta(i) = cellset(i).Beta
!    end do
!    
!    if(detect_type == ByCell)then
!        do i = 1,ncells
!            if(Cell_Beta(i)==0)cycle    !�����������ⵥԪ��ֱ������
!            do j = 1,4
!                cellset(cellset(i).nearcells(j)).Beta = 1           
!            end do
!            !write(*,*)i
!        end do
!    elseif(detect_type == ByDim)then
!        do i = 1,ncells         
!            if(Cell_Beta(i)==0)cycle    !�����������ⵥԪ��ֱ������
!            
!            !-----��������i�����ⵥԪ�����         
!            do j = 1,4  !4�����ڵ�Ԫ
!                indexNearCell  = cellset(i).nearcells(j)!���ڵ�Ԫ������
!                cellset(indexNearCell).Beta = 1
!                do k = 1,4
!                    if(i==cellset(indexNearCell).nearcells(k))then
!                        sideth(j) = k   !��¼����Ԫ�Ĳ�������ڵ�Ԫ �ĵڼ����
!                    end if
!                end do        
!            end do 
!            do k =1,nsp
!                do j = 1,4
!                    if(sideth(j)*j==2.or.sideth(j)*j==12.or.sideth(j)==j)then
!                        near_kk(j) = nsp+1-k    !�ж��ڵ�Ԫ�ǵڼ���
!                    else 
!                        near_kk(j) = k
!                    end if
!                end do
!                !x����
!                if(Cell_Beta_line(i,k,1) == 1)then
!                    if(sideth(4)==2.or.sideth(4)==4)then
!                        cellset(cellset(i).nearcells(4)).Beta_line(near_kk(4),1) = 1
!                    elseif(sideth(4)==1.or.sideth(4)==3)then
!                        cellset(cellset(i).nearcells(4)).Beta_line(near_kk(4),2) = 1
!                    end if
!                    if(sideth(2)==2.or.sideth(2)==4)then
!                        cellset(cellset(i).nearcells(2)).Beta_line(near_kk(2),1) = 1
!                    elseif(sideth(2)==1.or.sideth(2)==3)then
!                        cellset(cellset(i).nearcells(2)).Beta_line(near_kk(2),2) = 1
!                    end if
!                end if
!                !y����
!                if(Cell_Beta_line(i,k,2) == 1)then
!                    if(sideth(1)==2.or.sideth(1)==4)then
!                        cellset(cellset(i).nearcells(1)).Beta_line(near_kk(1),1) = 1
!                    elseif(sideth(1)==1.or.sideth(1)==3)then
!                        cellset(cellset(i).nearcells(1)).Beta_line(near_kk(1),2) = 1
!                    end if
!                    if(sideth(3)==2.or.sideth(3)==4)then
!                        cellset(cellset(i).nearcells(3)).Beta_line(near_kk(3),1) = 1
!                    elseif(sideth(3)==1.or.sideth(3)==3)then
!                        cellset(cellset(i).nearcells(3)).Beta_line(near_kk(3),2) = 1
!                    end if              
!                end if                
!            end do
!        end do   
!    end if
!   
!    deallocate(Cell_Beta,Cell_Beta_line)
!end subroutine add_buffer_cell
!

!
!void RoeSolver::RoeAverage(
!        double rhoL, double rhouL, double rhovL, double rhowL, double EL,
!        double rhoR, double rhouR, double rhovR, double rhowR, double ER,
!        double &uRoe, double &vRoe, double &wRoe, double &hRoe, double &URoe,
!        double &cRoe, double &ocRoe)
!    {        
!        static NekDouble gamma = m_params["gamma"]();
!        
!        // Left and right velocities
!        NekDouble orhoL = 1.0/rhoL;
!        NekDouble orhoR = 1.0/rhoR;
!        NekDouble uL = rhouL * orhoL;
!        NekDouble vL = rhovL * orhoL;
!        NekDouble wL = rhowL * orhoL;
!        NekDouble uR = rhouR * orhoR;
!        NekDouble vR = rhovR * orhoR;
!        NekDouble wR = rhowR * orhoR;
!
!        // Left and right pressures
!        NekDouble pL = (gamma - 1.0) *
!            (EL - 0.5 * (rhouL * uL + rhovL * vL + rhowL * wL));
!        NekDouble pR = (gamma - 1.0) *
!            (ER - 0.5 * (rhouR * uR + rhovR * vR + rhowR * wR));
!        
!        // Left and right enthalpy
!        NekDouble hL = (EL + pL) * orhoL;
!        NekDouble hR = (ER + pR) * orhoR;
!
!        // Square root of rhoL and rhoR.
!        NekDouble srL  = sqrt(rhoL);
!        NekDouble srR  = sqrt(rhoR);
!        NekDouble srLR = srL + srR;
!        NekDouble osrLR = 1.0/srLR;
!        
!        // Velocity, enthalpy and sound speed Roe averages (equation 11.60).
!        uRoe   = (srL * uL + srR * uR) * osrLR;
!        vRoe   = (srL * vL + srR * vR) * osrLR;
!        wRoe   = (srL * wL + srR * wR) * osrLR;
!        hRoe   = (srL * hL + srR * hR) * osrLR;
!        URoe   = (uRoe * uRoe + vRoe * vRoe + wRoe * wRoe);
!        cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe));
!        ocRoe  = 1.0/cRoe;
!    }
!
!if(case_comp == DoubleMach_case)then
!        ˫��շ��䣬��ȡ�޴�͸���ȱ߽粿����Ҫÿʱ�䲽���¸�ֵ
!        do i = 1,nbdsides
!            do j = 1,4               
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                if(indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then!�±߽�
!                    ���ж�λ�ã�ȡ�߽��
!                    p1 = j
!                    p2 = p1+1
!                    if(p1 == 4)p2=1
!                    if(xy_coor(cellset(indexCell).nodes(p1),1)>1.0_prec/6.0_prec .OR. xy_coor(cellset(indexCell).nodes(p2),1)>1.0_prec/6.0_prec)then                       
!                        �߽�Ϊ�޴�͸�������������öԳ��档���±߽�Գ�
!                        if(j == 1)then
!                            cellset(indexNearCell) = cellset(cellset(indexNearCell).index)   
!                            index = cellset(indexNearCell).index
!                            cellset(indexNearCell).Beta = cellset(index).Beta
!                            do k = 1,4
!                                cellset(indexNearCell).nodes(k) = cellset(index).nodes(5-k) 
!                            end do
!                            cellset(indexNearCell).sides(1) = cellset(index).sides(3) 
!                            cellset(indexNearCell).sides(3) = cellset(index).sides(1) 
!                            do k = 1,nsp
!                                do l = 1,nsp
!                                    cellset(indexNearCell).det_J(l,k) = cellset(index).det_J(nsp+1-l,k)
!                                    cellset(indexNearCell).MJacobi(l,k,:) = cellset(index).MJacobi(nsp+1-l,k,:) 
!                                    cellset(indexNearCell).Mdirect(l,k,:) = cellset(index).Mdirect(nsp+1-l,k,:) 
!                                    cellset(indexNearCell).sp_coor(l,k,:) = cellset(index).sp_coor(nsp+1-l,k,:) 
!                                    cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(index).spvalue_ori(nsp+1-l,k,:) 
!                                end do                         
!                            end do
!                            cellset(indexNearCell).MJacobi(:,:,2:3) = -cellset(indexNearCell).MJacobi(:,:,2:3) 
!                            cellset(indexNearCell).Mdirect(:,:,2:3) = -cellset(indexNearCell).Mdirect(:,:,2:3) 
!                            
!                            cellset(indexNearCell).spvalue_ori(:,:,3) = -cellset(indexNearCell).spvalue_ori(:,:,3)
!                            cellset(indexNearCell).spvalue_con(:,:,3) = -cellset(indexNearCell).spvalue_con(:,:,3)
!                        elseif(j == 2)then
!                            
!                        elseif(j == 3)then
!                            
!                        elseif(j == 4)then
!                            
!                        end if
!                    else   
!                        cellset(indexNearCell) = cellset(cellset(indexNearCell).index)  !�±߽�ǰ1/6��       
!                    end if
!                else if(indexNearCell > ncells+n_BCells_D)then!���ಿ�ֿ�������һ�㣬�����ϱ߽��漤��λ�ñ仯�ı߽�����Ҳ��������һ��   
!                    cellset(indexNearCell) = cellset(cellset(indexNearCell).index)                        
!                end if         
!            end do        
!        end do
    !elseif(case_comp == Riemann2D_case)then
    !!n_BCells_D+n_BCells_R+1,n_BCells_D+n_BCells_R+n_BCells_U
    !    do i = ncells+1,ncells+nbdsides
    !        !allocate(cellset(i).spvalue_ori(nsp,nsp,4),cellset(i).spvalue_con(nsp,nsp,4))
    !        cellset(i) = cellset(cellset(i).index)    !��ֵ
    !    end do
    !    !Down
    !    do i = 1,n_BCells_D
    !        !write(*,*)5551
    !        indexCell = BoundCells_index_set(i)
    !        do j = 1,4                               
    !            indexNearCell = cellset(indexCell).nearcells(j)
    !            if(indexNearCell > ncells .and. indexNearCell < ncells + n_BCells_D +1)then!�߽�
    !                th1 = j
    !                exit
    !            end if
    !        end do
    !       
    !        if(th1 == 1)then     
    !            !write(*,*) th1,indexCell,indexNearCell
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(l,:,:) = cellset(indexCell).spvalue_ori(1,:,:)                             
    !                cellset(indexNearCell).spvalue_con(l,:,:) = cellset(indexCell).spvalue_con(1,:,:) 
    !                !cellset(indexNearCell).spvalue_con(l,:,:) = cellset(indexCell).spvalue_con(1,:,:)
    !            end do                                 
    !        elseif(th1 == 2)then
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(l,:,:) = cellset(indexCell).spvalue_ori(:,nsp,:)                             
    !                cellset(indexNearCell).spvalue_con(l,:,:) = cellset(indexCell).spvalue_con(:,nsp,:) 
    !            end do 
    !        elseif(th1 == 3)then
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp,nsp+1-k,:)                             
    !                    cellset(indexNearCell).spvalue_con(l,k,:) = cellset(indexCell).spvalue_con(nsp,nsp+1-k,:) 
    !                end do
    !            end do 
    !        elseif(th1 == 4)then
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-k,1,:)                             
    !                    cellset(indexNearCell).spvalue_con(l,k,:) = cellset(indexCell).spvalue_con(nsp+1-k,1,:) 
    !                end do
    !            end do 
    !        end if
    !        !write(*,*)5552
    !    end do
    !    
    !    !Right
    !    do i = n_BCells_D+1,n_BCells_D+n_BCells_R
    !        indexCell = BoundCells_index_set(i)
    !        do j = 1,4                               
    !            indexNearCell = cellset(indexCell).nearcells(j)
    !            if(indexNearCell > ncells+n_BCells_D .and. indexNearCell < ncells+n_BCells_D+n_BCells_R+1 )then!�߽�
    !                th1 = j
    !                exit
    !            end if
    !        end do           
    !        if(th1 == 1)then           
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(:,l,:) = cellset(indexCell).spvalue_ori(1,:,:)                             
    !                cellset(indexNearCell).spvalue_con(:,l,:) = cellset(indexCell).spvalue_con(1,:,:) 
    !            end do                                 
    !        elseif(th1 == 2)then
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(:,1,:) = cellset(indexCell).spvalue_ori(:,nsp,:)                             
    !                cellset(indexNearCell).spvalue_con(:,1,:) = cellset(indexCell).spvalue_con(:,nsp,:) 
    !            end do 
    !        elseif(th1 == 3)then
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp,nsp+1-k,:)                             
    !                    cellset(indexNearCell).spvalue_con(k,l,:) = cellset(indexCell).spvalue_con(nsp,nsp+1-k,:) 
    !                end do
    !            end do 
    !        elseif(th1 == 4)then
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-k,1,:)                             
    !                    cellset(indexNearCell).spvalue_con(k,l,:) = cellset(indexCell).spvalue_con(nsp+1-k,1,:) 
    !                end do
    !            end do 
    !        end if
    !    end do   
    !    !Up
    !    do i = n_BCells_D+n_BCells_R+1,n_BCells_D+n_BCells_R+n_BCells_U
    !        indexCell = BoundCells_index_set(i)
    !        do j = 1,4                               
    !            indexNearCell = cellset(indexCell).nearcells(j)
    !            if(indexNearCell > ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells+n_BCells_D+n_BCells_R+n_BCells_U+1 )then!�߽�
    !                th1 = j
    !                exit
    !            end if
    !        end do
    !       
    !        if(th1 == 1)then           
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(1,nsp+1-k,:)                             
    !                    cellset(indexNearCell).spvalue_con(l,k,:) = cellset(indexCell).spvalue_con(1,nsp+1-k,:) 
    !                end do
    !            end do                                 
    !        elseif(th1 == 2)then
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp,:)                             
    !                    cellset(indexNearCell).spvalue_con(l,k,:) = cellset(indexCell).spvalue_con(nsp+1-k,nsp,:) 
    !                end do
    !            end do 
    !        elseif(th1 == 3)then
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(l,:,:) = cellset(indexCell).spvalue_ori(nsp,:,:)                             
    !                cellset(indexNearCell).spvalue_con(l,:,:) = cellset(indexCell).spvalue_con(nsp,:,:) 
    !            end do 
    !        elseif(th1 == 4)then
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(l,:,:) = cellset(indexCell).spvalue_ori(:,1,:)                             
    !                cellset(indexNearCell).spvalue_con(l,:,:) = cellset(indexCell).spvalue_con(:,1,:) 
    !            end do 
    !        end if
    !    end do 
    !    !Left
    !    do i = n_BCells_D+n_BCells_R+n_BCells_U+1,n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
    !        indexCell = BoundCells_index_set(i)
    !        do j = 1,4                               
    !            indexNearCell = cellset(indexCell).nearcells(j)
    !            if(indexNearCell > ncells+n_BCells_D+n_BCells_R+n_BCells_U .and. indexNearCell < ncells+n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L+1 )then!�߽�
    !                th1 = j
    !                exit
    !            end if
    !        end do
    !       
    !        if(th1 == 1)then           
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(1,nsp+1-k,:)                             
    !                    cellset(indexNearCell).spvalue_con(k,l,:) = cellset(indexCell).spvalue_con(1,nsp+1-k,:) 
    !                end do
    !            end do                                 
    !        elseif(th1 == 2)then
    !            do l = 1,nsp
    !                do k = 1,nsp
    !                    cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp,:)                             
    !                    cellset(indexNearCell).spvalue_con(k,l,:) = cellset(indexCell).spvalue_con(nsp+1-k,nsp,:) 
    !                end do
    !            end do 
    !        elseif(th1 == 3)then
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(:,l,:) = cellset(indexCell).spvalue_ori(nsp,:,:)                             
    !                cellset(indexNearCell).spvalue_con(:,l,:) = cellset(indexCell).spvalue_con(nsp,:,:) 
    !            end do 
    !        elseif(th1 == 4)then
    !            do l = 1,nsp
    !                cellset(indexNearCell).spvalue_ori(:,l,:) = cellset(indexCell).spvalue_ori(:,1,:)                             
    !                cellset(indexNearCell).spvalue_con(:,l,:) = cellset(indexCell).spvalue_con(:,1,:) 
    !            end do 
    !        end if
    !    end do 
        !stop

!if(j==1)then
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(4),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(1),2)-xy_coor(cellset(indexCell).nodes(4),2)
!                    cellset(indexNearCell).nodes(1) = pCounter
!                    pCounter = pCounter + 1
!                    
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(3),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(2),2)-xy_coor(cellset(indexCell).nodes(3),2)
!                    cellset(indexNearCell).nodes(2) = pCounter
!                    pCounter = pCounter + 1
!                    cellset(indexNearCell).nodes(3) = cellset(indexCell).nodes(1)
!                    cellset(indexNearCell).nodes(4) = cellset(indexCell).nodes(2)
!                    !write(*,*)cellset(indexCell).nodes
!                    !write(*,*)indexCell,cellset(indexNearCell).nodes
!                    !write(*,*)indexNearCell
!                    
!                elseif(j==2)then
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(1),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(2),2)-xy_coor(cellset(indexCell).nodes(1),2)
!                    cellset(indexNearCell).nodes(1) = pCounter
!                    pCounter = pCounter + 1
!                    
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(4),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(3),2)-xy_coor(cellset(indexCell).nodes(4),2)
!                    cellset(indexNearCell).nodes(2) = pCounter
!                    pCounter = pCounter + 1
!                    cellset(indexNearCell).nodes(3) = cellset(indexCell).nodes(3)
!                    cellset(indexNearCell).nodes(4) = cellset(indexCell).nodes(2)
!                elseif(j==3)then
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(2),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(3),2)-xy_coor(cellset(indexCell).nodes(2),2)
!                    cellset(indexNearCell).nodes(1) = pCounter
!                    pCounter = pCounter + 1
!                    
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(1),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(4),2)-xy_coor(cellset(indexCell).nodes(1),2)
!                    cellset(indexNearCell).nodes(2) = pCounter
!                    pCounter = pCounter + 1
!                    cellset(indexNearCell).nodes(3) = cellset(indexCell).nodes(4)
!                    cellset(indexNearCell).nodes(4) = cellset(indexCell).nodes(3)
!                elseif(j==4)then
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(2),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(3),2)-xy_coor(cellset(indexCell).nodes(2),2)
!                    cellset(indexNearCell).nodes(1) = pCounter
!                    pCounter = pCounter + 1
!                    
!                    xy_coor(pCounter,1) = xy_coor(cellset(indexCell).nodes(1),1)
!                    xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexCell).nodes(4),2)-xy_coor(cellset(indexCell).nodes(1),2)
!                    cellset(indexNearCell).nodes(2) = pCounter
!                    pCounter = pCounter + 1
!                    cellset(indexNearCell).nodes(3) = cellset(indexCell).nodes(4)
!                    cellset(indexNearCell).nodes(4) = cellset(indexCell).nodes(3)
!                end if



!
!subroutine indicator_TVB
!    !TVB��⣬�����غ����
!    ! ÿһά�����бȽϣ����ֲ�����minmod�����ı���Ϊ���ⵥԪ��ֹͣ����
!    !Refs:
!    !   [1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29. 
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m
!    integer :: sideth1,sideth2,sideth3,sideth4,near_i1,near_i2,near_i3,near_i4,indexNearCell,indexNearCell1, indexNearCell2, indexNearCell3, indexNearCell4 
!    integer :: index
!    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh2
!    real(prec),dimension(4) :: aveNearCell
!    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
!    real(prec),dimension(nsp,nsp,4) :: varOriCell
!    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
!    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa
!    do index = 1,ncells      
!    
!        !!����Ԫ -- ������ԭʼ����
!        varOriCell(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!
!        !varOriCell(:,:,1) = varOriCell(:,:,1)*varOriCell(:,:,4) !rho*p
!        !varOriCell(:,:,1) = varOriCell(:,:,4)/varOriCell(:,:,1) !p/rho
!        !�ڵ�Ԫ
!        do i = 1,4
!            indexNearCell  = cellset(index).nearcells(i)
!            varOriCellNearCell(i,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,1)*varOriCellNearCell(i,:,:,4) !rho*p
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,4)/varOriCellNearCell(i,:,:,1) !p/rho
!        end do 
!        !��ά���
!        !x����
!        loop1:do m = 1,1!������
!            !write(*,*)m
!            loop2:do i = 1,nsp!ÿ�зֱ����
!                dh2 = (2.0_prec*cellset(index).MJacobi(i,1,1))**2
!                !write(*,*)2.0_prec*cellset(index).MJacobi(i,1,1)
!                !����Ԫ�����άƽ��ֵ��ԭʼ����//�����㲻ΪGauss�㣬����Ҫ�Ȳ�ֵ��Gauss�����ٻ��֣������Ժ����ӵ�Ԫ���������
!                aveCell = Gauss_integral_SPs(varOriCell(i,:,m))/2.0_prec
!                varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(i,:,m),nsp)
!                varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(i,:,m),nsp)
!                !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell1  = cellset(index).nearcells(4)
!                !write(*,*) cellset(index).fpdet_J_F(i,1),cellset(indexNearCell).fpdet_J_F(i,nsp+1)
!                do j = 1,4
!                    if(cellset(indexNearCell1).nearcells(j)==index)then
!                        sideth1 = j
!                    end if
!                end do
!                if((sideth1==3).or.(sideth1==4))then
!                    near_i1 = nsp+1-i
!                else
!                    near_i1 = i
!                end if
!                if(sideth1==2.OR.sideth1==4)then
!                    aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,near_i1,:,m))/2.0_prec
!                else
!                    aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,:,near_i1,m))/2.0_prec
!                end if
!                !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell2  = cellset(index).nearcells(2)
!                do j = 1,4
!                    if(cellset(indexNearCell2).nearcells(j)==index)then
!                        sideth2 = j
!                    end if
!                end do
!                if((sideth2==1).or.(sideth2==2))then
!                    near_i2 = nsp+1-i
!                else
!                    near_i2 = i
!                end if
!                if(sideth2==2.OR.sideth2==4)then
!                    aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,near_i2,:,m))/2.0_prec
!                else
!                    aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,:,near_i2,m))/2.0_prec
!                end if
!    
!                deta_L =aveCell-aveNearCell(4)
!                deta_R = aveNearCell(2)-aveCell           
!                deta_CL = aveCell-varFluxL
!                deta_CR = varFluxR-aveCell
!
!                varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh2)      
!                varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh2)
!            
!                if(varmodL /= deta_CL .or. varmodR /= deta_CR)then
!                    if(detect_type == ByDim )then
!                        cellset(index).Beta_line(i,1)=1
!                        cellset(index).Beta = 1
!                    elseif(detect_type == ByCell)then
!                        cellset(index).Beta_line = 1
!                        cellset(index).Beta = 1
!                        
!                        exit loop1    
!                    end if       
!                end if     
!            
!                !y����
!                dh2 = (2.0_prec*cellset(index).MJacobi(1,i,4))**2
!                !write(*,*)2.0_prec*cellset(index).MJacobi(1,i,4)
!                aveCell = Gauss_integral_SPs(varOriCell(:,i,m))/2.0_prec
!                varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(:,i,m),nsp)
!                varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(:,i,m),nsp)        
!            
!                !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell3 = cellset(index).nearcells(1)
!                do j = 1,4
!                    if(cellset(indexNearCell3).nearcells(j)==index)then
!                        sideth3 = j
!                    end if
!                end do
!                if((sideth3==2).or.(sideth3==1))then
!                    near_i3 = nsp+1-i
!                else
!                    near_i3 = i
!                end if
!                if(sideth3==2.OR.sideth3==4)then
!                    aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,near_i3,:,m))/2.0_prec
!                else
!                    aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,:,near_i3,m))/2.0_prec
!                end if
!                !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell4 = cellset(index).nearcells(3)
!                do j = 1,4
!                    if(cellset(indexNearCell4).nearcells(j)==index)then
!                        sideth4 = j
!                    end if
!                end do
!                if((sideth4==3).or.(sideth4==4))then
!                    near_i4 = nsp+1-i
!                else
!                    near_i4 = i
!                end if
!                if(sideth4==2.OR.sideth4==4)then
!                    aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,near_i4,:,m))/2.0_prec
!                else
!                    aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,:,near_i4,m))/2.0_prec
!                end if
! 
!                deta_L = aveCell-aveNearCell(1)
!                deta_R = aveNearCell(3)-aveCell
!                deta_CL = aveCell-varFluxL
!                deta_CR = varFluxR-aveCell
!            
!                varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh2)           
!                varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh2)
!            
!                if(varmodL /= deta_CL .or. varmodR /= deta_CR)then        
!                    if(detect_type == ByDim )then
!                        cellset(index).Beta_line(i,2)=1
!                    elseif(detect_type == ByCell)then
!                        cellset(index).Beta_line = 1
!                        cellset(index).Beta = 1                   
!                        exit loop1         
!                    end if
!                end if  
!                
!            end do loop2
!        end do loop1
!    end do
!end subroutine indicator_TVB
!
!subroutine indicator_TVBCell
!    !TVB��⣬�����غ����
!    ! ÿһά�����бȽϣ����ֲ�����minmod�����ı���Ϊ���ⵥԪ��ֹͣ����
!    !Refs:
!    !   [1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29. 
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m
!    integer :: sideth1,sideth2,sideth3,sideth4,sideth(4),near_i1,near_i2,near_i3,near_i4,indexNearCell(4),indexNearCell1, indexNearCell2, indexNearCell3, indexNearCell4 
!    integer :: index
!    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh
!    real(prec),dimension(4) :: aveNearCell
!    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
!    real(prec),dimension(nsp,nsp,4) :: varOriCell
!    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
!    real(prec),dimension(nsp) :: varCellL,varCellR
!    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa,aveCellL,aveCellR
!    do index = 1,ncells      
!    
!        !!����Ԫ -- ������ԭʼ����
!        varOriCell(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!
!        !�ڵ�Ԫ
!        do i = 1,4
!            indexNearCell(i)  = cellset(index).nearcells(i)!���ڵ�Ԫ������
!            varOriCellNearCell(i,:,:,:) = cellset(indexNearCell(i)).spvalue_ori(:,:,:)        
!        end do 
!        
!        do j = 1,4  !4�����ڵ�Ԫ
!            do k = 1,4
!                if(index==cellset(indexNearCell(j)).nearcells(k))then
!                    sideth(j) = k   !��¼����Ԫ�Ĳ�������ڵ�Ԫ �ĵڼ����
!                end if
!            end do
!        end do
!        
!        !x����
!        loop1:do m = 1,1!������
!            !write(*,*)m
!            dh = (2.0_prec*cellset(index).MJacobi(3,1,1))!��Ҫ��
!            aveCell = Gauss_double_integral(varOriCell(:,:,:))/4.0_prec !����Ԫ��ֵ
!            do i = 1,nsp
!                varCellL(i) = LaI_nPs(kesi_l,SPs,varOriCell(i,:,m),nsp)
!                varCellR(i) = LaI_nPs(kesi_r,SPs,varOriCell(i,:,m),nsp)
!            end do
!            aveCellL = Gauss_integral_SPs(varCellL(:))/2.0_prec     !��Ԫ��߽��ֵ
!            aveCellR = Gauss_integral_SPs(varCellR(:))/2.0_prec     !��Ԫ�ұ߽��ֵ           
!
!             
!            
!            loop2:do i = 1,nsp!ÿ�зֱ����
!                dh = (2.0_prec*cellset(index).MJacobi(i,1,1))
!                !write(*,*)2.0_prec*cellset(index).MJacobi(i,1,1)
!                !����Ԫ�����άƽ��ֵ��ԭʼ����//�����㲻ΪGauss�㣬����Ҫ�Ȳ�ֵ��Gauss�����ٻ��֣������Ժ����ӵ�Ԫ���������
!                aveCell = Gauss_integral_SPs(varOriCell(i,:,m))/2.0_prec
!                
!                varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(i,:,m),nsp)
!                varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(i,:,m),nsp)
!                !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell1  = cellset(index).nearcells(4)
!                !write(*,*) cellset(index).fpdet_J_F(i,1),cellset(indexNearCell).fpdet_J_F(i,nsp+1)
!                do j = 1,4
!                    if(cellset(indexNearCell1).nearcells(j)==index)then
!                        sideth1 = j
!                    end if
!                end do
!                if((sideth1==3).or.(sideth1==4))then
!                    near_i1 = nsp+1-i
!                else
!                    near_i1 = i
!                end if
!                if(sideth1==2.OR.sideth1==4)then
!                    aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,near_i1,:,m))/2.0_prec
!                else
!                    aveNearCell(4) = Gauss_integral_SPs(varOriCellNearCell(4,:,near_i1,m))/2.0_prec
!                end if
!                !!!----�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell2  = cellset(index).nearcells(2)
!                do j = 1,4
!                    if(cellset(indexNearCell2).nearcells(j)==index)then
!                        sideth2 = j
!                    end if
!                end do
!                if((sideth2==1).or.(sideth2==2))then
!                    near_i2 = nsp+1-i
!                else
!                    near_i2 = i
!                end if
!                if(sideth2==2.OR.sideth2==4)then
!                    aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,near_i2,:,m))/2.0_prec
!                else
!                    aveNearCell(2) = Gauss_integral_SPs(varOriCellNearCell(2,:,near_i2,m))/2.0_prec
!                end if
!    
!                deta_L =aveCell-aveNearCell(4)
!                deta_R = aveNearCell(2)-aveCell           
!                deta_CL = aveCell-varFluxL
!                deta_CR = varFluxR-aveCell
!
!                varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh)      
!                varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh)
!            
!                if(varmodL /= deta_CL .or. varmodR /= deta_CR)then
!                    if(detect_type == ByDim )then
!                        cellset(index).Beta_line(i,1)=1
!                        cellset(index).Beta = 1
!                    elseif(detect_type == ByCell)then
!                        cellset(index).Beta_line = 1
!                        cellset(index).Beta = 1
!                        
!                        exit loop1    
!                    end if       
!                end if     
!            
!                !y����
!                dh = (2.0_prec*cellset(index).MJacobi(1,i,4))
!                !write(*,*)2.0_prec*cellset(index).MJacobi(1,i,4)
!                aveCell = Gauss_integral_SPs(varOriCell(:,i,m))/2.0_prec
!                varFluxL = LaI_nPs(kesi_l,SPs,varOriCell(:,i,m),nsp)
!                varFluxR = LaI_nPs(kesi_r,SPs,varOriCell(:,i,m),nsp)        
!            
!                !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell3 = cellset(index).nearcells(1)
!                do j = 1,4
!                    if(cellset(indexNearCell3).nearcells(j)==index)then
!                        sideth3 = j
!                    end if
!                end do
!                if((sideth3==2).or.(sideth3==1))then
!                    near_i3 = nsp+1-i
!                else
!                    near_i3 = i
!                end if
!                if(sideth3==2.OR.sideth3==4)then
!                    aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,near_i3,:,m))/2.0_prec
!                else
!                    aveNearCell(1) = Gauss_integral_SPs(varOriCellNearCell(1,:,near_i3,m))/2.0_prec
!                end if
!                !�ڵ�Ԫ�뱾��Ԫ�غϱ߽紦 -- ��
!                indexNearCell4 = cellset(index).nearcells(3)
!                do j = 1,4
!                    if(cellset(indexNearCell4).nearcells(j)==index)then
!                        sideth4 = j
!                    end if
!                end do
!                if((sideth4==3).or.(sideth4==4))then
!                    near_i4 = nsp+1-i
!                else
!                    near_i4 = i
!                end if
!                if(sideth4==2.OR.sideth4==4)then
!                    aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,near_i4,:,m))/2.0_prec
!                else
!                    aveNearCell(3) = Gauss_integral_SPs(varOriCellNearCell(3,:,near_i4,m))/2.0_prec
!                end if
! 
!                deta_L = aveCell-aveNearCell(1)
!                deta_R = aveNearCell(3)-aveCell
!                deta_CL = aveCell-varFluxL
!                deta_CR = varFluxR-aveCell
!            
!                varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh)           
!                varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh)
!            
!                if(varmodL /= deta_CL .or. varmodR /= deta_CR)then        
!                    if(detect_type == ByDim )then
!                        cellset(index).Beta_line(i,2)=1
!                    elseif(detect_type == ByCell)then
!                        cellset(index).Beta_line = 1
!                        cellset(index).Beta = 1                   
!                        exit loop1         
!                    end if
!                end if  
!                
!            end do loop2
!        end do loop1
!    end do
!end subroutine indicator_TVBCell
!
!subroutine indicator_MV
!    !Refs:[1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29.
!    !   ʵ���ϣ���ATV����С�ĸĶ���ԭ�����ǲ�ȡAverage Total Variation,�������������ء��ڴ˽�ƽ���ܱ���Ϊ����
!    !   ������ˣ�Ҳ�������⣬�Ǿ��ǹ⻬����������ⵥԪ������ǿ�������ڵ�����£����ܶ��������ⲻ������
!    !   ������֤������1����ȡȫ��������ֵ*theta��Ϊ�ж����ݣ�2����ȡ��ȫ��������ֵ-��Сֵ��*theta��Ϊ�ж����ݡ�
!    use real_precision
!    use type_module
!    use parameter_setting
!    use global_var
!    implicit none
!    integer :: i,j,k,l,m,near_k,ikesi1,ikesi2,ieta1,ieta2,ikesi,ieta
!    integer :: sideth(4),near_kk(4),sidethL,sidethR,indexNearCell,indexNearCell1, indexNearCell2, indexNearCell3, indexNearCell4 
!    integer :: index,indexCellL,indexCellR,Beta,Beta_line_x,Beta_line_y
!    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh2
!    real(prec),dimension(4) :: aveNearCell,TV,ATV,Max_TV,Max_v,Min_v
!    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
!    real(prec),dimension(nsp,nsp,4) :: varOriCell,varOriCellL,varOriCellR
!    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
!    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa,temp,theta,Vtemp_x(4),Vtemp_y(4)
!    theta = 0.1_prec
!    ! �ȼ���TV��the total variation of the solution��.��������֮�䶼��Ҫ���㡣
!    TV = 0.0_prec
!    Max_TV = 0.0_prec
!    Max_v = 0.0_prec
!    Min_v = 0.0_prec
!    !���ȼ���һ����Ԫ�ڲ������TV
!    do i = 1,ncells               
!        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)!!����ȡ����Ԫ -- ������ԭʼ����             
!        do m = 1,1  !������
!            Max_v = MaxVal(varOriCell(:,:,m))
!            Min_v = MinVal(varOriCell(:,:,m))
!            do k = 1,nsp
!                do l = 1,nsp-1
!                    temp = abs(varOriCell(k,l+1,m) - varOriCell(k,l,m)) !x direction
!                    TV(m) = TV(m) + temp
!                    Max_TV(m) = max(Max_TV(m),temp)
!                    temp = abs(varOriCell(l+1,k,m) - varOriCell(l,k,m)) !y
!                    TV(m) = TV(m) + temp
!                    Max_TV(m) = max(Max_TV(m),temp)
!                end do
!            end do      
!        end do
!    end do        
!    !write(*,*)'0000'
!    !Ȼ����ݱ߼���ߵ����ڵ�Ԫ�ڽ��ӱߵ��������
!    do i = 1,nsides        
!        indexCellL = sideset(i).nearcells(1)
!        indexCellR = sideset(i).nearcells(2)
!        if(indexCellL == 0 .OR. indexCellR == 0)cycle
!        varOriCellL(:,:,:) = cellset(indexCellL).spvalue_ori(:,:,:)
!        varOriCellR(:,:,:) = cellset(indexCellR).spvalue_ori(:,:,:)
!        do j = 1,4
!            if(i == cellset(indexCellL).sides(j)) sidethL = j            !���������ڵ�Ԫ�ĵڼ����
!            if(i == cellset(indexCellR).sides(j)) sidethR = j
!        end do
!                
!        do m = 1,1
!            do k = 1,nsp
!                if((sidethL*sidethR==2).or.(sidethL*sidethR==12).or.(sidethL==sidethR))then!������������������������෴
!                    near_k = nsp+1-k
!                else
!                    near_k = k
!                end if
!                if(sidethL==1)then
!                    ikesi1 = 1;     ieta1 = k;
!                elseif(sidethL==2)then
!                    ikesi1 = k;     ieta1 = nsp
!                elseif(sidethL==3)then
!                    ikesi1 = nsp;   ieta1 = k;
!                elseif(sidethL==4)then
!                    ikesi1 = k;     ieta1 = 1
!                end if
!                if(sidethR==1)then
!                    ikesi2 = 1;     ieta2 = near_k
!                elseif(sidethR==2)then
!                    ikesi2 = near_k;ieta2 = nsp
!                elseif(sidethR==3)then
!                    ikesi2 = nsp;   ieta2 = near_k
!                elseif(sidethR==4)then
!                    ikesi2 = near_k;ieta2 = 1
!                end if
!                temp = abs(varOriCellL(ikesi1,ieta1,m)-varOriCellR(ikesi2,ieta2,m))
!                TV(m) = TV(m) + temp
!                Max_TV(m) = max(Max_TV(m),temp)
!            end do
!        end do
!    end do
!    ATV = TV/(ncells*nsp*(nsp-1)*2+(nsides-nbdsides)*nsp)
!    !write(*,*)ATV(1),Max_TV(1)
!    !-------------------------------------------------------
!    !
!    !-------------------------------------------------------
!    !write(*,*)'0001'
!    !!����Ԫ -- ������ԭʼ����
!    do i = 1,ncells
!        Vtemp_x = 0.0_prec
!        Vtemp_y = 0.0_prec
!        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)    
!        !varOriCell(:,:,1) = varOriCell(:,:,1)*varOriCell(:,:,4) !rho*p
!        !varOriCell(:,:,1) = varOriCell(:,:,4)/varOriCell(:,:,1) !p/rho
!        !�ڵ�Ԫ
!        do j = 1,4  !4�����ڵ�Ԫ
!            indexNearCell  = cellset(i).nearcells(j)!���ڵ�Ԫ������
!            !write(*,*)j,indexNearCell,cellset(indexNearCell).nearcells
!            do k = 1,4
!                if(i==cellset(indexNearCell).nearcells(k))then
!                    sideth(j) = k   !��¼����Ԫ�Ĳ�������ڵ�Ԫ �ĵڼ����
!                end if
!            end do
!        
!            varOriCellNearCell(j,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,1)*varOriCellNearCell(i,:,:,4) !rho*p
!            !varOriCellNearCell(i,:,:,1) = varOriCellNearCell(i,:,:,4)/varOriCellNearCell(i,:,:,1) !p/rho
!        end do 
!        !write(*,*)i,sideth
!        !----x direction --------------------------------------------------------------------
!        beta = 0
!        do m = 1,1
!            !write(*,*)'12'
!            do k =1,nsp
!                if(sideth(4)==1)then    !x���򣬿������4�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                    ikesi = 1;          ieta = k;
!                elseif(sideth(4)==2)then
!                    ikesi = k;    ieta = nsp
!                elseif(sideth(4)==3)then
!                    ikesi = nsp;   ieta = nsp+1-k;
!                elseif(sideth(4)==4)then
!                    ikesi = nsp+1-k;     ieta = 1
!                end if  
!                !write(*,*)'13',ikesi,ieta
!                
!                temp =  abs(varOriCell(k,1,m)-varOriCellNearCell(4,ikesi,ieta,m))
!                Vtemp_x(m) = max(Vtemp_x(m),temp)                
!                
!                do l = 2,nsp
!                    temp =  abs(varOriCell(k,l,m)-varOriCell(k,l-1,m))
!                    Vtemp_x(m) = max(Vtemp_x(m),temp)
!                    if(temp>theta*Max_TV(m))then    !MV : theta*Max_TV(m);  ATV : theta*ATV(m)
!                        cellset(i).Beta_line(k,1) = 1
!                        cellset(i).Beta = 1
!                    end if
!                    !write(*,*)'2',temp,theta*Max_TV(m)
!                end do 
!                !x���򣬿������2�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                if(sideth(2)==1)then    
!                    ikesi = 1;          ieta = nsp+1-k;
!                elseif(sideth(2)==2)then
!                    ikesi = nsp+1-k;    ieta = nsp
!                elseif(sideth(2)==3)then
!                    ikesi = nsp;   ieta = k;
!                elseif(sideth(2)==4)then
!                    ikesi = k;     ieta = 1
!                end if  
!                temp =  abs(varOriCell(k,nsp,m)-varOriCellNearCell(2,ikesi,ieta,m))
!                Vtemp_x(m) = max(Vtemp_x(m),temp)
!
!            !--------------------------------------------------------------------------------------------------       
!    
!            !----y direction ----------------------------------------------------------------------------------
!                if(sideth(1)==1)then    !x���򣬿������1�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                    ikesi = 1;          ieta = nsp+1-k;
!                elseif(sideth(1)==2)then
!                    ikesi = nsp+1-k;    ieta = nsp
!                elseif(sideth(1)==3)then
!                    ikesi = nsp;   ieta = k;
!                elseif(sideth(1)==4)then
!                    ikesi = k;     ieta = 1
!                end if  
!                temp =  abs(varOriCell(1,k,m)-varOriCellNearCell(1,ikesi,ieta,m))
!                Vtemp_y(m) = max(Vtemp_y(m),temp)
!        
!                do l = 2,nsp
!                    temp =  abs(varOriCell(l,k,m)-varOriCell(l-1,k,m))
!                    Vtemp_y(m) = max(Vtemp_y(m),temp)
!                end do 
!                !x���򣬿������3�ĵ�һ���㣬��Ҫ���ڵ�Ԫ��Եĵ����ȡ���õ�ĵ�������
!                if(sideth(3)==1)then    
!                    ikesi = 1;          ieta = k;
!                elseif(sideth(3)==2)then
!                    ikesi = k;    ieta = nsp
!                elseif(sideth(3)==3)then
!                    ikesi = nsp;   ieta = nsp+1-k;
!                elseif(sideth(3)==4)then
!                    ikesi = nsp+1-k;     ieta = 1
!                end if  
!                temp =  abs(varOriCell(nsp,k,m)-varOriCellNearCell(3,ikesi,ieta,m))
!                Vtemp_y(m) = max(Vtemp_y(m),temp)
!                
!                if(Vtemp_x(m)>theta*Max_TV(m))then
!                    cellset(i).Beta_line(k,1) = 1
!                    cellset(i).Beta = 1
!                end if
!                if(Vtemp_y(m)>theta*Max_TV(m))then
!                    cellset(i).Beta_line(k,2) = 1
!                    cellset(i).Beta = 1
!                end if
!                
!            end do     
!        end do        
!    end do
!    
!
!end subroutine indicator_MV
    
    
!do RKstage = 1,3          
        !    write(*,*)'RKstage',RKstage   
        !    ������ⵥԪ
        !    call system_clock(int_time1)
        !    call mark_TroCell      
        !    call system_clock(int_time2)
        !    ���±߽����ⵥԪ�����
        !    call update_Boundary_Detect  
        !    call system_clock(int_time3)
        !    �任������ռ�
        !    call phy_to_com
        !    call system_clock(int_time4)
        !    ����߽�ͨ��
        !    call Face_Flux_upw   
        !    call system_clock(int_time5)
        !    RK �ƽ�
        !    call RK_advancing2(RKstage) !---------------------
        !    call system_clock(int_time6)
        !    �غ����->ԭʼ����
        !    call con_to_ori      
        !    !���±߽�
        !    call update_Boundary 
        !    !����ܶȣ�ѹ���Ƿ���ָ�ֵ����������
        !    call correct_density_pressure
        !    ԭʼ����->����ͨ��          
        !    call ori_to_flu
        !    call system_clock(int_time7)
        !    call print_num_data
        !    pause
        !    write(*,*)int_time2 - int_time1
        !    write(*,*)int_time3 - int_time2
        !    write(*,*)int_time4 - int_time3
        !    write(*,*)int_time5 - int_time4
        !    write(*,*)int_time6 - int_time5
        !    write(*,*)int_time7 - int_time6
        !    write(*,*)'-------------'
        !end do
    
    !subroutine time_advancing
!    !ʱ���ƽ�
!    !3rd TVD Runge Kutta method
!    use parameter_setting
!    use global_var
!    use type_module
!    use time_test
!    implicit none
!    integer nt,i,j,k,l,RKstage
!    integer :: int_time1,int_time2,int_time3,int_time4,int_time5,int_time6,int_time7
!    write(*,*) 'Time advance.'
!    !��������Ҫ�õ�����������ڴ�
!    do i = 1,ncells+nbdsides
!        allocate(cellset(i).spvalue_con(nsp,nsp,4),cellset(i).spvalue_fluF(nsp,nsp,4),cellset(i).spvalue_fluG(nsp,nsp,4))
!        allocate(cellset(i).spvalue_con_loc(nsp,nsp,4),cellset(i).spvalue_fluF_loc(nsp,nsp,4),cellset(i).spvalue_fluG_loc(nsp,nsp,4))
!        allocate(cellset(i).spvalue_con_tem(nsp,nsp,4))
!        allocate(cellset(i).fluxF_innerfp(nsp,nsp+1,4),cellset(i).fluxG_innerfp(nsp,nsp+1,4))
!        allocate(cellset(i).spvalue_ori_exa(nsp,nsp,4))
!        allocate(cellset(i).spvalue_ori_old(nsp,nsp,4))
!        allocate(cellset(i).Beta_line(nsp,2))
!        cellset(i).Beta_line = 0
!    end do
!    
!    do i = 1,nsides
!        allocate(sideset(i).fpvalue_upw(nsp,4))
!        sideset(i).fpvalue_upw(:,:) = 0.0_prec !��ʼ��
!    end do
!
!    !���±߽�
!    call update_Boundary_fluent  
!    !ԭʼ����->�غ���� 
!    call ori_to_con         
!    !ԭʼ����->����ͨ��
!    call ori_to_flu    
!
!    !�����غ���
!    open(1110,file = 'error_of_global_con_law.plt')
!    sum_q_initial = 0.0_prec
!    !call exact_solution(0.0_prec)
!    !call solve_errorGCL(0) !error of global conservation law
!    !���ⵥԪ��¼
!    open(1120,file = 'troubleCellRecord.plt')
!    
!    do nt = 1,nstep
!        ! ���ʱ�䲽��
!        !call solve_dt
!        if(nt*dt>T .and. dt_switch == dt_adapt)then
!            dt = T - (nt-1)*dt
!        end if
!        nt_temp = nt
!        if(nt == 1 .OR.nt == 21)then
!            !��ż���һ����Ҫ��ʱ������20������
!            call com_time(nt)
!        end if
!        
!        do RKstage = 1,3                     
!            !������ⵥԪ
!            call mark_TroCell      
!            !!��ӡ���ⵥԪ 
!            !call Trouble_dis
!            !���±߽����ⵥԪ�����
!            call update_Boundary_Detect  
!            !�任������ռ�
!            call phy_to_com
!            !����߽�ͨ��
!            call Face_Flux_upw   
!            !RK �ƽ�
!            call RK_advancing(RKstage) !---------------------
!            !�غ����->ԭʼ����
!            call con_to_ori      
!            !!���±߽�
!            call update_Boundary_fluent   
!            !!����ܶȣ�ѹ���Ƿ���ָ�ֵ����������
!            call correct_density_pressure
!            !ԭʼ����->����ͨ��          
!            call ori_to_flu  
!        end do
!        
!        if(nstep>10)then
!            if(mod(nt,print_step)==0)then
!                write(*,"('stepping'I8'/'I8)") nt,nstep
!                call print_num_data          !ÿnstep/10�� ���һ����ֵ��
!                !call exact_solution(nt*dt)
!                !call exact_solution(T)     !׼ȷ��
!                !call print_exa_data    
!                !stop
!                call now_time
!                !call renewal_program
!                !pause
!                !call solve_errorGCL(nt_temp)
!            end if
!        else
!            write(*,"('stepping'I8'/'I8)") nt,nstep
!        endif
!    end do
!    if(dt_switch == dt_adapt)then
!        dt = dt_temp
!    end if
!    !call solve_errorGCL(nstep)
!    close(1110)
!    close(1120)
!end subroutine time_advancing
    
    !subroutine RK_advancing(stage)
!    !��ϲ���: ���ݵ�Ԫ�Ƿ����ⵥԪ��ȡ��ͬ����
!    use real_precision
!    use global_var
!    use type_module
!    implicit none
!    integer i,stage
!    integer,parameter :: stage1 = 1,stage2 = 2,stage3 = 3
!    real(prec),dimension(:,:,:) :: Lsub(nsp,nsp,4)
!    Lsub = 0.0_prec
!    
!    !ѡ��ʱ�䲽
!    select case(stage)
!    case(stage1)
!        do i = 1,ncells
!            !RK1
!            if(cellset(i).Beta == 0)then
!                call RHS_CPR(i,Lsub)
!            elseif(cellset(i).Beta == 1)then
!                !write(*,*)'22'
!                call RHS_NNW2(i,Lsub)               
!            end if
!            cellset(i).spvalue_con_tem(:,:,:) = cellset(i).spvalue_con(:,:,:)
!            cellset(i).spvalue_con(:,:,:) = cellset(i).spvalue_con_tem(:,:,:) + dt*Lsub(:,:,:)
!            if(isnan(cellset(i).spvalue_con(1,1,1)))then
!                write(*,"('error! NAN occur in',I6,'-th step,',I6,'-th cell')")nt_temp,i       
!                stop
!            end if
!        end do
!
!    case(stage2)
!        do i = 1,ncells
!            !RK2
!            if(cellset(i).Beta == 0)then
!                call RHS_CPR(i,Lsub)
!            elseif(cellset(i).Beta == 1)then
!                call RHS_NNW2(i,Lsub)
!            end if
!            cellset(i).spvalue_con(:,:,:) = 3.0_prec/4.0_prec*cellset(i).spvalue_con_tem(:,:,:) + 1.0_prec/4.0_prec*(cellset(i).spvalue_con(:,:,:)+dt*Lsub(:,:,:))
!        end do        
!    case(stage3)
!        do i = 1,ncells
!            !RK3
!            if(cellset(i).Beta == 0)then
!                call RHS_CPR(i,Lsub)
!            elseif(cellset(i).Beta == 1)then
!                call RHS_NNW2(i,Lsub)
!            end if
!            cellset(i).spvalue_con(:,:,:) = 1.0_prec/3.0_prec*cellset(i).spvalue_con_tem(:,:,:) + 2.0_prec/3.0_prec*(cellset(i).spvalue_con(:,:,:)+dt*Lsub(:,:,:))
!        end do
!    case default
!        write(*,*)'Error! Occur in subroutine RK_advancing.'
!        stop
!    end select
!    
!end subroutine RK_advancing
!subroutine set_BC
!    !�����������ñ߽磬��δ�����������������ã�Ŀǰ����ض�����
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer i,j,k,l
!    integer :: index1 = 1,index2 = 1,index3 = 1,index4 = 1,index0 = 1,index_bd1,index_bd2
!    real(prec) :: dis1,dis2,d_D,d_R,d_U,d_L,p_c1,p_c2,dis_p_c
!    
!    integer :: i_D,i_R,i_U,i_L,O_sideth
!    integer :: indexCell,indexCell2,indexNearCell,indexNearCell2,p1,p2
!    index1 = 1;index2 = 1;index3 = 1;index4 = 1;index0 = 1;
!    n_BCells_D = 0;n_BCells_R = 0;n_BCells_U = 0;n_BCells_L = 0 ;
!    do i = 1,nnodes
!        !���Ǽ�����ı߽磬������ģ�ͣ���ʣ�µ������浥Ԫ
!        d_D = abs(xy_coor(i,2)-yl)!Down
!        d_R = abs(xy_coor(i,1)-xr)!Right
!        d_U = abs(xy_coor(i,2)-yr)!Up
!        d_L = abs(xy_coor(i,1)-xl)!Left
!        if(d_D < 1.0e-10)then
!            n_BCells_D = n_BCells_D + 1 !�²�߽ڵ������
!        end if
!        if(d_R < 1.0e-10)then
!            n_BCells_R = n_BCells_R + 1 !�Ҳ�߽ڵ������
!        end if
!        if(d_U < 1.0e-10)then
!            n_BCells_U = n_BCells_U + 1 !�ϲ�߽ڵ������
!        end if
!        if(d_L < 1.0e-10)then
!            n_BCells_L = n_BCells_L + 1 !���߽ڵ������
!        end if
!    end do
!    
!    !�����߽��Ͻڵ���Ŀ����������
!    n_BCells_D = n_BCells_D - 1 !�ߵ�����
!    n_BCells_R = n_BCells_R - 1 
!    n_BCells_U = n_BCells_U - 1 
!    n_BCells_L = n_BCells_L - 1 
!    !write(*,*) ncells,n_BCells_D, n_BCells_R, n_BCells_U, n_BCells_L
!    nbdsides = n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
!    allocate(BoundCells_index_set(nbdsides))    !��¼��ߵ�Ԫ������
!    !stop
!    do i = 1,ncells    
!        !write(*,*)i
!        i_D = 0
!        do j =1,4
!            d_D = abs(xy_coor(cellset(i).nodes(j),2)-yl)!Down
!            if(d_D < 1.0e-10)then
!                i_D = i_D+1
!            end if
!        end do
!        if(i_D==2)then   !�ױ��ϵĵ�Ԫ
!            O_sideth = 1 !�߽絥Ԫ�ڵױ��ϵĲ�߱��һ��Ϊ1
!            if(cellset(i).nearcells(1)/=0)then  !�Ų��Ų�Ϊ1�ĵ�Ԫ
!                do j = 1,4
!                    if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ1
!                        O_sideth = j 
!                    end if
!                end do
!            end if
!            cellset(i).nearcells(O_sideth) = ncells+index1!�����Ϊncells+index1�Ľṹ�帳Ϊi���ڵ�Ԫ
!            cellset(ncells+index1).nearcells(3) = i
!            BoundCells_index_set(index1)=i      !��¼�߽絥Ԫ��ŵ�����
!            index1 = index1+1
!        end if
!        !write(*,*)'22'
!        i_R = 0
!        do j =1,4
!            d_R = abs(xy_coor(cellset(i).nodes(j),1)-xr)!Right
!            if(d_R < 1.0e-10)then
!                i_R = i_R+1
!            end if
!        end do
!        if(i_R==2)then   !�ױ��ϵĵ�Ԫ
!            O_sideth = 2 !�߽絥Ԫ���ұ��ϵĲ�߱��һ��Ϊ2
!            if(cellset(i).nearcells(2)/=0)then  !�Ų��Ų�Ϊ2�ĵ�Ԫ
!                do j = 1,4
!                    if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ2
!                        O_sideth = j 
!                    end if
!                end do
!            end if
!            cellset(i).nearcells(O_sideth) = ncells+n_BCells_D+index2
!            cellset(ncells+n_BCells_D+index2).nearcells(4) = i
!            BoundCells_index_set(n_BCells_D+index2)=i
!            index2 = index2+1
!        end if
!        !write(*,*)'24'
!        i_U = 0
!        do j =1,4
!            d_U = abs(xy_coor(cellset(i).nodes(j),2)-yr)!Up
!            if(d_U < 1.0e-10)then
!                i_U = i_U+1
!            end if
!        end do
!        if(i_U==2)then   !�ױ��ϵĵ�Ԫ
!            O_sideth = 3 !�߽絥Ԫ���ϱ��ϵĲ�߱��һ��Ϊ3
!            if(cellset(i).nearcells(3)/=0)then  !�Ų��Ų�Ϊ3�ĵ�Ԫ
!                do j = 1,4
!                    if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ3
!                        O_sideth = j 
!                    end if
!                end do
!            end if
!            cellset(i).nearcells(O_sideth) = ncells+n_BCells_D+n_BCells_R+index3
!            cellset(ncells+n_BCells_D+n_BCells_R+index3).nearcells(1) = i
!            BoundCells_index_set(n_BCells_D+n_BCells_R+index3)=i
!            index3 = index3+1
!        end if
!        i_L = 0
!        do j =1,4
!            d_L = abs(xy_coor(cellset(i).nodes(j),1)-xl)!Up
!            if(d_L < 1.0e-10)then
!                i_L = i_L+1
!            end if
!        end do
!        !write(*,*)'25'
!        if(i_L==2)then   !�ױ��ϵĵ�Ԫ
!            O_sideth = 4 !�߽絥Ԫ������ϵĲ�߱��һ��Ϊ4
!            if(cellset(i).nearcells(4)/=0)then  !�Ų��Ų�Ϊ4�ĵ�Ԫ
!                do j = 1,4
!                    if(cellset(i).nearcells(j)==0)then  !���迼�ǹսǴ��ĵ�Ԫ����Ϊ�սǴ���Ϊ4
!                        O_sideth = j 
!                    end if
!                end do
!            end if
!            cellset(i).nearcells(O_sideth) = ncells+n_BCells_D+n_BCells_R+n_BCells_U+index4
!            !write(*,*)'27',ncells,ncells+nbdsides,ncells+n_BCells_D+n_BCells_R+n_BCells_U+index4
!            cellset(ncells+n_BCells_D+n_BCells_R+n_BCells_U+index4).nearcells(2) = i
!            !write(*,*)'28'
!            BoundCells_index_set(n_BCells_D+n_BCells_R+n_BCells_U+index4)=i
!            index4 = index4+1
!        end if
!       
!    end do
!    index0 = index0-1       !�߽����
!    !write(*,*)index1,index2,index3,index4,index0
!    !write(*,*)BoundCells_index_set
!    !do i = ncells+1,ncells+nbdsides
!    !    write(*,*)cellset(i).nearcells(:)    
!    !end do
!  
!    !do i = 1,nbdsides
!    !    write(*,*) BoundCells_index_set(i),cellset(BoundCells_index_set(i)).nearcells(:)
!    !end do
!    !���������趨��ӵ���Щ��չ��Ԫ��ָ��λ��
!    if(case_comp==equEntropy_case.or.case_comp==test_case.OR. case_comp == CompositeVortexShock_case)then      !�����в�ȡѭ���߽�         
!        Loop1:do i = 1,n_BCells_D                 !�±߽�
!        !write(*,*)'i',i
!        !����BoundCells_index_set(i)���ң���С��Χ
!            Loop2:do j = 1,4                
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                !write(*,*)'j',j,indexCell,indexNearCell
!                if(indexNearCell > ncells .and.indexNearCell < ncells+n_BCells_D+1 )then!��λ���߽��ϵĲ��,Down
!                    p1 = j          !��ߵ�Ԫ�����
!                    p2 = p1+1        !��ߵ�Ԫ�յ���
!                    if(p1 == 4)p2=1
!                    !write(*,*)p1,p2
!                    !if(j==1)p2=2
!                    p_c1 = (xy_coor(cellset(indexCell).nodes(p1),1)+xy_coor(cellset(indexCell).nodes(p2),1))*0.5_prec
!                    !write(*,*) p_c1 
!                    Loop3:do k = n_BCells_D+n_BCells_R+1,n_BCells_D+n_BCells_R+n_BCells_U
!                        !write(*,*)'k',k
!                        Loop4:do l = 1,4
!                            !write(*,*)'l',l
!                            indexCell2 = BoundCells_index_set(k)
!                            indexNearCell2 = cellset(indexCell2).nearcells(l)
!                            !write(*,*)indexCell2,l,indexNearCell2 ,ncells+n_BCells_D+n_BCells_R,ncells+n_BCells_D+n_BCells_R+n_BCells_U
!                            if((indexNearCell2 > ncells+n_BCells_D+n_BCells_R) .and.(indexNearCell2 < ncells+n_BCells_D+n_BCells_R+n_BCells_U+1) )then!��λ���߽��ϵĲ��,Up
!                                p1 = l          !��ߵ�Ԫ�����
!                                p2 = p1+1!��ߵ�Ԫ�յ���
!                                if(p1 == 4)p2=1
!                                !write(*,*)'000',p1,p2
!                                p_c2 = (xy_coor(cellset(indexCell2).nodes(p1),1)+xy_coor(cellset(indexCell2).nodes(p2),1))*0.5_prec
!                                dis_p_c = abs(p_c1 - p_c2)
!                                !write(*,*) dis_p_c
!                                if(dis_p_c<1.0e-10)then!���е�������غϣ����Ӧ������ֻ�����ڹ�ص������ϣ�����������бһ���Ƕ��޷�����
!                                    cellset(indexNearCell).index = indexCell2
!                                    cellset(indexNearCell2).index = indexCell
!                                    !write(*,*) '---------',i,l,cellset(indexNearCell).index ,cellset(indexNearCell2).index 
!                                   exit Loop2                                   
!                                end if
!                            end if
!                        end do Loop4
!                    end do Loop3
!                end if
!            end do Loop2
!        end do Loop1
!        Loop5:do i = n_BCells_D+1,n_BCells_D+n_BCells_R                 !�ұ߽�
!        !����BoundCells_index_set(i)���ң���С��Χ
!            Loop6:do j = 1,4
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                if((indexNearCell > ncells+n_BCells_D) .and.indexNearCell < (ncells+n_BCells_D+n_BCells_R+1) )then!��λ���߽��ϵĲ��,Right
!                    p1 = j          !��ߵ�Ԫ�����
!                    p2 = p1+1        !��ߵ�Ԫ�յ���
!                    if(p1 == 4)p2=1
!                    p_c1 = (xy_coor(cellset(indexCell).nodes(p1),2)+xy_coor(cellset(indexCell).nodes(p2),2))*0.5_prec
!                    Loop7:do k = n_BCells_D+n_BCells_R+n_BCells_U+1,n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
!                        Loop8:do l = 1,4
!                            indexCell2 = BoundCells_index_set(k)
!                            indexNearCell2 = cellset(indexCell2).nearcells(l)
!                            if(indexNearCell2 > ncells+n_BCells_D+n_BCells_R+n_BCells_U .and.indexNearCell2 < ncells+n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L+1 )then!��λ���߽��ϵĲ��,Up
!                                p1 = l          !��ߵ�Ԫ�����
!                                p2 = p1+1        !��ߵ�Ԫ�յ���
!                                if(p1 == 4)p2=1
!                                p_c2 = (xy_coor(cellset(indexCell2).nodes(p1),2)+xy_coor(cellset(indexCell2).nodes(p2),2))*0.5_prec
!                                dis_p_c = abs(p_c1 - p_c2)
!                                if(dis_p_c<1.0e-10)then!���е�������غϣ����Ӧ������ֻ�����ڹ�ص������ϣ�����������бһ���Ƕ��޷�����
!                                    cellset(indexNearCell).index = indexCell2
!                                    cellset(indexNearCell2).index = indexCell
!                                    !write(*,*) i,cellset(indexNearCell).index ,cellset(indexNearCell2).index 
!                                    exit Loop6
!                                end if
!                            end if
!                        end do Loop8
!                    end do Loop7
!                end if
!            end do Loop6
!        end do Loop5
!    elseif(case_comp == Riemann2D_case .OR. case_comp == DoubleMach_case .OR. case_comp == VortexShock_case )then   !��άRiemann���⣬����һ�㣻˫������⣬����һ�㣬�ڸ��±߽����������и���
!        do i = 1,nbdsides                
!            do j = 1,4                
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                !write(*,*)'j',j,indexCell,indexNearCell
!                if(indexNearCell > ncells )then!��λ���߽��ϵĲ��,Down
!                    cellset(indexNearCell).index = indexCell                   
!                end if
!            end do
!        end do 
!    elseif((case_comp==SodShockTube_case.OR.case_comp==LaxShockTube_case .OR. case_comp==ShuOsher_case).AND.dire_shock	== 0)then
!        !�������ڱ߽磬���ҽ�֧�߽�
!        Loop11:do i = 1,n_BCells_D                 !�±߽�
!        !����BoundCells_index_set(i)���ң���С��Χ
!            Loop12:do j = 1,4                
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                if(indexNearCell > ncells .and.indexNearCell < ncells+n_BCells_D+1 )then!��λ���߽��ϵĲ��,Down
!                    p1 = j          !��ߵ�Ԫ�����
!                    p2 = p1+1        !��ߵ�Ԫ�յ���
!                    if(p1 == 4)p2=1                 
!                    p_c1 = (xy_coor(cellset(indexCell).nodes(p1),1)+xy_coor(cellset(indexCell).nodes(p2),1))*0.5_prec
!                    Loop13:do k = n_BCells_D+n_BCells_R+1,n_BCells_D+n_BCells_R+n_BCells_U
!                        Loop14:do l = 1,4
!                            indexCell2 = BoundCells_index_set(k)
!                            indexNearCell2 = cellset(indexCell2).nearcells(l)
!                            if((indexNearCell2 > ncells+n_BCells_D+n_BCells_R) .and.(indexNearCell2 < ncells+n_BCells_D+n_BCells_R+n_BCells_U+1) )then!��λ���߽��ϵĲ��,Up
!                                p1 = l          !��ߵ�Ԫ�����
!                                p2 = p1+1!��ߵ�Ԫ�յ���
!                                if(p1 == 4)p2=1
!                                p_c2 = (xy_coor(cellset(indexCell2).nodes(p1),1)+xy_coor(cellset(indexCell2).nodes(p2),1))*0.5_prec
!                                dis_p_c = abs(p_c1 - p_c2)                                
!                                if(dis_p_c<1.0e-10)then!���е�������غϣ����Ӧ������ֻ�����ڹ�ص������ϣ�����������бһ���Ƕ��޷�����
!                                    cellset(indexNearCell).index = indexCell2
!                                    cellset(indexNearCell2).index = indexCell                                    
!                                   exit Loop12                                   
!                                end if
!                            end if
!                        end do Loop14
!                    end do Loop13
!                end if
!            end do Loop12
!        end do Loop11
!        
!        Loop15:do i = n_BCells_D+1,n_BCells_D+n_BCells_R                 !�ұ߽�
!        !����BoundCells_index_set(i)���ң���С��Χ
!            Loop16:do j = 1,4
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                if((indexNearCell > ncells+n_BCells_D) .and.indexNearCell < (ncells+n_BCells_D+n_BCells_R+1) )then!��λ���߽��ϵĲ��,Right
!                    cellset(indexNearCell).index = indexCell 
!                end if
!            end do Loop16
!        end do Loop15
!        Loop17:do k = n_BCells_D+n_BCells_R+n_BCells_U+1,n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
!            Loop18:do l = 1,4
!                indexCell2 = BoundCells_index_set(k)
!                indexNearCell2 = cellset(indexCell2).nearcells(l)
!                if((indexNearCell2 > ncells+n_BCells_D+n_BCells_R+n_BCells_U) .and.(indexNearCell2 < ncells+n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L+1) )then!��λ���߽��ϵĲ��,Up
!                    cellset(indexNearCell2).index = indexCell2
!                end if
!            end do Loop18
!        end do Loop17
!    elseif((case_comp==SodShockTube_case.OR.case_comp==LaxShockTube_case .OR. case_comp==ShuOsher_case).AND.dire_shock	== 1)then
!        !���������������������֤��y���򴫲��������õ��Ĳ���
!        Loop21:do i = 1,n_BCells_D                 !�±߽�
!            Loop22:do j = 1,4                
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                if(indexNearCell > ncells .and.indexNearCell < ncells+n_BCells_D+1 )then!��λ���߽��ϵĲ��,Down
!                   cellset(indexNearCell).index = indexCell             
!                end if
!            end do Loop22
!        end do Loop21
!        Loop23:do k = n_BCells_D+n_BCells_R+1,n_BCells_D+n_BCells_R+n_BCells_U
!            Loop24:do l = 1,4
!                indexCell2 = BoundCells_index_set(k)
!                indexNearCell2 = cellset(indexCell2).nearcells(l)
!                if((indexNearCell2 > ncells+n_BCells_D+n_BCells_R) .and.(indexNearCell2 < ncells+n_BCells_D+n_BCells_R+n_BCells_U+1) )then!��λ���߽��ϵĲ��,Up
!                    cellset(indexNearCell2).index = indexCell2
!                end if
!            end do Loop24
!        end do Loop23
!        Loop25:do i = n_BCells_D+1,n_BCells_D+n_BCells_R                 !�ұ߽�
!
!            Loop26:do j = 1,4
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)
!                if((indexNearCell > ncells+n_BCells_D) .and.indexNearCell < (ncells+n_BCells_D+n_BCells_R+1) )then!��λ���߽��ϵĲ��,Right
!                    p1 = j          !��ߵ�Ԫ�����
!                    p2 = p1+1        !��ߵ�Ԫ�յ���
!                    if(p1 == 4)p2=1
!                    p_c1 = (xy_coor(cellset(indexCell).nodes(p1),2)+xy_coor(cellset(indexCell).nodes(p2),2))*0.5_prec
!                    Loop27:do k = n_BCells_D+n_BCells_R+n_BCells_U+1,n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L
!                        Loop28:do l = 1,4
!                            indexCell2 = BoundCells_index_set(k)
!                            indexNearCell2 = cellset(indexCell2).nearcells(l)
!                            if(indexNearCell2 > ncells+n_BCells_D+n_BCells_R+n_BCells_U .and.indexNearCell2 < ncells+n_BCells_D+n_BCells_R+n_BCells_U+n_BCells_L+1 )then!��λ���߽��ϵĲ��,Up
!                                p1 = l          !��ߵ�Ԫ�����
!                                p2 = p1+1        !��ߵ�Ԫ�յ���
!                                if(p1 == 4)p2=1
!                                p_c2 = (xy_coor(cellset(indexCell2).nodes(p1),2)+xy_coor(cellset(indexCell2).nodes(p2),2))*0.5_prec
!                                dis_p_c = abs(p_c1 - p_c2)
!                                if(dis_p_c<1.0e-10)then!���е�������غϣ����Ӧ������ֻ�����ڹ�ص������ϣ�����������бһ���Ƕ��޷�����
!                                    cellset(indexNearCell).index = indexCell2
!                                    cellset(indexNearCell2).index = indexCell
!                                    exit Loop26
!                                end if
!                            end if
!                        end do Loop28
!                    end do Loop27
!                end if
!            end do Loop26
!        end do Loop25  
!    end if  
!    !write(*,*)nbdsides
!    !do i = 1,nbdsides!��֤ÿ��������ĵ�Ԫ����ȷ�����ڵ�Ԫ
!    !    write(*,*) BoundCells_index_set(i)!,cellset(BoundCells_index_set(i)).nearcells
!    !    write(*,*)  cellset(cellset(BoundCells_index_set(i)).nearcells(1)).index,&
!    !                cellset(cellset(BoundCells_index_set(i)).nearcells(2)).index,&
!    !                cellset(cellset(BoundCells_index_set(i)).nearcells(3)).index,&
!    !                cellset(cellset(BoundCells_index_set(i)).nearcells(4)).index
!    !    write(*,*)
!    !end do
!    !stop
!   
!           
!   !do i  = ncells+1,ncells+20
!   ! write(*,*) i,cellset(i).index
!   !end do
!   !stop
!end subroutine set_BC
!
!subroutine pre_BoundaryCell
!    !   ������߽絥ԪԤ�������õ�Ԫ���꣬������ʵ���꣬Jacobi�Ȳ���ʱ��仯��ֵ
!    !   ��Ԫ��״�����ڵ�Ԫ�ر߽��߶Գơ�����ѭ���߽�������case֮��������¸�ֵ������صı��������ص��ĵ�Ԫ��״����������˵����
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer :: i,j,k,l,indexCell,indexNearCell,index_temp,p1,p2,index,th1,nearcells(4),pCounter,sideth1,vertexth(4)
!    real(prec),dimension(:,:) :: quad_vertex(4,2)
!    real(prec),dimension(:,:,:) :: weight_sps(nsp,nsp,4)
!    real(prec) :: coor_C(2),M_Jacobi(4),M_direct(4)
!    pCounter = nnodes+1
!    do i = 1,nbdsides    
!    !write(*,*)i
!        do j = 1,4         
!            
!            ! �߽������ڵ�Ԫ
!            indexCell = BoundCells_index_set(i)
!            
!            ! �߽������ڵ�Ԫ ���ڵ�Ԫ
!            indexNearCell = cellset(indexCell).nearcells(j)
!
!            ! �ױ� i�����Ƿ��ǵױߣ�indexNearCell�����Ƿ������ⵥԪ
!            if(i>0 .and. i < n_BCells_D +1 .and. indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then!�±߽�
!                ! sideth1��ʾ�߽絥Ԫ�Ǳ߽絥Ԫ�ĵڼ���ߣ� vertexth(1)��ʾ���ⵥԪ�����1.���ݶԳƹ�ϵ�Ƶ�
!                ! 4 _____________3
!                ! |               |
!                ! |               |
!                ! |1_____________2|
!                !------------------sym
!                ! 1 _____________2
!                ! |               |
!                ! |               |
!                ! |4_____________3|
!                sideth1 = j
!                vertexth(1) = sideth1-1
!                if(vertexth(1)==0) vertexth(1) = 4
!                do k = 2,4
!                    vertexth(k) = vertexth(k-1)-1
!                    if(vertexth(k)==0) vertexth(k) = 4
!                end do
!                do k = 1,4
!                    cellset(indexNearCell).nodes(k) = cellset(indexCell).nodes(vertexth(k)) 
!                end do
!                
!                xy_coor(pCounter,1) = xy_coor(cellset(indexNearCell).nodes(1),1)
!                xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(4),2)-xy_coor(cellset(indexNearCell).nodes(1),2)
!                cellset(indexNearCell).nodes(1) = pCounter
!                pCounter = pCounter + 1
!                xy_coor(pCounter,1) = xy_coor(cellset(indexNearCell).nodes(2),1)
!                xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(3),2)-xy_coor(cellset(indexNearCell).nodes(2),2)
!                cellset(indexNearCell).nodes(2) = pCounter
!                pCounter = pCounter + 1
!
!                !do k =1,4
!                !    write(*,*)(xy_coor(cellset(indexNearCell).nodes(k),l),l=1,2)
!                !end do 
!            elseif(i>n_BCells_D  .and. i < n_BCells_D + n_BCells_R +1 .and.indexNearCell>ncells+n_BCells_D .and. indexNearCell < ncells + n_BCells_D +n_BCells_R +1)then!�ұ߽�
!                sideth1 = j
!                vertexth(1) = sideth1
!                if(vertexth(1)==0) vertexth(1) = 4
!                do k = 2,4
!                    vertexth(k) = vertexth(k-1)-1
!                    if(vertexth(k)==0) vertexth(k) = 4
!                end do
!                do k = 1,4
!                    cellset(indexNearCell).nodes(k) = cellset(indexCell).nodes(vertexth(k)) 
!                end do
!                
!                xy_coor(pCounter,1) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(1),1)-xy_coor(cellset(indexNearCell).nodes(2),1)
!                xy_coor(pCounter,2) = xy_coor(cellset(indexNearCell).nodes(2),2)
!                cellset(indexNearCell).nodes(2) = pCounter
!                pCounter = pCounter + 1
!                xy_coor(pCounter,1) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(4),1)-xy_coor(cellset(indexNearCell).nodes(3),1)
!                xy_coor(pCounter,2) = xy_coor(cellset(indexNearCell).nodes(3),2)
!                cellset(indexNearCell).nodes(3) = pCounter
!                pCounter = pCounter + 1
!
!                !do k =1,4
!                !    write(*,*)(xy_coor(cellset(indexNearCell).nodes(k),l),l=1,2)
!                !end do 
!                     
!            elseif(i>n_BCells_D+ n_BCells_R  .and. i < n_BCells_D + n_BCells_R + n_BCells_U +1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+1)then!�ϱ߽�
!                sideth1 = j
!                vertexth(1) = sideth1+1
!                if(vertexth(1)==0) vertexth(1) = 4
!                do k = 2,4
!                    vertexth(k) = vertexth(k-1)-1
!                    if(vertexth(k)==0) vertexth(k) = 4
!                end do
!                do k = 1,4
!                    cellset(indexNearCell).nodes(k) = cellset(indexCell).nodes(vertexth(k)) 
!                end do
!                
!                xy_coor(pCounter,1) = xy_coor(cellset(indexNearCell).nodes(4),1)
!                xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(1),2)-xy_coor(cellset(indexNearCell).nodes(4),2)
!                cellset(indexNearCell).nodes(4) = pCounter
!                pCounter = pCounter + 1
!                xy_coor(pCounter,1) = xy_coor(cellset(indexNearCell).nodes(3),1)
!                xy_coor(pCounter,2) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(2),2)-xy_coor(cellset(indexNearCell).nodes(3),2)
!                cellset(indexNearCell).nodes(3) = pCounter
!                pCounter = pCounter + 1
!
!                !do k =1,4
!                !    write(*,*)(xy_coor(cellset(indexNearCell).nodes(k),l),l=1,2)
!                !end do 
!            elseif(i>n_BCells_D+ n_BCells_R + n_BCells_U .and. i < n_BCells_D + n_BCells_R + n_BCells_U + n_BCells_L+1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R+ n_BCells_U .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+ n_BCells_L+1)then!��߽�
!                sideth1 = j
!                vertexth(1) = sideth1-2
!                if(vertexth(1)==0) vertexth(1) = 4
!                do k = 2,4
!                    vertexth(k) = vertexth(k-1)-1
!                    if(vertexth(k)==0) vertexth(k) = 4
!                end do
!                do k = 1,4
!                    cellset(indexNearCell).nodes(k) = cellset(indexCell).nodes(vertexth(k)) 
!                end do
!                
!                xy_coor(pCounter,1) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(2),1)-xy_coor(cellset(indexNearCell).nodes(1),1)
!                xy_coor(pCounter,2) = xy_coor(cellset(indexNearCell).nodes(1),2)
!                cellset(indexNearCell).nodes(1) = pCounter
!                pCounter = pCounter + 1
!                xy_coor(pCounter,1) = 2.0_prec*xy_coor(cellset(indexNearCell).nodes(3),1)-xy_coor(cellset(indexNearCell).nodes(4),1)
!                xy_coor(pCounter,2) = xy_coor(cellset(indexNearCell).nodes(4),2)
!                cellset(indexNearCell).nodes(4) = pCounter
!                pCounter = pCounter + 1
!
!                !do k =1,4
!                !    write(*,*)(xy_coor(cellset(indexNearCell).nodes(k),l),l=1,2)
!                !end do         
!            end if         
!        end do     
!    end do  
!    
!    !ѭ���������߽絥ԪJacobi�ͽ��ȫ������
!    do i = ncells+1,ncells+nbdsides   
!        !--------------------------------------------------------------------------------------------
!        !ȡ���ı��ε�Ԫ��������ֵ
!        do j = 1,4
!            do k = 1,2
!                quad_vertex(j,k) = xy_coor(cellset(i).nodes(j),k)
!            end do
!        end do
!
!        !���ڲ��������
!        allocate(cellset(i).sp_coor(nsp,nsp,2))!�����ڴ棬�ڲ����ȫ������    (nsp,nsp)     
!        do j = 1,nsp
!            do k =1,nsp
!                call quad_C2Phy(quad_vertex,SPs_local(j,k,:),cellset(i).sp_coor(j,k,:))
!            end do
!        end do        
!        !--------------------------------------------------------------------------------------------        
!        !����Jacobi ����  [xkesi,ykesi;xeta,yeta]
!        allocate(cellset(i).MJacobi(nsp,nsp,4),cellset(i).Mdirect(nsp,nsp,4),cellset(i).det_J(nsp,nsp))!�����ڴ��MJacobi
!        do j = 1,nsp
!            do k = 1,nsp
!                call solve_Jacobi(quad_vertex,SPs_local(j,k,:),cellset(i).MJacobi(j,k,:),cellset(i).Mdirect(j,k,:),cellset(i).det_J(j,k))                   
!            end do
!        end do
!        
!        !-------------------------------------------------------------------------------------------
!        !-------------------------------------------------------------------------------------------
!        !ͨ����Jacobi ����  [xkesi,ykesi;xeta,yeta]
!        allocate(cellset(i).fpMdirect_F(nsp,nsp+1,4),cellset(i).fpMdirect_G(nsp,nsp+1,4),cellset(i).fpdet_J_F(nsp,nsp+1),cellset(i).fpdet_J_G(nsp,nsp+1))!�����ڴ��MJacobi
!        do j = 1,nsp
!            coor_C(2) = SPs(j)
!            do k = 1,nsp+1
!                coor_C(1) = FPs(k)               
!                call solve_Jacobi(quad_vertex,coor_C(:),M_Jacobi,cellset(i).fpMdirect_F(j,k,:),cellset(i).fpdet_J_F(j,k))  
!                !write(*,"(F15.6)") M_Jacobi(1)
!            end do 
!        end do
!        do j = 1,nsp
!            coor_C(1) = SPs(j)
!            do k = 1,nsp+1
!                coor_C(2) = FPs(k)               
!                call solve_Jacobi(quad_vertex,coor_C(:),M_Jacobi,cellset(i).fpMdirect_G(j,k,:),cellset(i).fpdet_J_G(j,k))         
!                !write(*,"(F15.6)") M_Jacobi(3)
!            end do 
!        end do  
!    end do
!end subroutine pre_BoundaryCell
!
!subroutine update_Boundary
!    !ÿһʱ�䲽��ÿһʱ��㶼��Ҫ���±߽�
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none
!    integer :: i,j,k,l,m,indexCell,indexNearCell,index_temp,p1,p2,index,th1,nearcells(4)
!    integer :: ikesi,ieta
!    real(prec) :: x,y,ss,tt,ruvp(4),coor(2)
!    if(case_comp == DoubleMach_case)then
!        !˫��շ��䣬�±߽��ȡ�޴�͸���ȱ߽粿����Ҫÿʱ�䲽���¸�ֵ
!        do i = 1,nbdsides    
!            indexCell = BoundCells_index_set(i)
!            do j = 1,4                 
!                indexNearCell = cellset(indexCell).nearcells(j)          
!                if(i>0 .and. i < n_BCells_D +1 .and.indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then
!                !�±߽�
!                    !write(*,*)j,indexCell,indexNearCell
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp+1-k,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,k,:)                   
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(k,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!                            end do
!                        end if                                                  
!                    end do
!                    cellset(indexNearCell).spvalue_ori(:,:,3) = -cellset(indexNearCell).spvalue_ori(:,:,3)
!                    do k = 1,nsp                       
!                        do l = 1,nsp
!                            !if(j==1)then      
!                            !    ikesi  = nsp+1-l;     ieta = k;
!                            !elseif(j==2)then  
!                            !    ikesi = k;     ieta = nsp+1-l;
!                            !elseif(j==3)then  
!                            !    ikesi = nsp+1-l;   ieta = nsp+1-k;
!                            !elseif(j==4)then  
!                            !    ikesi = nsp+1-k;     ieta = nsp+1-l; 
!                            !end if                        
!                            !
!                            !cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(ikesi,ieta,:) 
!                            !cellset(indexNearCell).spvalue_ori(l,k,3) = -cellset(indexNearCell).spvalue_ori(l,k,3)
!                            
!                            x = cellset(indexNearCell).sp_coor(l,k,1)
!                            y = cellset(indexNearCell).sp_coor(l,k,2)
!                            !write(*,*) x,y
!                            if(x < 1.0_prec/6.0_prec)then
!                                !write(*,*) x,y
!                                ruvp(1)  = 8.0_prec
!                                ruvp(2)  = 7.145_prec
!                                ruvp(3)  = -4.125_prec
!                                ruvp(4)  = 116.5_prec    
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = ruvp
!                            end if
!                        end do        
!                    end do                    
!                                             
!                elseif(i>n_BCells_D+ n_BCells_R  .and. i < n_BCells_D + n_BCells_R + n_BCells_U +1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+1)then
!                !�ϱ߽�               
!                    
!                    do k = 1, nsp
!                        do l = 1,nsp
!                            x = cellset(indexNearCell).sp_coor(k,l,1)
!                            y = cellset(indexNearCell).sp_coor(k,l,2)
!                            ss = 10.0_prec
!                            tt = nt_temp*dt
!                            if(y - sqrt(3.0_prec)*(x - 1.0_prec/6.0_prec) > -ss*tt*2.0_prec  )then
!                                ruvp(1) = 8.0_prec
!                                ruvp(2)  = 7.145_prec
!                                ruvp(3)  = -4.125_prec
!                                ruvp(4)  = 116.5_prec                         
!                            else
!                                ruvp(1)  = 1.4_prec
!                                ruvp(2)  = 0.0_prec
!                                ruvp(3)  = 0.0_prec
!                                ruvp(4)  = 1.0_prec
!                            endif
!                            cellset(indexNearCell).spvalue_ori(k,l,:) = ruvp
!                        end do 
!                    end do
!                                        
!                elseif(indexNearCell > ncells)then
!                !�� �ұ߽�             
!                    cellset(indexNearCell).spvalue_ori = cellset(indexCell).spvalue_ori 
!                end if         
!            end do     
!        end do
!    elseif(case_comp == Riemann2D_case)then
!        !----xy��Ϊ��֧�߽�------------------------
!        !do i = ncells+1,ncells+nbdsides    
!        !    indexCell = i
!        !    indexNearCell = cellset(indexCell).index
!        !    do j = 1,4
!        !        if(cellset(indexNearCell).nearcells(j)==indexCell)then
!        !            exit
!        !        end if
!        !    end do
!        !    !write(*,*)j,indexNearCell,indexCell
!        !    if(i>ncells+0 .and. i < ncells+n_BCells_D +1)then
!        !    !�±߽�
!        !        do k = 1,nsp                       
!        !            if(j==1)then      
!        !                cellset(indexCell).spvalue_ori(k,:,:) = cellset(indexNearCell).spvalue_ori(nsp+1-k,:,:)                           
!        !            elseif(j==2)then  
!        !                cellset(indexCell).spvalue_ori(k,:,:) = cellset(indexNearCell).spvalue_ori(:,k,:)                   
!        !            elseif(j==3)then  
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(k,l,:) = cellset(indexNearCell).spvalue_ori(k,nsp+1-l,:) 
!        !                end do
!        !            elseif(j==4)then                              
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(k,l,:) = cellset(indexNearCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!        !                end do
!        !            end if                                                
!        !        end do
!        !    elseif(i>ncells+n_BCells_D  .and. i <ncells+ n_BCells_D + n_BCells_R +1 )then
!        !    !�ұ߽�
!        !
!        !        do k = 1,nsp                       
!        !            if(j==1)then      
!        !                cellset(indexCell).spvalue_ori(:,k,:) = cellset(indexNearCell).spvalue_ori(k,:,:)                           
!        !            elseif(j==2)then  
!        !                cellset(indexCell).spvalue_ori(:,k,:) = cellset(indexNearCell).spvalue_ori(:,nsp+1-k,:)                 
!        !            elseif(j==3)then  
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(l,k,:) = cellset(indexNearCell).spvalue_ori(nsp+1-k,nsp+1-l,:) 
!        !                end do
!        !            elseif(j==4)then                              
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(l,k,:) = cellset(indexNearCell).spvalue_ori(nsp+1-l,k,:)
!        !                end do
!        !            end if                                                
!        !        end do
!        !
!        !    elseif(i>ncells+n_BCells_D+ n_BCells_R  .and. i < ncells+n_BCells_D + n_BCells_R + n_BCells_U +1 )then
!        !    !�ϱ߽�
!        !
!        !        do k = 1,nsp                       
!        !            if(j==1)then      
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(k,l,:) = cellset(indexNearCell).spvalue_ori(k,nsp+1-l,:) 
!        !                end do                           
!        !            elseif(j==2)then  
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(k,l,:) = cellset(indexNearCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!        !                end do            
!        !            elseif(j==3)then  
!        !                cellset(indexCell).spvalue_ori(k,:,:) = cellset(indexNearCell).spvalue_ori(nsp+1-k,:,:)                    
!        !            elseif(j==4)then   
!        !                cellset(indexCell).spvalue_ori(k,:,:) = cellset(indexNearCell).spvalue_ori(:,k,:)                        
!        !            end if                                                
!        !        end do
!        !
!        !    elseif(i>ncells+n_BCells_D+ n_BCells_R + n_BCells_U .and. i < ncells+n_BCells_D + n_BCells_R + n_BCells_U + n_BCells_L+1 )then!��߽�
!        !    !��߽�
!        !
!        !        do k = 1,nsp                       
!        !            if(j==1)then      
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(l,k,:) = cellset(indexNearCell).spvalue_ori(nsp+1-k,nsp+1-l,:) 
!        !                end do                    
!        !            elseif(j==2)then  
!        !                do l = 1,nsp
!        !                    cellset(indexCell).spvalue_ori(l,k,:) = cellset(indexNearCell).spvalue_ori(nsp+1-l,k,:)
!        !                end do           
!        !            elseif(j==3)then  
!        !                cellset(indexCell).spvalue_ori(:,k,:) = cellset(indexNearCell).spvalue_ori(k,:,:)                     
!        !            elseif(j==4)then    
!        !                cellset(indexCell).spvalue_ori(:,k,:) = cellset(indexNearCell).spvalue_ori(:,nsp+1-k,:)                       
!        !            end if                                                
!        !        end do   
!        !
!        !    end if         
!        !end do
!        !do i = 1,nbdsides    
!        !    indexCell = BoundCells_index_set(i)
!        !    do j = 1,4                         
!        !        indexNearCell = cellset(indexCell).nearcells(j)
!        !        !write(*,*)j,indexCell,indexNearCell
!        !        if(i>0 .and. i < n_BCells_D +1 .and.indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then
!        !        !�±߽�
!        !            do m=1,4
!        !                do k = 1,nsp                       
!        !                    if(j==1)then      
!        !                        cellset(indexNearCell).spvalue_ori(k,:,m) = cellset(indexCell).spvalue_ori(nsp+1-k,:,m)                           
!        !                    elseif(j==2)then  
!        !                        cellset(indexNearCell).spvalue_ori(k,:,m) = cellset(indexCell).spvalue_ori(:,k,m)                   
!        !                    elseif(j==3)then  
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(k,l,m) = cellset(indexCell).spvalue_ori(k,nsp+1-l,m) 
!        !                        end do
!        !                    elseif(j==4)then                              
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(k,l,m) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,m)
!        !                        end do
!        !                    end if                                                
!        !                end do
!        !            end do
!        !        elseif(i>n_BCells_D  .and. i < n_BCells_D + n_BCells_R +1 .and.indexNearCell>ncells+n_BCells_D .and. indexNearCell < ncells + n_BCells_D +n_BCells_R +1)then
!        !        !�ұ߽�
!        !            do m = 1,4
!        !                do k = 1,nsp                       
!        !                    if(j==1)then      
!        !                        cellset(indexNearCell).spvalue_ori(:,k,m) = cellset(indexCell).spvalue_ori(k,:,m)                           
!        !                    elseif(j==2)then  
!        !                        cellset(indexNearCell).spvalue_ori(:,k,m) = cellset(indexCell).spvalue_ori(:,nsp+1-k,m)                 
!        !                    elseif(j==3)then  
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(l,k,m) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp+1-l,m) 
!        !                        end do
!        !                    elseif(j==4)then                              
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(l,k,m) = cellset(indexCell).spvalue_ori(nsp+1-l,k,m)
!        !                        end do
!        !                    end if                                                
!        !                end do
!        !            end do
!        !        elseif(i>n_BCells_D+ n_BCells_R  .and. i < n_BCells_D + n_BCells_R + n_BCells_U +1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+1)then
!        !        !�ϱ߽�
!        !            do m = 1,4
!        !                do k = 1,nsp                       
!        !                    if(j==1)then      
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(k,l,m) = cellset(indexCell).spvalue_ori(k,nsp+1-l,m) 
!        !                        end do                           
!        !                    elseif(j==2)then  
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(k,l,m) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,m)
!        !                        end do            
!        !                    elseif(j==3)then  
!        !                        cellset(indexNearCell).spvalue_ori(k,:,m) = cellset(indexCell).spvalue_ori(nsp+1-k,:,m)                    
!        !                    elseif(j==4)then   
!        !                        cellset(indexNearCell).spvalue_ori(k,:,m) = cellset(indexCell).spvalue_ori(:,k,m)                        
!        !                    end if                                                
!        !                end do
!        !            end do
!        !        elseif(i>n_BCells_D+ n_BCells_R + n_BCells_U .and. i < n_BCells_D + n_BCells_R + n_BCells_U + n_BCells_L+1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R+ n_BCells_U .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+ n_BCells_L+1)then!��߽�
!        !        !��߽�
!        !            do m =1,4
!        !                do k = 1,nsp                       
!        !                    if(j==1)then      
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(l,k,m) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp+1-l,m) 
!        !                        end do                    
!        !                    elseif(j==2)then  
!        !                        do l = 1,nsp
!        !                            cellset(indexNearCell).spvalue_ori(l,k,m) = cellset(indexCell).spvalue_ori(nsp+1-l,k,m)
!        !                        end do           
!        !                    elseif(j==3)then  
!        !                        cellset(indexNearCell).spvalue_ori(:,k,m) = cellset(indexCell).spvalue_ori(k,:,m)                     
!        !                    elseif(j==4)then    
!        !                        cellset(indexNearCell).spvalue_ori(:,k,m) = cellset(indexCell).spvalue_ori(:,nsp+1-k,m)                       
!        !                    end if                                                
!        !                end do   
!        !            end do
!        !        end if         
!        !    end do     
!        !end do
!        do i = 1,nbdsides    
!            indexCell = BoundCells_index_set(i)
!            do j = 1,4                         
!                indexNearCell = cellset(indexCell).nearcells(j)
!                !write(*,*)j,indexCell,indexNearCell
!                if(i>0 .and. i < n_BCells_D +1 .and.indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then
!                !�±߽�
!                              
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(1,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,nsp,:)                   
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,1,:)
!                            end do
!                        end if                                                
!                    end do
!                elseif(i>n_BCells_D  .and. i < n_BCells_D + n_BCells_R +1 .and.indexNearCell>ncells+n_BCells_D .and. indexNearCell < ncells + n_BCells_D +n_BCells_R +1)then
!                !�ұ߽�
!                         
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(1,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(:,nsp,:)                 
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-l,1,:)
!                            end do
!                        end if                                                
!                    end do
!                     
!                elseif(i>n_BCells_D+ n_BCells_R  .and. i < n_BCells_D + n_BCells_R + n_BCells_U +1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+1)then
!                !�ϱ߽�
!                     do k = 1,nsp                       
!                        if(j==1)then      
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(1,nsp+1-l,:) 
!                            end do                           
!                        elseif(j==2)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp,:)
!                            end do            
!                        elseif(j==3)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp,:,:)                    
!                        elseif(j==4)then   
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,1,:)                        
!                        end if                                                
!                    end do
!                elseif(i>n_BCells_D+ n_BCells_R + n_BCells_U .and. i < n_BCells_D + n_BCells_R + n_BCells_U + n_BCells_L+1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R+ n_BCells_U .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+ n_BCells_L+1)then!��߽�
!                !��߽�
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(1,nsp+1-l,:) 
!                            end do                    
!                        elseif(j==2)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp,:)
!                            end do           
!                        elseif(j==3)then  
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(nsp,:,:)                     
!                        elseif(j==4)then    
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(:,1,:)                       
!                        end if                                                
!                    end do          
!                end if         
!            end do     
!        end do        
!                        
!    elseif(case_comp == SodShockTube_case .or. case_comp == LaxShockTube_case .or. case_comp == ShuOsher_case)then
!        !����ѭ���߽����������ң�ShockTube_case��֧�߽������������ڱ߽��ϵ�ֵ����ShuOsher_caseΪ���֧������߽��������߽����أ���ĿǰҲ��ֱ�ӵ��ڱ߽�ֵ��
!        !��Щ��������ֱ����߽絥Ԫָ���Ӧλ�ã�֮�󣬲���Ҫ������ֵ
!        do i = 1,nbdsides
!            do j = 1,4
!                indexCell = BoundCells_index_set(i)
!                indexNearCell = cellset(indexCell).nearcells(j)                                              
!                if(indexNearCell>ncells)then
!                    if(case_comp == ShuOsher_case .and. indexNearCell>ncells+n_BCells_D .and. indexNearCell < ncells + n_BCells_D +n_BCells_R +1)then
!                        !�Ҳ�߱��ֳ�ʼֵ
!                    elseif(indexNearCell>ncells+n_BCells_D+n_BCells_R+ n_BCells_U .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+ n_BCells_L+1)then!��߽�
!                        !���߱��ֳ�ʼֵ
!                    else
!                        !���࣬����߽絥Ԫindex��ָ����ͬ��ѭ���߽�
!                        nearcells = cellset(indexNearCell).nearcells
!                        cellset(indexNearCell) = cellset(cellset(indexNearCell).index)    !��ֵ
!                        cellset(indexNearCell).nearcells = nearcells
!                    end if
!                end if            
!            end do        
!        end do
!    elseif(case_comp == equEntropy_case)then
!        do i = ncells+1,ncells+nbdsides
!            !index = cellset(i).index
!            !cellset(i) = cellset(cellset(i).index)
!            !cellset(i).index = index
!            cellset(i).spvalue_ori = cellset(cellset(i).index).spvalue_ori   
!        end do
!    elseif(case_comp == VortexShock_case)then
!        !----------------------------
!        do i = 1,nbdsides    
!            indexCell = BoundCells_index_set(i)
!            do j = 1,4                         
!                indexNearCell = cellset(indexCell).nearcells(j)
!                !write(*,*)j,indexCell,indexNearCell
!                if(i>0 .and. i < n_BCells_D +1 .and.indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then
!                !�±߽�
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp+1-k,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,k,:)                   
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(k,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!                            end do
!                        end if                                                  
!                    end do
!                    cellset(indexNearCell).spvalue_ori(:,:,3) = -cellset(indexNearCell).spvalue_ori(:,:,3)
!                elseif(i>n_BCells_D  .and. i < n_BCells_D + n_BCells_R +1 .and.indexNearCell>ncells+n_BCells_D .and. indexNearCell < ncells + n_BCells_D +n_BCells_R +1)then
!                !�ұ߽�
!                   
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(k,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(:,nsp+1-k,:)                 
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-l,k,:)
!                            end do
!                        end if                                                
!                    end do
!            
!                elseif(i>n_BCells_D+ n_BCells_R  .and. i < n_BCells_D + n_BCells_R + n_BCells_U +1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+1)then
!                !�ϱ߽�
!    
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(k,nsp+1-l,:) 
!                            end do                           
!                        elseif(j==2)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!                            end do            
!                        elseif(j==3)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp+1-k,:,:)                    
!                        elseif(j==4)then   
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,k,:)                        
!                        end if                                                
!                    end do
!                    cellset(indexNearCell).spvalue_ori(:,:,3) = -cellset(indexNearCell).spvalue_ori(:,:,3)
!                elseif(i>n_BCells_D+ n_BCells_R + n_BCells_U .and. i < n_BCells_D + n_BCells_R + n_BCells_U + n_BCells_L+1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R+ n_BCells_U .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+ n_BCells_L+1)then!��߽�
!                !��߽�
!                    !do k = 1,nsp                       
!                    !    if(j==1)then      
!                    !        do l = 1,nsp
!                    !            cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp+1-l,:) 
!                    !        end do                    
!                    !    elseif(j==2)then  
!                    !        do l = 1,nsp
!                    !            cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-l,k,:)
!                    !        end do           
!                    !    elseif(j==3)then  
!                    !        cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(k,:,:)                     
!                    !    elseif(j==4)then    
!                    !        cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(:,nsp+1-k,:)                       
!                    !    end if                                                
!                    !end do   
!                    ruvp(1) = 1.0_prec                    
!                    ruvp(2) = Ms*sqrt(gamma)
!                    ruvp(3) = 0.0_prec
!                    ruvp(4) = 1.0_prec
!                    do k = 1, nsp
!                        do l = 1,nsp
!                            !coor(1) = cellset(indexNearCell).sp_coor(l,k,1)
!                            !coor(2) = cellset(indexNearCell).sp_coor(l,k,2)
!                            !call VortexShock_init(coor,ruvp)
!                            !
!                            cellset(indexNearCell).spvalue_ori(k,l,:) = ruvp
!                        end do 
!                    end do
!                    
!                end if         
!            end do     
!        end do
!    elseif(case_comp == CompositeVortexShock_case )then
!        do i = 1,nbdsides    
!            indexCell = BoundCells_index_set(i)
!            do j = 1,4                         
!                indexNearCell = cellset(indexCell).nearcells(j)
!                !write(*,*)j,indexCell,indexNearCell
!                if(i>0 .and. i < n_BCells_D +1 .and.indexNearCell>ncells .and. indexNearCell < ncells + n_BCells_D +1)then
!                !�±߽�
!                    !����
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp+1-k,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,k,:)                   
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(k,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!                            end do
!                        end if                                                  
!                    end do
!                    cellset(indexNearCell).spvalue_ori(:,:,3) = -cellset(indexNearCell).spvalue_ori(:,:,3)
!                    !ѭ��
!                    !cellset(indexNearCell).spvalue_ori = cellset(cellset(indexNearCell).index).spvalue_ori 
!                    !����
!                    !do k = 1,nsp                       
!                    !    if(j==1)then      
!                    !        cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(1,:,:)                           
!                    !    elseif(j==2)then  
!                    !        cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,nsp,:)                   
!                    !    elseif(j==3)then  
!                    !        do l = 1,nsp
!                    !            cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp,nsp+1-l,:) 
!                    !        end do
!                    !    elseif(j==4)then                              
!                    !        do l = 1,nsp
!                    !            cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,1,:)
!                    !        end do
!                    !    end if                                                
!                    !end do
!                elseif(i>n_BCells_D  .and. i < n_BCells_D + n_BCells_R +1 .and.indexNearCell>ncells+n_BCells_D .and. indexNearCell < ncells + n_BCells_D +n_BCells_R +1)then
!                !�ұ߽�
!                   
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(k,:,:)                           
!                        elseif(j==2)then  
!                            cellset(indexNearCell).spvalue_ori(:,k,:) = cellset(indexCell).spvalue_ori(:,nsp+1-k,:)                 
!                        elseif(j==3)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-k,nsp+1-l,:) 
!                            end do
!                        elseif(j==4)then                              
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(l,k,:) = cellset(indexCell).spvalue_ori(nsp+1-l,k,:)
!                            end do
!                        end if                                                
!                    end do
!            
!                elseif(i>n_BCells_D+ n_BCells_R  .and. i < n_BCells_D + n_BCells_R + n_BCells_U +1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+1)then
!                !�ϱ߽�
!                    !����
!                    do k = 1,nsp                       
!                        if(j==1)then      
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(k,nsp+1-l,:) 
!                            end do                           
!                        elseif(j==2)then  
!                            do l = 1,nsp
!                                cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp+1-k,:)
!                            end do            
!                        elseif(j==3)then  
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp+1-k,:,:)                    
!                        elseif(j==4)then   
!                            cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,k,:)                        
!                        end if                                                
!                    end do
!                    cellset(indexNearCell).spvalue_ori(:,:,3) = -cellset(indexNearCell).spvalue_ori(:,:,3)
!                    !ѭ��
!                    !cellset(indexNearCell).spvalue_ori = cellset(cellset(indexNearCell).index).spvalue_ori 
!                    !����
!                    !do k = 1,nsp                       
!                    !    if(j==1)then      
!                    !        do l = 1,nsp
!                    !            cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(1,nsp+1-l,:) 
!                    !        end do                           
!                    !    elseif(j==2)then  
!                    !        do l = 1,nsp
!                    !            cellset(indexNearCell).spvalue_ori(k,l,:) = cellset(indexCell).spvalue_ori(nsp+1-l,nsp,:)
!                    !        end do            
!                    !    elseif(j==3)then  
!                    !        cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(nsp,:,:)                    
!                    !    elseif(j==4)then   
!                    !        cellset(indexNearCell).spvalue_ori(k,:,:) = cellset(indexCell).spvalue_ori(:,1,:)                        
!                    !    end if                                                
!                    !end do
!                elseif(i>n_BCells_D+ n_BCells_R + n_BCells_U .and. i < n_BCells_D + n_BCells_R + n_BCells_U + n_BCells_L+1 .and.indexNearCell>ncells+n_BCells_D+n_BCells_R+ n_BCells_U .and. indexNearCell < ncells + n_BCells_D +n_BCells_R+n_BCells_U+ n_BCells_L+1)then!��߽�
!                !��߽�
!  
!                    ruvp(1) = 1.0_prec                    
!                    ruvp(2) = Ms*sqrt(gamma)
!                    ruvp(3) = 0.0_prec
!                    ruvp(4) = 1.0_prec
!                    do k = 1, nsp
!                        do l = 1,nsp
!                            cellset(indexNearCell).spvalue_ori(k,l,:) = ruvp
!                        end do 
!                    end do
!                    
!                end if         
!            end do     
!        end do
!    end if
!    
!    !do i = 1,nbdsides!��֤ÿ��������ĵ�Ԫ����ȷ�����ڵ�Ԫ
!    !    write(*,*) BoundCells_index_set(i)!,cellset(BoundCells_index_set(i)).nearcells
!    !    write(*,*)  cellset(cellset(BoundCells_index_set(i)).nearcells(1)).index,&
!    !                cellset(cellset(BoundCells_index_set(i)).nearcells(2)).index,&
!    !                cellset(cellset(BoundCells_index_set(i)).nearcells(3)).index,&
!    !                cellset(cellset(BoundCells_index_set(i)).nearcells(4)).index
!    !    write(*,*)
!    !end do
!    !stop
!end subroutine update_Boundary
!    
!    
!subroutine get_oriL(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,oriL)
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none 
!    integer ::index,nearCellIndex,LC_sideth,RC_sideth,l,m,k,RC_l
!    real(prec),dimension(:,:) :: u(3,4),conL(4),ll(4,4),rr(4,4),chara_con(3,4),ui_all(8,4)
!    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4),varConL(nsp,nsp,4),varConR(nsp,nsp,4)
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),oriL(4),chara_oriL(4),chara_conL(4),uA(4)
!    real(prec) :: realDis1,realDis2,detJ_fp,dP1P2
!    integer :: beta_dim
!    !ȡ�����ҵ�Ԫ�ϵ�ԭʼ����   
!    varOriL(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!    varOriR(:,:,:) = cellset(nearCellIndex).spvalue_ori(:,:,:)
!    !ȡ�����ҵ�Ԫ�ϵ��غ����
!    varConL = cellset(index).spvalue_con  !��������غ����
!    varConR = cellset(nearCellIndex).spvalue_con
!    !write(*,*)'0211'
!    if(LC_sideth == 2 .or.LC_sideth == 4)then
!        beta_dim = cellset(index).Beta_line(l,1)
!    elseif(LC_sideth == 1 .or.LC_sideth == 3)then
!        beta_dim = cellset(index).Beta_line(l,2)
!    end if
!    !!�������������ǿ������ⵥԪ�ı߽�ֵ����2�ײ�ֵ�õ�
!    !if(RC_sideth == 2 .or.RC_sideth == 4)then
!    !    beta_dim = max(beta_dim ,cellset(nearCellIndex).Beta_line(RC_l,1))
!    !elseif(RC_sideth == 1 .or.RC_sideth == 3)then
!    !    beta_dim = max(beta_dim ,cellset(nearCellIndex).Beta_line(RC_l,2))
!    !end if
!    !write(*,*)'0212'
!    !if(cellset(index).Beta == 0)then
!    if(beta_dim == 0)then
!    
!        if(LC_sideth == 1)then                                      
!            call LaI_nPs_arr(kesi_l,SPs,varOriL(:,l,:),nsp,oriL)
!        else if(LC_sideth == 2)then
!            call LaI_nPs_arr(kesi_r,SPs,varOriL(l,:,:),nsp,oriL)
!        else if(LC_sideth == 3)then
!            call LaI_nPs_arr(kesi_r,SPs,varOriL(:,l,:),nsp,oriL)
!        else if(LC_sideth == 4)then
!            call LaI_nPs_arr(kesi_l,SPs,varOriL(l,:,:),nsp,oriL)
!        end if
!       
!    else if(beta_dim == 1)then 
!        if(var_type == ori_type .OR. var_type == character_type)then
!            !!----ԭʼ����---------------------!!----��������---------------------
!            if(LC_sideth == 1)then    
!                u(2,:) = varOriL(1,l,:)
!                u(3,:) = varOriL(2,l,:)
!            else if(LC_sideth == 2)then
!                u(2,:) = varOriL(l,nsp,:)
!                u(3,:) = varOriL(l,nsp-1,:)
!            else if(LC_sideth == 3)then
!                u(2,:) = varOriL(nsp,l,:)
!                u(3,:) = varOriL(nsp-1,l,:)
!            else if(LC_sideth == 4)then
!                u(2,:) = varOriL(l,1,:)
!                u(3,:) = varOriL(l,2,:)
!            end if     
!
!            call get_nearCellInfo(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,varOriL,varOriR,u,uA,ui_all,dP1P2)
!
!            call FaceFluxC2NNW2(index,nearCellIndex,LC_sideth,l,u,uA,ui_all,dP1P2,oriL(:),innerDisOri)
!            !!----ԭʼ����---------------------!!----��������---------------------  
!        end if
!        !write(*,*)'24'
!!write(*,*)'0213'
!        if(LC_sideth==1)then    !�ӵ�Ԫ�Ȱѽ���ԭʼ������������֮��絥Ԫ��������һ��      
!            cellset(index).fluxG_innerfp(l,2,:) = innerDisOri
!        elseif(LC_sideth==2)then
!            cellset(index).fluxF_innerfp(l,nsp,:) = innerDisOri
!        elseif(LC_sideth==3)then
!            cellset(index).fluxG_innerfp(l,nsp,:) = innerDisOri
!        elseif(LC_sideth==4)then
!            cellset(index).fluxF_innerfp(l,2,:) = innerDisOri
!        end if
!        !if(index == 20 .and. nearCellIndex==21)then
!        !    write(*,*)u(2,4),uA(4)
!        !    write(*,*)oriL(4),innerDisOri(4)
!        !    !stop
!        !end if
!        
!        !write(*,*)'0214'
!    end if
!end subroutine get_oriL
!
!subroutine get_oriR(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,oriR)
!    !�ǽṹ�ı������в����򣬿��Ǽ������Ӧ�Ĳ�߱�Ź�ϵ������ͨ����������
!    !����һ�ַ�������һ���м�������ڵ�Ԫ���С���ת��ʹ֮��Ӧ������1-3��2-4��Ӧ������
!    use global_var
!    use parameter_setting
!    use type_module
!    implicit none 
!    integer :: index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,m,j,k
!    real(prec),dimension(:,:) :: u(3,4),ll(4,4),rr(4,4),chara_u(3,4),chara_con(3,4),ui_all(8,4)
!    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4),varConL(nsp,nsp,4),varConR(nsp,nsp,4)
!    real(prec),external :: LaI_nPs
!    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),conR(4),oriR(4),chara_conR(4),uA(4)
!    real(prec) :: realDis,detJ_fp,dP1P2
!    integer :: beta_dim
!    !ȡ�����ҵ�Ԫ�ϵ�ԭʼ����   
!    varOriL(:,:,:) = cellset(index).spvalue_ori(:,:,:)
!    varOriR(:,:,:) = cellset(nearCellIndex).spvalue_ori(:,:,:)
!
!    !ȡ�����ҵ�Ԫ�ϵ��غ���� 
!    varConL(:,:,:) = cellset(index).spvalue_con(:,:,:)
!    varConR(:,:,:) = cellset(nearCellIndex).spvalue_con(:,:,:)
!    !����Ҳ��㴦��Jacobi����ʽֵ���任������ռ�
!    !if(RC_sideth == 1)then
!    !    detJ_fp = cellset(nearCellIndex).fpdet_J_G(RC_l,1)
!    !elseif(RC_sideth == 2)then
!    !    detJ_fp = cellset(nearCellIndex).fpdet_J_F(RC_l,nsp+1)
!    !else if(RC_sideth == 3)then
!    !    detJ_fp = cellset(nearCellIndex).fpdet_J_G(RC_l,nsp+1)
!    !else if(RC_sideth == 4)then
!    !    detJ_fp = cellset(nearCellIndex).fpdet_J_F(RC_l,1)
!    !end if
!    if(RC_sideth == 2 .or.RC_sideth == 4)then
!        beta_dim = cellset(nearCellIndex).Beta_line(RC_l,1)
!    elseif(RC_sideth == 1 .or. RC_sideth == 3)then
!        beta_dim = cellset(nearCellIndex).Beta_line(RC_l,2)
!    end if
!    !!!�������������ǿ������ⵥԪ�ı߽�ֵ����2�ײ�ֵ�õ�
!    !if(LC_sideth == 2 .or.LC_sideth == 4)then
!    !    beta_dim = max(beta_dim ,cellset(index).Beta_line(l,1))
!    !elseif(LC_sideth == 1 .or.LC_sideth == 3)then
!    !    beta_dim = max(beta_dim ,cellset(index).Beta_line(l,2))
!    !end if
!    !if(cellset(nearCellIndex).Beta == 0)then
!    if(beta_dim == 0)then
!
!        if(RC_sideth == 1)then
!            call LaI_nPs_arr(kesi_l,SPs,varOriR(:,RC_l,:),nsp,oriR)
!        elseif(RC_sideth == 2)then
!            call LaI_nPs_arr(kesi_r,SPs,varOriR(RC_l,:,:),nsp,oriR)
!        else if(RC_sideth == 3)then
!            call LaI_nPs_arr(kesi_r,SPs,varOriR(:,RC_l,:),nsp,oriR)
!        else if(RC_sideth == 4)then
!            call LaI_nPs_arr(kesi_l,SPs,varOriR(RC_l,:,:),nsp,oriR)
!        end if      
!    else if(beta_dim == 1)then
!        if(var_type==ori_type .OR. var_type==character_type)then
!            !----ԭʼ����!��������------------------------------------- 
!            if(RC_sideth == 1)then
!                u(2,:) = varOriR(1,RC_l,:)
!                u(3,:) = varOriR(2,RC_l,:)   
!            else if(RC_sideth == 2)then
!                u(2,:) = varOriR(RC_l,nsp,:)
!                u(3,:) = varOriR(RC_l,nsp-1,:)               
!            else if(RC_sideth == 3)then
!                u(2,:) = varOriR(nsp,RC_l,:)
!                u(3,:) = varOriR(nsp-1,RC_l,:)
!            else if(RC_sideth == 4)then
!                u(2,:) = varOriR(RC_l,1,:)
!                u(3,:) = varOriR(RC_l,2,:)
!            end if           
!            !write(*,*)'RR'
!            call get_nearCellInfo(nearCellIndex,index,RC_sideth,LC_sideth,RC_l,l,varOriR,varOriL,u(:,:),uA,ui_all,dP1P2)
!            call FaceFluxC2NNW2(nearCellIndex,index,RC_sideth,RC_l,u,uA,ui_all,dP1P2,oriR(:),innerDisOri)!��ֵ���̿���ѡ��������������ԭʼ����
!            !----ԭʼ����-------------------------------------
!        end if        
!
!        if(RC_sideth==1)then    !�ӵ�Ԫ�Ȱѽ���ԭʼ������������֮��絥Ԫ��������һ��      
!            cellset(nearCellIndex).fluxG_innerfp(RC_l,2,:) = innerDisOri
!        elseif(RC_sideth==2)then
!            cellset(nearCellIndex).fluxF_innerfp(RC_l,nsp,:) = innerDisOri
!        elseif(RC_sideth==3)then
!            cellset(nearCellIndex).fluxG_innerfp(RC_l,nsp,:) = innerDisOri
!        elseif(RC_sideth==4)then
!            cellset(nearCellIndex).fluxF_innerfp(RC_l,2,:) = innerDisOri
!        end if
!        
!     
!    end if
!end subroutine get_oriR
    
!subroutine clearMechineZero(Array)
!    use real_precision
!    use parameter_setting
!    implicit none
!    real(prec) ::Array(nsp,nsp)
!    integer :: k,l,m
!    
!    do k = 1,nsp
!        do l = 1,nsp
!            if(abs(Array(k,l))<1.0e-15)then
!                Array(k,l) = 0.0_prec
!            end if                
!        end do
!    end do
!    
!end subroutine clearMechineZero