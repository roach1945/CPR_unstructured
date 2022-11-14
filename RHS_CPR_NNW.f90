module RHS_module
     
    !-----------------------------------------------------------------------------
    !
    !   ������RHS �õ��Ĺ��ñ���
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use real_precision
    use parameter_setting,only:nsp
    implicit none
    
    real(prec),dimension(4) :: fu_sub_kesi,gu_sub_eta
    real(prec),dimension(:,:),  allocatable :: detJ
    real(prec),dimension(:,:,:),allocatable :: FG_upw
    real(prec),dimension(:,:,:),allocatable :: O1,O2,Q1,F1,G1 
    real(prec),dimension(:,:,:),allocatable :: fp_fu_kesi,fp_gu_eta
    
end module RHS_module
    
subroutine RHS_CPR_NNW(index,L_sub_global)

    use global_var
    use parameter_setting
    use type_module
    use time_test
    use RHS_module
    implicit none
    
    integer k,l,m,index,i,j,ii
    real(prec),dimension(4)   :: f_l,f_r,g_l,g_r
    real(prec),dimension(:,:) :: Q_fp(nsp+1,4),O_fp(nsp+1,4),F_fp(nsp+1,4),G_fp(nsp+1,4)
    real(prec),dimension(:,:,:) :: L_sub_global(nsp,nsp,4),F_sub_kesi(nsp,nsp,4),G_sub_eta(nsp,nsp,4)
        
    allocate(FG_upw(4,nsp,4),O1(nsp,nsp,4),O2(nsp,nsp,4),Q1(nsp,nsp,4),F1(nsp,nsp,4),G1(nsp,nsp,4),detJ(nsp,nsp))
         
    !-----------------------------------------------------------------------------
    !
    !   �����������ӵ�Ԫͨ��
    !
    !-----------------------------------------------------------------------------
     
    call get_subcellFlux(index)

    !-----------------------------------------------------------------------------
    !
    !   ���������� F_kesi
    !
    !-----------------------------------------------------------------------------
    
    do k = 1,nsp
        if(cellset(index).Beta_line(k,1)==0)then
   
            if(fluxDer_switch == fluxDer_SPs)then
                
                ! nsp�����㴦ͨ��ֵ��ֵ,ͨ������������ֵ����
                
                call LaI_nPs_arr(kesi_l,SPs,F1(k,:,:),nsp,f_l(:))
                call LaI_nPs_arr(kesi_r,SPs,F1(k,:,:),nsp,f_r(:))

                do l = 1,nsp              
                    call LaI_nPs_deri_arr(SPs(l),SPs,F1(k,:,:),nsp,fu_sub_kesi)
                    F_sub_kesi(k,l,:) = fu_sub_kesi(:)+(FG_upw(4,k,:)-f_l(:))*gl_coe(l) + (FG_upw(2,k,:)-f_r(:))*gr_coe(l)                 
                end do
                
            elseif(fluxDer_switch == fluxDer_FPs)then
                
                ! nsp+1��ͨ���㴦ͨ��ֵ��ֵ,ͨ��������ͨ����ֵ����
                
                F_fp(:,:) = cellset(index).fluxF_innerfp(k,:,:) 
                f_l(:) = F_fp(1,:)
                f_r(:) = F_fp(nsp+1,:)
                do l = 1,nsp              
                    call LaI_nPs_deri_arr(SPs(l),FPs,F_fp,nsp+1,fu_sub_kesi)
                    F_sub_kesi(k,l,:) = fu_sub_kesi(:)+(FG_upw(4,k,:)-f_l(:))*gl_coe(l) + (FG_upw(2,k,:)-f_r(:))*gr_coe(l)                 
                end do       
            end if
            
        elseif(cellset(index).Beta_line(k,1)==1)then
            
            if(method_adv_Der == Operator_2Ps)then
                
                ! �����������  1st order
                do l =1,nsp
                    F_sub_kesi(k,l,:) = (cellset(index).fluxF_innerfp(k,l+1,:)-cellset(index).fluxF_innerfp(k,l,:))/(dis_sp_fp(2*l-1)+dis_sp_fp(2*l))    
                end do
            elseif(method_adv_Der == Operator_3Ps)then
                
                ! 2nd order �������Ȩ
                call solve_adv_der(cellset(index).fluxF_innerfp(k,:,:),cellset(index).spvalue_fluF_loc(k,:,:),F_sub_kesi(k,:,:))!�������
            end if       
            
        end if
    end do
    
    !-----------------------------------------------------------------------------
    !
    !   ���������� G_eta
    !
    !-----------------------------------------------------------------------------
    
    do k = 1,nsp
        if(cellset(index).Beta_line(k,2)==0)then        
            
            if(fluxDer_switch == fluxDer_SPs)then
                
                ! nsp�����㴦ͨ��ֵ��ֵ,ͨ������������ֵ����
                call LaI_nPs_arr(kesi_l,SPs,G1(:,k,:),nsp,g_l(:))
                call LaI_nPs_arr(kesi_r,SPs,G1(:,k,:),nsp,g_r(:))
                do l = 1,nsp              
                    call LaI_nPs_deri_arr(SPs(l),SPs,G1(:,k,:),nsp,gu_sub_eta)
                    G_sub_eta(l,k,:) = gu_sub_eta(:) +(FG_upw(1,k,:)-g_l(:))*gl_coe(l) + (FG_upw(3,k,:)-g_r(:))*gr_coe(l)                 
                end do
                
            elseif(fluxDer_switch == fluxDer_FPs)then
                
                ! nsp+1��ͨ���㴦ͨ��ֵ��ֵ,ͨ��������ͨ����ֵ����
                G_fp(:,:) = cellset(index).fluxG_innerfp(k,:,:) 
                g_l(:) = G_fp(1,:)
                g_r(:) = G_fp(nsp+1,:)
                do l = 1,nsp              
                    call LaI_nPs_deri_arr(SPs(l),FPs,G_fp,nsp+1,gu_sub_eta)
                    G_sub_eta(l,k,:) = gu_sub_eta(:) +(FG_upw(1,k,:)-g_l(:))*gl_coe(l) + (FG_upw(3,k,:)-g_r(:))*gr_coe(l)                
                end do
            end if
                  
        elseif(cellset(index).Beta_line(k,2)==1)then
            
            if(method_adv_Der == Operator_2Ps)then
                
                ! �����������   
                do l =1,nsp         
                    G_sub_eta(l,k,:) = (cellset(index).fluxG_innerfp(k,l+1,:)-cellset(index).fluxG_innerfp(k,l,:))/(dis_sp_fp(2*l-1)+dis_sp_fp(2*l))            
                end do
            elseif(method_adv_Der == Operator_3Ps)then
                
                ! 2nd order �������Ȩ
                call solve_adv_der(cellset(index).fluxG_innerfp(k,:,:),cellset(index).spvalue_fluG_loc(:,k,:),G_sub_eta(:,k,:))
            end if       
            
        end if
    end do
    
    detJ = cellset(index).det_J
    do m =1,4    
        L_sub_global(:,:,m) = -(F_sub_kesi(:,:,m)+G_sub_eta(:,:,m))/detJ(:,:)
    end do
    !
    !do k = 1,nsp
    !    !write(*,"('L',5F15.9)") L_sub_global(k,:,1)
    !    do l =1,nsp     
    !        do m =1,4
    !            if(isnan(F_sub_kesi(k,l,m)))then
    !                write(*,*)'F',xy_coor(cellset(index).nodes(1),:)
    !                write(*,*)index,cellset(index).nearcells(1),k,l,m,cellset(index).Beta  
    !                write(*,"(6F10.5)")cellset(index).fluxF_innerfp(1,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxF_innerfp(2,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxF_innerfp(3,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxF_innerfp(4,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxF_innerfp(5,:,m)
    !                write(*,*)FG_upw(:,:,1)
    !                call print_num_data(T_temp)
    !                stop
    !            end if
    !            if(isnan(G_sub_eta(k,l,m)))then
    !                write(*,*)'G',xy_coor(cellset(index).nodes(1),:)
    !                write(*,*)index,k,l,m 
    !                write(*,*)cellset(index).fluxG_innerfp(k,l,m)
    !                write(*,"(6F10.5)")cellset(index).fluxG_innerfp(1,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxG_innerfp(2,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxG_innerfp(3,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxG_innerfp(4,:,m)
    !                write(*,"(6F10.5)")cellset(index).fluxG_innerfp(5,:,m)
    !                write(*,*)FG_upw(:,:,1)
    !                call print_num_data(T_temp)
    !                stop
    !            end if
    !        end do
    !    end do
    !end do
    ! do m =1,4
    !    if(index == 390 )then
    !        write(*,*)'F',xy_coor(cellset(index).nodes(1),:)
    !        !write(*,*)index,cellset(index).nearcells(1),k,l,m,cellset(index).Beta  
    !        write(*,"(6F20.10)")cellset(index).fluxF_innerfp(1,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxF_innerfp(2,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxF_innerfp(3,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxF_innerfp(4,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxF_innerfp(5,:,m)
    !        write(*,*)FG_upw(2,:,m)
    !        !call print_num_data
    !        !stop
    !    end if
    !    if(index == 390 )then
    !        write(*,*)'G',xy_coor(cellset(index).nodes(1),:)
    !        !write(*,*)index,k,l,m 
    !        !write(*,*)cellset(index).fluxG_innerfp(k,l,m)
    !        write(*,"(6F20.10)")cellset(index).fluxG_innerfp(1,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxG_innerfp(2,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxG_innerfp(3,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxG_innerfp(4,:,m)
    !        write(*,"(6F20.10)")cellset(index).fluxG_innerfp(5,:,m)
    !        write(*,*)FG_upw(3,:,m)
    !        !call print_num_data
    !        !stop
    !    end if
    !end do  
    
    deallocate(FG_upw,O1,O2,Q1,F1,G1,detJ)
    
end subroutine RHS_CPR_NNW
    
subroutine get_subcellFlux(index)
 
    !-----------------------------------------------------------------------------
    !
    !   �������ӵ�Ԫ����ͨ������
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    use RHS_module
    implicit none
    
    real(prec),external :: NonLinear_WCNS,LaI_nPs_deri,LaI_nPs
    
    integer i,j,k,m,l,index,sideIndex,nearCellIndex,LC_sideth,RC_sideth
    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2,x,norm1,norm(2)
    
    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4),u_cell_L(4),u_cell_R(4),xk(3)
    real(prec),dimension(4) :: f_l,f_r,g_l,g_r,conser,ori
    real(prec),dimension(:,:) :: u(3,4),quad_vertex(4,2),ll(4,4),rr(4,4),chara_con(3,4)  
    real(prec),dimension(:,:,:) :: oriDis(nsp+1,2,4)
    real(prec),dimension(:,:,:) :: con_phy(nsp,nsp,4),Ff_phy(nsp,nsp+1,4),Fg_phy(nsp,nsp+1,4),Gf_phy(nsp,nsp+1,4),Gg_phy(nsp,nsp+1,4)
    
    !-----------------------------------------------------------------------------
    !
    !   ������1��ʾ��Ԫ�ڣ�2��ʾ�ٵ�Ԫ
    !
    !-----------------------------------------------------------------------------
    
    O1 = cellset(index).spvalue_ori
    Q1 = cellset(index).spvalue_con_loc
    F1 = cellset(index).spvalue_fluF_loc 
    G1 = cellset(index).spvalue_fluG_loc
    con_phy = cellset(index).spvalue_con
      
    !-----------------------------------------------------------------------------
    !
    !   ��������Ԫ����ϵ�ӭ��ͨ��FG_upw(i,:,:)�����ձߵı��ѭ��ȡ������sideset(sideIndex).fpvalue_upw(:,:)�Ľ���ͨ��
    !
    !-----------------------------------------------------------------------------
    
    do i =1,4
        
        sideIndex = cellset(index).sides(i)
        nearCellIndex = cellset(index).nearcells(i)
  
        if(index == sideset(sideIndex).nearCells(1))then
            
            ! ��Ԫ�ڱߵ���࣬��ͨ����˳��һ�� 
            FG_upw(i,:,:)=sideset(sideIndex).fpvalue_upw(:,:)     
            
        elseif(index == sideset(sideIndex).nearCells(2))then
            
            ! ��Ԫ�ڱߵ���࣬��ͨ����˳����Ҫ�����жϣ�1-2 3-4 ���ڣ�ͨ����˳���෴���߽絥Ԫ���Ǳ߽�ߵ���Ԫ�����ÿ���
            RC_sideth = i
            do j = 1,4
                if(cellset(nearCellIndex).nearcells(j)==index)then
                    LC_sideth = j
                    exit
                end if
            end do
            if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then
                do l = 1,nsp
                    FG_upw(i,l,:)=sideset(sideIndex).fpvalue_upw(nsp+1-l,:) 
                end do   
            else            
                FG_upw(i,:,:) = sideset(sideIndex).fpvalue_upw(:,:) 
            end if
            
            ! �߽�ӭ��ͨ�������ݱߵ���Ԫ��⣬���Ե�����Ԫ�Ǳߵ��ҵ�Ԫʱ���˱����ҵ�Ԫ�ı��ȷ��ӭ��ͨ���ķ���������1-4 2-3����ʱ������ͻ��ȡ�෴��
            if(LC_sideth==RC_sideth .or. (LC_sideth*RC_sideth==4).or.(LC_sideth*RC_sideth==6))then
                FG_upw(i,:,:) = -1.0_prec*FG_upw(i,:,:)
            end if
        else    
            write(*,*)'error'
            stop
        end if     
    end do 
    

    do j = 1,4
        quad_vertex(j,1:2) = xy_coor(cellset(index).nodes(j),1:2)
    end do   
    
    !-----------------------------------------------------------------------------
    !
    !   ������F ͨ����ֵ��nsp+1��ͨ����Riemannͨ����cellset(index).fluxF_innerfp(k,l,:)
    !
    !----------------------------------------------------------------------------- 
    
    oriDis = 0.0_prec
    
    ! �������� Normal vector
    n1 = 1.0_prec
    n2 = 0.0_prec 
    
    ! ��ά���ԣ���������ά��ͨ�������ά
    do k = 1,nsp!��/��
        
        if(cellset(index).Beta_line(k,1)==0 .and. fluxDer_switch == fluxDer_FPs)then
            
            ! CPR������ͨ���������㲻һ��ʱ���Ҳ�ȡLP����
            do l = 1,nsp+1
                call LaI_nPs_arr(FPs(l),SPs(:),O1(k,:,:),nsp,ori(:))
                call Func_ori_to_fluxF(ori(:),Ff_phy(k,l,:)) 
                call Func_ori_to_fluxG(ori(:),Fg_phy(k,l,:))
                kesix = cellset(index).fpMdirect_F(k,l,1)
                kesiy = cellset(index).fpMdirect_F(k,l,2)
                FPdetJ = cellset(index).fpdet_J_F(k,l) 
                cellset(index).fluxF_innerfp(k,l,:) = FPdetJ*(kesix*Ff_phy(k,l,:) + kesiy*Fg_phy(k,l,:))
            end do
        elseif(cellset(index).Beta_line(k,1)==1)then
        
            ! CNNW2����.
            
            !�߽紦������ͨ����������1
            oriDis(2,1,:) = cellset(index).fluxF_innerfp(k,2,:)
            
            ! �ڲ���
            do l = 2,nsp-1                  
                select case(var_type)
                case(ori_type)
                    u = O1(k,l-1:l+1,:)
                    call FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)
                    oriDis(l,2,:)   = u_cell_L
                    oriDis(l+1,1,:) = u_cell_R
                case(character_type)
                    u = O1(k,l-1:l+1,:)
                    call proj_matrix(u(2,:),ll,rr,0)  
                    do j = 1,3
                        call Characteristic_projection(u(j,:),ll,chara_con(j,:)) 
                    end do            
                    call FaceFluxC2NNW2inner(index,l,chara_con,u_cell_L,u_cell_R)
                    call Inverse_Characteristic_projection(u_cell_L,rr,oriDis(l,2,:))
                    call Inverse_Characteristic_projection(u_cell_R,rr,oriDis(l+1,1,:))
                end select           
            end do
            
            ! �߽紦������ͨ����������nsp
            oriDis(nsp,2,:) = cellset(index).fluxF_innerfp(k,nsp,:)
                  
            ! �ӵ�Ԫ����Riemannͨ��
            do l = 2,nsp                 
                
                ! ͨ���㴦��Jacobi����ֱ�Ӷ�����
                kesix = cellset(index).fpMdirect_F(k,l,1)
                kesiy = cellset(index).fpMdirect_F(k,l,2)
                etax  = cellset(index).fpMdirect_F(k,l,3)
                etay  = cellset(index).fpMdirect_F(k,l,4)
                FPdetJ = cellset(index).fpdet_J_F(k,l) 
                
                norm(1) = kesix*n1+etax*n2
                norm(2) = kesiy*n1+etay*n2    
                
                ! Roemman Flux
                call getRiemannFlux(oriDis(l,1,:),oriDis(l,2,:),norm,cellset(index).fluxF_innerfp(k,l,:))
                
                ! ��¼
                cellset(index).fluxF_innerfp(k,l,:) = cellset(index).fluxF_innerfp(k,l,:)*FPdetJ    
            end do
            
            LC_sideth = 4    
            call FP_choose(index,LC_sideth,k,cellset(index).fluxF_innerfp(k,1,:))            
            LC_sideth = 2   
            call FP_choose(index,LC_sideth,k,cellset(index).fluxF_innerfp(k,nsp+1,:))

        end if
    end do   
    

    !-----------------------------------------------------------------------------
    !
    !   ������G ͨ����ֵ��nsp+1��ͨ����Riemannͨ��,cellset(index).fluxG_innerfp(k,l,:)
    !
    !----------------------------------------------------------------------------- 
    
    oriDis = 0.0_prec
    
    !�������� Normal vector
    n1 = 0.0_prec
    n2 = 1.0_prec 
    
    do k = 1,nsp!��/��
        if(cellset(index).Beta_line(k,2)==0 .and. fluxDer_switch == fluxDer_FPs)then
            
            !CPR.�����ͨ���㴦��ͨ��ֵ
            do l = 1,nsp+1
                call LaI_nPs_arr(FPs(l),SPs,O1(:,k,:),nsp,ori(:))
                call Func_ori_to_fluxF(ori(:),Gf_phy(k,l,:)) 
                call Func_ori_to_fluxG(ori(:),Gg_phy(k,l,:))
                etax  = cellset(index).fpMdirect_G(k,l,3)
                etay  = cellset(index).fpMdirect_G(k,l,4)
                FPdetJ = cellset(index).fpdet_J_G(k,l) 
                cellset(index).fluxG_innerfp(k,l,:) = FPdetJ*(etax*Gf_phy(k,l,:) + etay*Gg_phy(k,l,:))                 
            end do
        elseif(cellset(index).Beta_line(k,2)==1)then
            
            ! CNNW2����.
            
            ! �߽紦������ͨ����������1
            oriDis(2,1,:) = cellset(index).fluxG_innerfp(k,2,:)
            
            ! �ڲ���
            do l = 2,nsp-1
                select case(var_type)
                case(ori_type)
                    u = O1(l-1:l+1,k,:)
                    call FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)
                    oriDis(l,2,:)   = u_cell_L
                    oriDis(l+1,1,:) = u_cell_R
                case(character_type)
                    u = O1(l-1:l+1,k,:)
                    call proj_matrix(u(2,:),ll,rr,1)  
                    do j = 1,3
                        call Characteristic_projection(u(j,:),ll,chara_con(j,:))
                    end do            
                    call FaceFluxC2NNW2inner(index,l,chara_con,u_cell_L,u_cell_R)
                    call Inverse_Characteristic_projection(u_cell_L,rr,oriDis(l,2,:))
                    call Inverse_Characteristic_projection(u_cell_R,rr,oriDis(l+1,1,:))
                end select
            end do
            
            ! �߽紦������ͨ����������nsp
            oriDis(nsp,2,:) = cellset(index).fluxG_innerfp(k,nsp,:)

            do l = 2,nsp
                
                !ͨ���㴦��Jacobi����ֱ�Ӷ�����
                kesix = cellset(index).fpMdirect_G(k,l,1)
                kesiy = cellset(index).fpMdirect_G(k,l,2)
                etax  = cellset(index).fpMdirect_G(k,l,3)
                etay  = cellset(index).fpMdirect_G(k,l,4)
                FPdetJ = cellset(index).fpdet_J_G(k,l)
                norm(1) = kesix*n1+etax*n2
                norm(2) = kesiy*n1+etay*n2  
                
                ! Roemman Flux
                call getRiemannFlux(oriDis(l,1,:),oriDis(l,2,:),norm,cellset(index).fluxG_innerfp(k,l,:))
                cellset(index).fluxG_innerfp(k,l,:) = cellset(index).fluxG_innerfp(k,l,:)*FPdetJ
            end do
            
            LC_sideth = 1      
            call FP_choose(index,LC_sideth,k,cellset(index).fluxG_innerfp(k,1,:))  
            LC_sideth = 3 
            call FP_choose(index,LC_sideth,k,cellset(index).fluxG_innerfp(k,nsp+1,:)) 
        end if
    end do

end subroutine get_subcellFlux

subroutine solve_adv_der(flux_fp,flux_sp,flux_Der)!�������
    !���׸�ʽ�Ķ��ײ�����ӣ�����֤����
    use real_precision
    use parameter_setting
    use type_module
    use global_var
    implicit none
    integer :: i,j,k,l
    real(prec),dimension(:,:) :: flux_fp(nsp+1,4),flux_sp(nsp,4),flux_Der(nsp,4)
    real(prec) :: d1,d2,d1_h,d2_h,f_a(4),f_b(4)
    do k = 1,nsp
        d1 = dis_sp_fp(2*k-1)
        d2 = dis_sp_fp(2*k)
        !d1_h = d1*0.5_prec
        !d2_h = d2*0.5_prec
        f_a = (flux_sp(k,:)-flux_fp(k,:))/d1
        f_b = (flux_fp(k+1,:)-flux_sp(k,:))/d2
        flux_Der(k,:) = f_a*d2/(d1+d2) + f_b*d1/(d1+d2)               
    end do
      
end subroutine solve_adv_der
    
subroutine FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)

    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    real(prec),external :: lim,lim2,vanLeerLimiter 

    integer    :: index,nearCellIndex,Cell_sideth,m,l,j,ikesi2,ieta2 
    real(prec) :: kesifp,etafp
    real(prec) :: d1,d2,d3,d4,dd1,dd2,dd3,dd4
    real(prec) :: cc,w1,w2,w3,w4,w5,w6,Fai,varMax,varMin,limA,limB
    
    real(prec),dimension(:,:) :: u(3,4),cell_vertex(4,2)
    real(prec),dimension(2)   :: coor_C,coor_P1,coor_P2,coor_Pfp
    real(prec),dimension(4)   :: uA1,uB1,uA2,uB2,du1,du2,du,u_cell_L,u_cell_R,limu,theta0
    
    ! ͨ��������֮��ľ���
    d1 = dis_sp_fp(2*(l-1))
    d2 = dis_sp_fp(2*(l-1)+1)  
    d3 = dis_sp_fp(2*l)
    d4 = dis_sp_fp(2*l+1)

    ! ������Ȩ
    cc = 1.0_prec!����1������д���������滻
    dd1 = cc/d1
    dd2 = cc/d2    
    dd3 = cc/d3
    dd4 = cc/d4
    
    w1 = dd1/(dd1+dd2)
    w2 = dd2/(dd1+dd2)
    w3 = dd3/(dd3+dd4)
    w4 = dd4/(dd3+dd4)
    
    !ͨ���㴦ֵ
    uA1 = w1*u(1,:)+w2*u(2,:)
    uB1 = w3*u(2,:)+w4*u(3,:)
    
    !��Ԫ����
    !! 1
    w5 = dd2/(dd2+dd3)
    w6 = dd3/(dd2+dd3)
    du1 = (u(2,:)-uA1)/d2
    du2 = (uB1-u(2,:))/d3
    du = w5*du1 + w6*du2 !old
    !! 2
    !du=(uB1-uA1)/(d2+d3)
    
    !! 3
    !du1 = (u(2,:)-u(1,:))/(d1+d2)    !��Ԫ�ڲ�.����ʹ��1,2,3������ݶȣ�
    !du2 = (u(3,:)-u(2,:))/(d3+d4)
    !du = (d3+d4)/(d1+d2+d3+d4)*du1 + (d1+d2)/(d1+d2+d3+d4)*du2 !��Ԫ�ڲ�����ʹ��1,2,3������ݶȣ�
    !
    !!-------------------------------------------
    !do m =1,4
    !    if(du2(m) == 0)then
    !        theta0(m) = 1.0_prec
    !    else
    !        theta0(m) = du1(m)/du2(m)
    !    end if
    !    limu(m) = (abs(theta0(m))+(theta0(m)))/(1.0_prec+abs(theta0(m)))
    !
    !    u_cell_L(m) = u(2,m)-limu(m)*d2*du2(m)
    !    u_cell_R(m) = u(2,m)+limu(m)*d3*du2(m)
    !end do
    
    !----------------------------------------------
    
    !���¼���uA,uB
    uA2 = u(2,:)-du*d2
    uB2 = u(2,:)+du*d3
    
    do m = 1,4
        varMax = maxVal(u(:,m))
        varMin = minVal(u(:,m))
        limA = lim(u(2,m),uA2(m),varMax,varMin)
        limB = lim(u(2,m),uB2(m),varMax,varMin)
        Fai = min(limA,limB)    
        
        if(method_subcell == method_NNW .and. limiter_switch == limiter_close)then
            !�ر������������Բ�ֵ
            Fai = 1.0_prec
        elseif(method_subcell == method_Godunov)then        
            !Godunovs' scheme
            Fai = 0.0_prec
        end if
        
        !��Ԫ��ֵ�ǽ�����ֵ����Ԫ��ֵ�ǽ�����ֵ
        u_cell_L(m) = u(2,m)-Fai*du(m)*d2
        u_cell_R(m) = u(2,m)+Fai*du(m)*d3
        !u_cell_L(m) = u(2,m)-limA*du(m)*d2
        !u_cell_R(m) = u(2,m)+limB*du(m)*d3
    end do

end subroutine FaceFluxC2NNW2inner

subroutine FP_choose(index,sideth,l,FP_upw_temp)
    use real_precision
    use type_module
    use parameter_setting
    implicit none
    integer :: index,nearCellIndex,sideIndex,sideth,near_sideth,i,j,k,l
    real(prec),dimension(:) :: FP_upw_temp(4)
    real(prec) :: norm
    
    nearCellIndex = cellset(index).nearCells(sideth)
    sideIndex = cellset(index).sides(sideth)
   
    if(index == sideset(sideIndex).nearCells(1))then
        
        ! ��Ԫ�ڱߵ����,��ͨ����˳��һ�� 
        FP_upw_temp(:)=sideset(sideIndex).fpvalue_upw(l,:)      
    elseif(index == sideset(sideIndex).nearCells(2))then
        
        do j = 1,4
            if(cellset(nearCellIndex).nearcells(j)==index)then
                near_sideth = j
            end if
        end do         
        
        ! 1-2 3-4 ���ڣ�ͨ����˳���෴���߽絥Ԫ���Ǳ߽�ߵ���Ԫ�����ÿ���  
        if((sideth*near_sideth==2).or.(sideth*near_sideth==12).or.(sideth==near_sideth))then 
            
            ! ͨ����˳���෴ 
            FP_upw_temp(:)=sideset(sideIndex).fpvalue_upw(nsp+1-l,:)                   
        else
            
            ! ͨ����˳��һ�� 
            FP_upw_temp(:) = sideset(sideIndex).fpvalue_upw(l,:)      
        endif
        
        ! �߽�ӭ��ͨ�������ݱߵ���Ԫ��⣬���Ե�����Ԫ�Ǳߵ��ҵ�Ԫʱ���˱����ҵ�Ԫ�ı��ȷ��ӭ��ͨ���ķ���������1-4 2-3����ʱ������ͻ��ȡ�෴��
        if(sideth==near_sideth .or. (sideth*near_sideth==4).or.(sideth*near_sideth==6))then
            FP_upw_temp(:) = -1.0_prec*FP_upw_temp(:)
        end if
    else    
        write(*,*)'error'
        stop
    end if

end subroutine FP_choose
    
subroutine get_fp_fluxDer
!���ͨ���㴦�������
    use RHS_module
    use real_precision
    use parameter_setting
    use global_var
    implicit none
    real(prec),external :: LaI_nPs,LaI_nPs_deri
    integer m,k,l

    if(fluxD_type == LP)then
        do k = 1,nsp
            do l = 1,nsp+1
                do m = 1,4
                    fp_fu_kesi(k,l,m) = LaI_nPs_deri(FPs(l),SPs,F1(k,:,m),nsp)
                    fp_gu_eta(l,k,m) = LaI_nPs_deri(FPs(l),SPs,G1(:,k,m),nsp)
                end do
            end do
        end do
    end if
end subroutine get_fp_fluxDer
    
subroutine get_fluxDer(k,l)
    !���������
    use RHS_module
    use real_precision
    use parameter_setting
    use global_var
    implicit none
    real(prec),external :: LaI_nPs,LaI_nPs_deri
    integer m,k,l
    real(prec) :: r,u,v,e
    real(prec),dimension(4,4) :: Fdq,Gdq
    real(prec),dimension(4)   :: Qder1,Qder2 
    if(fluxD_type == LP)then
        do m = 1,4
            !fu_sub_kesi(m) = LaI_nPs_deri(SPs(l),SPs,F1(k,:,m),nsp)
            !gu_sub_eta(m)  = LaI_nPs_deri(SPs(k),SPs,G1(:,l,m),nsp)
            fu_sub_kesi(m) = LaI_nPs(SPs(l),FPs,fp_fu_kesi(k,:,m),nsp+1)
            gu_sub_eta(m)  = LaI_nPs(SPs(k),FPs,fp_gu_eta(:,l,m),nsp+1)
        end do
    elseif(fluxD_type == CR)then
        r = Q1(k,l,1)
        u = Q1(k,l,2)/r
        v = Q1(k,l,3)/r
        e = Q1(k,l,4)
        !A(U)
        Fdq(1,1) = 0.0_prec
        Fdq(1,2) = 1.0_prec
        Fdq(1,3) = 0.0_prec
        Fdq(1,4) = 0.0_prec
        Fdq(2,1) = 0.5_prec*((gamma-3.0_prec)*u**2+(gamma-1.0_prec)*v**2)
        Fdq(2,2) = (3.0_prec-gamma)*u
        Fdq(2,3) = (1.0_prec-gamma)*v
        Fdq(2,4) = gamma-1.0_prec
        Fdq(3,1) = -u*v
        Fdq(3,2) = v
        Fdq(3,3) = u
        Fdq(3,4) = 0.0_prec
        Fdq(4,1) = -gamma*e*u/r+(gamma-1.0_prec)*(u**3+u*v**2)
        Fdq(4,2) = gamma*e/r-0.5_prec*(gamma-1.0_prec)*(3.0_prec*u**2+v**2)
        Fdq(4,3) = (1.0_prec-gamma)*u*v
        Fdq(4,4) = gamma*u
        !B(U)
        Gdq(1,1) = 0.0_prec
        Gdq(1,2) = 0.0_prec
        Gdq(1,3) = 1.0_prec
        Gdq(1,4) = 0.0_prec
        Gdq(2,1) = -u*v
        Gdq(2,2) = v
        Gdq(2,3) = u
        Gdq(2,4) = 0.0_prec
        Gdq(3,1) = 0.5_prec*((gamma-3.0_prec)*v**2+(gamma-1.0_prec)*u**2)
        Gdq(3,2) = (1.0_prec-gamma)*u
        Gdq(3,3) = (3.0_prec-gamma)*v
        Gdq(3,4) = gamma-1.0_prec
        Gdq(4,1) = -gamma*e*v/r+(gamma-1.0_prec)*(v**3+v*u**2)
        Gdq(4,2) = (1-gamma)*u*v
        Gdq(4,3) = gamma*e/r-0.5_prec*(gamma-1.0_prec)*(3.0_prec*v**2+u**2)
        Gdq(4,4) = gamma*v
        do m = 1,4
            Qder1(m) = LaI_nPs_deri(SPs(l),SPs,Q1(k,:,m),nsp)
            Qder2(m) = LaI_nPs_deri(SPs(k),SPs,Q1(:,l,m),nsp)
        end do
        do m = 1,4
            fu_sub_kesi(m) = dot_product(Fdq(m,:),Qder1(:))
            gu_sub_eta(m)  = dot_product(Gdq(m,:),Qder2(:))
        end do
    end if

end subroutine get_fluxDer
