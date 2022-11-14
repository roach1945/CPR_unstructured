subroutine RHS_NNW2(index,L_sub_global)
    use global_var
    use parameter_setting
    use type_module
    implicit none
    integer k,l,m,index,i,j
    real(prec),dimension(:,:,:) :: L_sub_global(nsp,nsp,4),F_sub_kesi(nsp,nsp,4),G_sub_eta(nsp,nsp,4)
    real(prec),dimension(:,:) :: detJ(nsp,nsp)
    real(prec),external :: LaI_nPs_deri
    detJ = cellset(index).det_J

    call get_subcellFlux2(index)

    do k = 1,nsp
        do l =1,nsp     
            !求解点解对流项导数            
            F_sub_kesi(k,l,:) = (cellset(index).fluxF_innerfp(k,l+1,:)-cellset(index).fluxF_innerfp(k,l,:))/(dis_sp_fp(2*l-1)+dis_sp_fp(2*l))
            G_sub_eta(k,l,:)  = (cellset(index).fluxG_innerfp(l,k+1,:)-cellset(index).fluxG_innerfp(l,k,:))/(dis_sp_fp(2*k-1)+dis_sp_fp(2*k))  

            do m =1,4
                if(isnan(F_sub_kesi(k,l,m)))then
                    write(*,*)'F'
                    write(*,*)index,cellset(index).nearcells(1),cellset(522).nearcells(3),k,l,m  
                    write(*,"(4F10.5)")cellset(index).fluxF_innerfp(k,:,1)
                    write(*,"(4F10.5)")cellset(index).fluxF_innerfp(k,:,2)
                    write(*,"(4F10.5)")cellset(index).fluxF_innerfp(k,:,3)
                    write(*,"(4F10.5)")cellset(index).fluxF_innerfp(k,:,4)
                    call print_num_data
                    stop
                end if
                if(isnan(G_sub_eta(k,l,m)))then
                    write(*,*)'G'
                    write(*,*)index,k,l,m 
                    write(*,*)cellset(index).fluxG_innerfp(l,k,m)
                    write(*,"(4F10.5)")cellset(index).fluxG_innerfp(l,:,1)
                    write(*,"(4F10.5)")cellset(index).fluxG_innerfp(l,:,2)
                    write(*,"(4F10.5)")cellset(index).fluxG_innerfp(l,:,3)
                    write(*,"(4F10.5)")cellset(index).fluxG_innerfp(l,:,4)
                    write(*,*)cellset(index).fluxG_innerfp
                    call print_num_data
                    stop
                end if
            end do
        end do
    end do
    do m =1,4
        L_sub_global(:,:,m) = -(F_sub_kesi(:,:,m)+G_sub_eta(:,:,m))/detJ(:,:)
    end do
    !if(index==729)then
    !    write(*,*)'Fkesi 1'
    !    write(*,"(5F15.6)")F_sub_kesi(1,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(2,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(3,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(4,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(5,:,1)
    !    write(*,*)'Geta 1'
    !    write(*,"(5F15.6)")G_sub_eta(1,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(2,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(3,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(4,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(5,:,1)
    !    write(*,*)'J 1'
    !    write(*,"(5F15.6)")detJ(1,:)
    !    write(*,"(5F15.6)")detJ(2,:)
    !    write(*,"(5F15.6)")detJ(3,:)
    !    write(*,"(5F15.6)")detJ(4,:)
    !    write(*,"(5F15.6)")detJ(5,:)
    !    write(*,*)'L_sub 1'
    !    write(*,"(5F15.6)")L_sub_global(1,:,1)
    !    write(*,"(5F15.6)")L_sub_global(2,:,1)
    !    write(*,"(5F15.6)")L_sub_global(3,:,1)
    !    write(*,"(5F15.6)")L_sub_global(4,:,1)
    !    write(*,"(5F15.6)")L_sub_global(5,:,1)
    !    !pause
    !end if
    !if(index==730)then
    !    write(*,*)'Fkesi 2'
    !    write(*,"(5F15.6)")F_sub_kesi(1,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(2,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(3,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(4,:,1)
    !    write(*,"(5F15.6)")F_sub_kesi(5,:,1)
    !    write(*,*)'Geta 2'
    !    write(*,"(5F15.6)")G_sub_eta(1,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(2,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(3,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(4,:,1)
    !    write(*,"(5F15.6)")G_sub_eta(5,:,1)
    !    write(*,*)'L_sub 2'
    !    write(*,"(5F15.6)")L_sub_global(1,:,1)
    !    write(*,"(5F15.6)")L_sub_global(2,:,1)
    !    write(*,"(5F15.6)")L_sub_global(3,:,1)
    !    write(*,"(5F15.6)")L_sub_global(4,:,1)
    !    write(*,"(5F15.6)")L_sub_global(5,:,1)
    !    !pause
    !end if
    
    !write(*,*) F_sub_kesi(1,1,1),G_sub_eta(1,1,1)
end subroutine RHS_NNW2

subroutine get_subcellFlux2(index)
    use global_var
    use parameter_setting
    use type_module
    implicit none
    integer i,j,k,m,l,index,sideIndex,nearCellIndex,LC_sideth,RC_sideth
    real(prec),dimension(:,:) :: u(3,4),quad_vertex(4,2),ll(4,4),rr(4,4),chara_con(3,4) 
    real(prec),dimension(4) :: f_l,f_r,g_l,g_r
    real(prec),dimension(:,:,:) :: oriDis(nsp+1,2,4)!间断通量
    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4),u_cell_L(4),u_cell_R(4),norm(2),xk(3)
    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2,x
    real(prec),external :: NonLinear_WCNS
    !write(*,*) index,'*********************************************************************'
    !取出左四边形单元顶点坐标值
    do j = 1,4
        do i = 1,2
            quad_vertex(j,i) = xy_coor(cellset(index).nodes(j),i)
        end do
    end do   

    !F没有用计算域上的变量值，
    oriDis = 0.0_prec
    !法向向量 Normal vector
    n1 = 1.0_prec
    n2 = 0.0_prec 
    do k = 1,nsp!行/列
        oriDis(2,1,:) = cellset(index).fluxF_innerfp(k,2,:)!边界处是黎曼通量，跳过点1
        do l = 2,nsp-1!内部点                  
            select case(var_type)
            case(ori_type)
            !原始变量
                u = cellset(index).spvalue_ori(k,l-1:l+1,:)
                call FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)
                oriDis(l,2,:)   = u_cell_L
                oriDis(l+1,1,:) = u_cell_R

            case(con_type)
            !守恒变量
                u = cellset(index).spvalue_con(k,l-1:l+1,:)
                call FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)
             
                call Func_con_to_ori(u_cell_L,oriDis(l,2,:))
                call Func_con_to_ori(u_cell_R,oriDis(l+1,1,:))
            case(character_type)
            !特征变量
                u = cellset(index).spvalue_ori(k,l-1:l+1,:)
                call proj_matrix(u(2,:),ll,rr,0)  
                do j = 1,3
                    call Characteristic_projection(u(j,:),ll,chara_con(j,:)) !第k个点 
                end do            
                call FaceFluxC2NNW2inner(index,l,chara_con,u_cell_L,u_cell_R)
                call Inverse_Characteristic_projection(u_cell_L,rr,oriDis(l,2,:))
                call Inverse_Characteristic_projection(u_cell_R,rr,oriDis(l+1,1,:))
            end select           
        end do
        oriDis(nsp,2,:) = cellset(index).fluxF_innerfp(k,nsp,:)
        
        do l = 2,nsp                 
            !通量点处的Jacobi矩阵，直接度量等
            coor_C(1) = FPs(l)
            coor_C(2) = SPs(k)             
            call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)     
            !write(*,*)FPdetJ
            kesix = FPdirect(1)
            kesiy = FPdirect(2)
            etax  = FPdirect(3)
            etay  = FPdirect(4)
            norm(1) = kesix*n1+etax*n2
            norm(2) = kesiy*n1+etay*n2          
            !Roemman Flux(  Lax Friedrich Flux)
            call laxf_flux(oriDis(l,1,:),oriDis(l,2,:),norm,cellset(index).fluxF_innerfp(k,l,:))
            cellset(index).fluxF_innerfp(k,l,:) = cellset(index).fluxF_innerfp(k,l,:)*FPdetJ
        end do
        
        LC_sideth = 4             
        do l = 1,nsp
            call FP_choose(index,LC_sideth,l,cellset(index).fluxF_innerfp(l,1,:))            
        end do
        
        LC_sideth = 2   
        do l = 1,nsp
            call FP_choose(index,LC_sideth,l,cellset(index).fluxF_innerfp(l,nsp+1,:))            
        end do
        !cellset(index).fluxF_innerfp(:,1,:) = sideset(cellset(index).sides(4)).fpvalue_upw(:,:) 
        !cellset(index).fluxF_innerfp(:,nsp+1,:) = sideset(cellset(index).sides(2)).fpvalue_upw(:,:)         
    end do
    !write(*,*)cellset(index).fluxF_innerfp(:,:,:)
    !G
    oriDis = 0.0_prec
    !法向向量 Normal vector
    n1 = 0.0_prec
    n2 = 1.0_prec 
    do k = 1,nsp!行/列
        oriDis(2,1,:) = cellset(index).fluxG_innerfp(k,2,:)!边界处是黎曼通量，跳过点1
        do l = 2,nsp-1!内部点
            select case(var_type)
            case(ori_type)
                !原始变量
                u = cellset(index).spvalue_ori(l-1:l+1,k,:)
                call FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)
                oriDis(l,2,:)   = u_cell_L
                oriDis(l+1,1,:) = u_cell_R
                !WCNS
                !do m = 1,4
                !    x = FPs(l)
                !    xk = SPs(l-1:l+1)
                !    oriDis(l,2,m) = NonLinear_WCNS(x,xk,u(:,m))
                !    x = FPs(l+1)
                !    oriDis(l+1,1,m) = NonLinear_WCNS(x,xk,u(:,m))
                !end do
            case(con_type)
                !守恒变量
                u = cellset(index).spvalue_con(l-1:l+1,k,:)
                call FaceFluxC2NNW2inner(index,l,u,u_cell_L,u_cell_R)
                call Func_con_to_ori(u_cell_L,oriDis(l,2,:))
                call Func_con_to_ori(u_cell_R,oriDis(l+1,1,:))            
            case(character_type)
                !特征变量
                u = cellset(index).spvalue_ori(l-1:l+1,k,:)
                call proj_matrix(u(2,:),ll,rr,1)  
                do j = 1,3
                    call Characteristic_projection(u(j,:),ll,chara_con(j,:)) !第k个点 
                end do            
                call FaceFluxC2NNW2inner(index,l,chara_con,u_cell_L,u_cell_R)
                call Inverse_Characteristic_projection(u_cell_L,rr,oriDis(l,2,:))
                call Inverse_Characteristic_projection(u_cell_R,rr,oriDis(l+1,1,:))
            end select

        end do
        oriDis(nsp,2,:) = cellset(index).fluxG_innerfp(k,nsp,:)
        !write(*,*)'oriDis',oriDis(2:nsp,:,:)
        do l = 2,nsp                    
            !通量点处的Jacobi矩阵，直接度量等
            coor_C(1) = SPs(k)
            coor_C(2) = FPs(l)             
            call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)           
            kesix = FPdirect(1)
            kesiy = FPdirect(2)
            etax  = FPdirect(3)
            etay  = FPdirect(4)
            norm(1) = kesix*n1+etax*n2
            norm(2) = kesiy*n1+etay*n2          
            !Roemman Flux(  Lax Friedrich Flux)
            call laxf_flux(oriDis(l,1,:),oriDis(l,2,:),norm,cellset(index).fluxG_innerfp(k,l,:))
            !write(*,*)'555',oriDis(l,1,:),oriDis(l,2,:),norm,cellset(index).fluxG_innerfp(k,l,:)
            cellset(index).fluxG_innerfp(k,l,:) = cellset(index).fluxG_innerfp(k,l,:)*FPdetJ
        end do
        LC_sideth = 1            
        do l = 1,nsp
            call FP_choose(index,LC_sideth,l,cellset(index).fluxG_innerfp(l,1,:))            
        end do
        
        LC_sideth = 3 
        do l = 1,nsp
            call FP_choose(index,LC_sideth,l,cellset(index).fluxG_innerfp(l,nsp+1,:))            
        end do
        !cellset(index).fluxG_innerfp(:,1,:) = sideset(cellset(index).sides(1)).fpvalue_upw(:,:) 
        !cellset(index).fluxG_innerfp(:,nsp+1,:) = sideset(cellset(index).sides(3)).fpvalue_upw(:,:) 
    end do
    !do k = 1,nsp
    !    do l = 1,nsp+1
    !        write(*,*) k,l
    !        write(*,"(4F15.6)")cellset(index).fluxG_innerfp(k,l,:)
    !        write(*,"(8F15.6)") oriDis(l,:,:)
    !    end do
    !end do
end subroutine get_subcellFlux2


