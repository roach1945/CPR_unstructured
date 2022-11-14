
    
subroutine RHS_CPR(index,Lsub)
    !这部分是求解一个cell里所有解点的右端项，放在Lsub里
    use real_precision
    use parameter_setting,only:nsp
    use RHS_module
    implicit none
    integer index,i,j,k,n1,n2    
    real(prec),dimension(:,:,:) :: Lsub(nsp,nsp,4)
    
    allocate(FG_upw(4,nsp,4),O1(nsp,nsp,4),O2(nsp,nsp,4),Q1(nsp,nsp,4),F1(nsp,nsp,4),G1(nsp,nsp,4),detJ(nsp,nsp))!O1,O2原始变量
    allocate(fp_fu_kesi(nsp,nsp+1,4),fp_gu_eta(nsp+1,nsp,4))
    
    call get_local(index)           !取出计算右端项需要的值
    call get_RHS_dif(index,Lsub)    !右端项
    
    deallocate(FG_upw,O1,O2,Q1,F1,G1,detJ)
    deallocate(fp_fu_kesi,fp_gu_eta)
end subroutine RHS_CPR

subroutine get_local(index)
    use parameter_setting,only:nsp
    use type_module
    use RHS_module
    use global_var
    implicit none
    integer index,i,j,k,l,m,n
    integer :: sideIndex,nearCellIndex,LC_sideth,RC_sideth  
    real(prec) :: norm
    real(prec) :: detJJ(nsp,nsp) 
    real(prec) :: kesix,kesiy,etax,etay
    
    detJJ  = cellset(index).det_J  
    !1表示单元内，2表示临单元
    O1 = cellset(index).spvalue_ori
    Q1 = cellset(index).spvalue_con_loc
    F1 = cellset(index).spvalue_fluF_loc 
    G1 = cellset(index).spvalue_fluG_loc
    do i =1,4 !side
        sideIndex = cellset(index).sides(i)
        nearCellIndex = cellset(index).nearcells(i)
        
        !计算域下边界间断通量
        !1表示单元内，2表示临单元
        O2 = cellset(nearCellIndex).spvalue_ori           
        if(index == sideset(sideIndex).nearCells(1))then
            FG_upw(i,:,:)=sideset(sideIndex).fpvalue_upw(:,:)      !若是单元在边的左侧则通量点顺序一致 
        elseif(index == sideset(sideIndex).nearCells(2))then
            RC_sideth = i
            do j = 1,4
                if(cellset(nearCellIndex).nearcells(j)==index)then
                    LC_sideth = j
                    exit
                end if
            end do
            
            if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then!1-2 3-4 相邻，通量点顺序相反。边界单元都是边界边的左单元，不用考虑
                do l = 1,nsp
                    FG_upw(i,l,:)=sideset(sideIndex).fpvalue_upw(nsp+1-l,:)      !通量点顺序相反
                end do   
                if(LC_sideth==RC_sideth)then
                    !write(*,*) 'RHS_CPR.f90 line66'
                    norm = -1.0_prec
                    FG_upw(i,:,:) = norm*FG_upw(i,:,:)
                end if
                !write(*,*)LC_sideth,RC_sideth
            else           
            !边界迎风通量是依据边的左单元求解，所以当所求单元是边的右单元时，此边在右单元的编号确定迎风通量的法向向量，1-4 2-3相邻时有所冲突，取相反数
                if((LC_sideth*RC_sideth==4).or.(LC_sideth*RC_sideth==6))then
                    norm = -1.0_prec
                else
                    norm = 1.0_prec
                end if
                FG_upw(i,:,:) = norm*sideset(sideIndex).fpvalue_upw(:,:)      !通量点顺序一致 
            end if
        else    
            write(*,*)'error'
            stop
        end if
        
    end do 

end subroutine get_local

subroutine get_RHS_dif(index,L_sub_global)
    use global_var
    use parameter_setting
    use type_module
    use RHS_module
    implicit none
    integer k,l,m,index
    real(prec),external :: LaI_nPs,LaI_nPs_deri,gl_sub_kesi,gr_sub_kesi
    real(prec),dimension(4) :: f_l,f_r,g_l,g_r
    real(prec),dimension(:,:,:) :: L_sub_global(nsp,nsp,4),F_sub_kesi(nsp,nsp,4),G_sub_eta(nsp,nsp,4)
    
    !求解通量点处对流项导数
    call get_fp_fluxDer
    detJ = cellset(index).det_J
    do k = 1,nsp
        do l =1,nsp      
            do m =1,4              
                f_l(m) = LaI_nPs(kesi_l,SPs,F1(k,:,m),nsp)!通量插值
                f_r(m) = LaI_nPs(kesi_r,SPs,F1(k,:,m),nsp)!同一单元
                g_l(m) = LaI_nPs(kesi_l,SPs,G1(:,l,m),nsp)
                g_r(m) = LaI_nPs(kesi_r,SPs,G1(:,l,m),nsp)
            end do  
            !call LaI_nPs_arr(kesi_l,SPs,F1(k,:,:),nsp,f_l)!通量插值
            !call LaI_nPs_arr(kesi_r,SPs,F1(k,:,:),nsp,f_r)!同一单元
            !call LaI_nPs_arr(kesi_l,SPs,G1(:,l,:),nsp,g_l)
            !call LaI_nPs_arr(kesi_r,SPs,G1(:,l,:),nsp,g_r)
            !求解对流项导数
            call get_fluxDer(k,l)
                    
            !F_sub_kesi(k,l,:) = fu_sub_kesi(:)+(FG_upw(4,k,:)-f_l(:))*gl_sub_kesi(SPs(l)) + (FG_upw(2,k,:)-f_r(:))*gr_sub_kesi(SPs(l))
            !G_sub_eta(k,l,:)  = gu_sub_eta(:) +(FG_upw(1,l,:)-g_l(:))*gl_sub_kesi(SPs(k)) + (FG_upw(3,l,:)-g_r(:))*gr_sub_kesi(SPs(k)) 
            F_sub_kesi(k,l,:) = fu_sub_kesi(:)+(FG_upw(4,k,:)-f_l(:))*gl_coe(l) + (FG_upw(2,k,:)-f_r(:))*gr_coe(l)
            G_sub_eta(k,l,:)  = gu_sub_eta(:) +(FG_upw(1,l,:)-g_l(:))*gl_coe(k) + (FG_upw(3,l,:)-g_r(:))*gr_coe(k)             
        end do
    end do
    do m =1,4
        L_sub_global(:,:,m) = -(F_sub_kesi(:,:,m)+G_sub_eta(:,:,m))/detJ(:,:)
    end do
end subroutine

subroutine get_fp_fluxDer
!求解通量点处对流项导数
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
                    fp_fu_kesi(k,l,m) = LaI_nPs_deri(LPs(l),SPs,F1(k,:,m),nsp)
                    fp_gu_eta(l,k,m) = LaI_nPs_deri(LPs(l),SPs,G1(:,k,m),nsp)
                end do
            end do
        end do
    end if
end subroutine get_fp_fluxDer
subroutine get_fluxDer(k,l)
    !求解对流项导数
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
            fu_sub_kesi(m) = LaI_nPs(SPs(l),LPs,fp_fu_kesi(k,:,m),nsp+1)
            gu_sub_eta(m)  = LaI_nPs(SPs(k),LPs,fp_gu_eta(:,l,m),nsp+1)
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

function gl_sub_kesi(kesi)
    !!gl = (-1)**K/2*(P_k-P_k-1),Pk is Legendre polynomial,左端矫正函数
    !!K=1,gl = -1/2*(x-1),gl_sub_kesi = -1/2
    !!K=2,gl = 1/4*(3*x**2-2*x-1),gl_sub_kesi = 1/4*(6*x-2)
    !!K=3,gl = -1/4*(5*x**3-3*x**2-3*x+1),gl_sub_kesi = -1/4*(15*x**2-6*x-3)
    !!K=4,gl = 1/16*(35*x**4-20*x**3-30*x**2+12*x+3),gl_sub_kesi = 1/16*(140*x**3-60*x**2-60*x+12)
    !!K=5,gl = -1/16*(63*x**5-35*x**4-70*x**3+30*x**2+15*x-3),gl_sub_kesi = -1/16*(315*x**4-140*x**3-210*x**2+60*x+15)
    use real_precision
    use parameter_setting ,only:nsp
    implicit none
    real(prec) :: gl_sub_kesi,kesi
    if (nsp .EQ. 1)then
        gl_sub_kesi = -0.5_prec
    elseif (nsp .EQ. 2)then
        gl_sub_kesi = 0.25*(6.0_prec*kesi - 2.0_prec)
    elseif (nsp .EQ. 3)then
        gl_sub_kesi = -0.25_prec*(15.0_prec*kesi**2.0_prec-6.0_prec*kesi-3.0_prec)
    elseif(nsp .EQ. 4)then
        gl_sub_kesi = (140.0_prec*kesi**3-60.0_prec*kesi**2-60.0_prec*kesi+12.0_prec)/16.0_prec
    elseif(nsp .EQ. 5)then
        gl_sub_kesi = (315.0_prec*kesi**4-140.0_prec*kesi**3-210.0_prec*kesi**2+60.0_prec*kesi+15.0_prec)*(-1.0_prec)/16.0_prec
    endif
    return
end function gl_sub_kesi

function gr_sub_kesi(kesi)
    !!gr = 1/2*(P_k+P_k-1),右端矫正函数
    !!K=1,gr = 1/2*(x+1),gr_sub_kesi = 1/2
    !!K=2,gr = 1/4*(3*x**2+2*x-1),gl_sub_kesi = 1/4*(6*x+2)
    !!K=3,gr = 1/4*(5*x**3+3*x**2-3*x-1),gr_sub_kesi = 1/4*(15*x**2+6*x-3)
    !!K=4,gr = 1/16*(35*x**4+20*x**3-30*x**2-12*x+3),gr_sub_kesi = 1/16*(140*x**3+60*x**2-60*x-12)
    !!K=5,gr = 1/16*(63*x**5+35*x**4-70*x**3-30*x**2+15*x+3),gr_sub_kesi = 1/16*(315*x**4+140*x**3-210*x**2-60*x+15)
    use real_precision
    use parameter_setting ,only:nsp
    implicit none
    real(prec) :: gr_sub_kesi,kesi

    if (nsp .EQ. 1)then
        gr_sub_kesi = 0.5_prec
    elseif (nsp .EQ. 2)then
        gr_sub_kesi = 0.25*(6.0_prec*kesi + 2.0_prec)    
    elseif (nsp .EQ. 3)then
        gr_sub_kesi = 0.25_prec*(15.0_prec*kesi**2.0_prec+6.0_prec*kesi-3.0_prec)
    elseif(nsp .EQ. 4)then
        gr_sub_kesi = (140.0_prec*kesi**3+60.0_prec*kesi**2-60.0_prec*kesi-12.0_prec)/16.0_prec
    elseif(nsp .EQ. 5)then
        gr_sub_kesi = (315.0_prec*kesi**4+140.0_prec*kesi**3-210.0_prec*kesi**2-60.0_prec*kesi+15.0_prec)/16.0_prec
    endif
    return
end function gr_sub_kesi