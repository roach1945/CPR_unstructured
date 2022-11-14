subroutine Face_Flux_upw
     
    !-----------------------------------------------------------------------------
    !
    !   ���������㵥Ԫ����ϵ�ӭ��ͨ��
    !   ������������Ԫ�غϵ�Ԫ�߽�ͨ������һ��
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    use global_var
    use type_module
    use real_precision
    use bc_module
    implicit none
    
    integer :: i,j,m,k,l,RC_l,LC_sideth,RC_sideth           !----LC_sideth  ��Ԫ�ĵ�-th����
    integer :: indexCellL,indexCellR                        !----indexCellL ��Ԫ�ı��

    real(prec) :: FPdetJ,kesix,kesiy,etax,etay,n1,n2
    real(prec),dimension(1:2) :: norm
    real(prec),dimension(1:4) :: ori1,ori2
    real(prec),dimension(nsp,4) :: oriL,oriR
    real(prec),dimension(:)   :: coor_C(2),FPJacobi(4),FPdirect(4)
    real(prec),dimension(:,:) :: quad_vertex(4,2)
    
    do i = 1,nsides   
             
        !-----------------------------------------------------------------------------
        !
        !   ������׼�����̣�ȷ��LC_sideth��RC_sideth
        !
        !-----------------------------------------------------------------------------
     
        indexCellL = sideset(i).nearcells(1)
        indexCellR = sideset(i).nearcells(2)

        if(indexCellL==0)then
            
            ! indexCellLΪ0��˵���Ǳ߽��ϵıߣ�����ԪΪ���ⵥԪ������ݱ߽�����ȡ��Ԫ
            do j = 1,4
                if(cellset(indexCellR).sides(j)==i)then
                    indexCellL = cellset(indexCellR).nearcells(j)
                    RC_sideth = j
                    exit
                end if
            end do  
            do j = 1,4
                if(cellset(indexCellL).sides(j)==i)then
                    LC_sideth = j
                    exit
                end if
            end do  

        else if(indexCellR==0)then
            
            ! indexCellRΪ0��˵���Ǳ߽��ϵıߣ����ҵ�ԪΪ���ⵥԪ������ݱ߽�����ȡ�ҵ�Ԫ
            do j = 1,4
                if(cellset(indexCellL).sides(j)==i)then
                    indexCellR = cellset(indexCellL).nearcells(j)
                    LC_sideth = j
                    exit
                end if
            end do
            do j = 1,4
                if(cellset(indexCellR).sides(j)==i)then
                    RC_sideth = j
                    exit
                end if
            end do
        else
            
            !�ж�����������Ԫ�ĵڼ�����ߣ����ڵ�Ԫ��������ͨ����������Ԫ��⣬�ڹ��������������ڵ�Ԫӭ��ͨ����һ�µ�
            do j = 1,4
                if(cellset(indexCellL).sides(j) == i)then
                    LC_sideth = j
                end if
                if(cellset(indexCellR).sides(j) == i)then
                    RC_sideth = j
                end if
            end do      
        end if

        !-----------------------------------------------------------------------------
        !
        !   ������Riemann Flux �����ͨ�����ӭ��ͨ�����
        !
        !-----------------------------------------------------------------------------

        do l = 1,nsp    
            
            ! ������ڵ�Ԫ��Ӧ����/�б��
            call get_RC_l(LC_sideth,RC_sideth,l,RC_l)
            
            ! ��ȡ��������ֵ,ԭʼ����       
            call get_oriL_weight(indexCellL,indexCellR,LC_sideth,RC_sideth,l,RC_l,oriL(l,:))
            call get_oriR_weight(indexCellL,indexCellR,LC_sideth,RC_sideth,l,RC_l,oriR(l,:))
            
            ! �������� Normal vector
            if(LC_sideth == 1)then
                n1 = 0.0_prec
                n2 = -1.0_prec         
            else if(LC_sideth == 2)then
                n1 = 1.0_prec
                n2 = 0.0_prec
            else if(LC_sideth == 3)then
                n1 = 0.0_prec
                n2 = 1.0_prec
            else if(LC_sideth == 4)then
                n1 = -1.0_prec
                n2 = 0.0_prec
            end if
            
            ! ��Ԫ�����ͨ���㴦��Jacobi����ֱ�Ӷ�����
            if(LC_sideth==1)then
                coor_C(1) = SPs(l)
                coor_C(2) = -1.0_prec
            else if(LC_sideth == 2)then
                coor_C(1) = 1.0_prec
                coor_C(2) = SPs(l)
            else if(LC_sideth == 3)then
                coor_C(1) = SPs(l)
                coor_C(2) = 1.0_prec
            else if(LC_sideth == 4)then
                coor_C(1) = -1.0_prec
                coor_C(2) = SPs(l)
            end if     
            
            do j = 1,4
                quad_vertex(j,:) = xy_coor(cellset(indexCellL).nodes(j),:)
            end do    
            
            call solve_Jacobi(quad_vertex,coor_C,FPJacobi,FPdirect,FPdetJ)
            
            kesix = FPdirect(1)
            kesiy = FPdirect(2)
            etax  = FPdirect(3)
            etay  = FPdirect(4)

            norm(1) = kesix*n1+etax*n2
            norm(2) = kesiy*n1+etay*n2
            
            ! Roemman Flux
            call getRiemannFlux(oriL(l,:),oriR(l,:),norm,sideset(i).fpvalue_upw(l,:))         
            sideset(i).fpvalue_upw(l,:) = sideset(i).fpvalue_upw(l,:)*FPdetJ    
            
        end do
                        
        if(LC_sideth == 1 .or. LC_sideth == 4)then
            sideset(i).fpvalue_upw(:,:) = -1.0_prec*sideset(i).fpvalue_upw(:,:)
        end if  
        
    end do
    
end subroutine Face_Flux_upw
    
subroutine get_RC_l(LC_sideth,RC_sideth,l,RC_l)
 
    !-----------------------------------------------------------------------------
    !
    !   ������ȷ�����ڵ�Ԫ��Ӧ���������ڵ�Ԫ�еı��
    !   ������1-2 3-4 ���ڣ�ͨ����˳���෴���߽絥Ԫ���Ǳ߽�ߵ���Ԫ�����ÿ���
    !
    !-----------------------------------------------------------------------------
 
    use parameter_setting
    implicit none
    
    integer :: LC_sideth,RC_sideth,l,RC_l
    
    if((LC_sideth*RC_sideth==2).or.(LC_sideth*RC_sideth==12).or.(LC_sideth==RC_sideth))then
        RC_l = nsp+1-l
    else           
        RC_l = l
    end if

end subroutine get_RC_l

!!## ��Ȩ��ֵ
subroutine get_oriL_weight(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,oriL)
        
    !-----------------------------------------------------------------------------
    !
    !   ��������ⵥԪ���ͨ���㴦Riemannͨ����ֵ
    !   ������NNW��ֵ��ʹ��ԭʼ����������������
    !         Larange��ֵ��NNW��ֵ ��Ȩ�������ǲ����ơ��ο�ʯ��Ȩ˶ʿ���ġ��ǽṹ����CPR�������ӵ�Ԫ���Ƽ�����
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none 
    
    integer ::index,nearCellIndex,LC_sideth,RC_sideth,l,m,k,RC_l
    real(prec),dimension(:,:) :: u(3,4),conL(4),ll(4,4),rr(4,4),chara_con(3,4),ui_all(8,4)
    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4)
    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),oriL(4),chara_oriL(4),chara_conL(4),uA(4)
    real(prec) :: realDis1,realDis2,detJ_fp,dP1P2
    integer :: beta_dim,beta_dim_near
    real(prec),dimension(:) :: oriL_CPR(4),oriL_CNNW2(4)
    real(prec) :: smoother
    
    ! ���ҵ�Ԫ�ϵ�ԭʼ����   
    varOriL(:,:,:) = cellset(index).spvalue_ori(:,:,:)
    varOriR(:,:,:) = cellset(nearCellIndex).spvalue_ori(:,:,:)
    
    !-----------------------------------------------------------------------------
    !
    !   ������smoother �ǵ�Ԫ�Ĺ⻬�����ӣ�oriL = smoother*oriL_CNNW2 + (1.0_prec - smoother)*oriL_CPR
    !         smoother = 0 ����CPR��Larange��ֵ��smoother = 1 ����CNNW2��NNW��ֵ
    !         ��MDH������smoother,�Ӷ����Ʋ�ȡ�ļ����ʽ��indicator_MDH3 �� d_smooth = 0 ����ʹ�ü�Ȩ�ķ���
    !         beta_dim = 1 ������Ԫ�����ⵥԪ,beta_dim_near = 1�����ڵ�Ԫ�����ⵥԪ
    !
    !-----------------------------------------------------------------------------
    
    beta_dim_near = 0
    if(LC_sideth == 2 .or.LC_sideth == 4)then
        beta_dim = cellset(index).Beta_line(l,1)       
        smoother = cellset(index).Smooth_line_x(l,1)
    elseif(LC_sideth == 1 .or.LC_sideth == 3)then
        beta_dim = cellset(index).Beta_line(l,2)
        smoother = cellset(index).Smooth_line_y(l,1)
    end if
    if(RC_sideth == 2 .or.RC_sideth == 4)then  
        beta_dim_near = cellset(nearCellIndex).Beta_line(RC_l,1)       
    elseif(RC_sideth == 1 .or.RC_sideth == 3)then
        beta_dim_near  = cellset(nearCellIndex).Beta_line(RC_l,2)                        
    end if

    oriL_CPR   = 0.0_prec
    oriL_CNNW2 = 0.0_prec
    
    ! �⻬��Ԫ��Larange��ֵ
    if(beta_dim == 0 .or. (beta_dim == 0 .and. beta_dim_near == 1))then
        if(LC_sideth == 1)then                                      
            call LaI_nPs_arr(kesi_l,SPs,varOriL(:,l,:),nsp,oriL)
        else if(LC_sideth == 2)then
            call LaI_nPs_arr(kesi_r,SPs,varOriL(l,:,:),nsp,oriL)
        else if(LC_sideth == 3)then
            call LaI_nPs_arr(kesi_r,SPs,varOriL(:,l,:),nsp,oriL)
        else if(LC_sideth == 4)then
            call LaI_nPs_arr(kesi_l,SPs,varOriL(l,:,:),nsp,oriL)
        end if
       oriL_CPR = oriL
    end if
    
    ! ���ⵥԪ��NNW��ֵ
    if(beta_dim == 1 .or. (beta_dim == 0 .and. beta_dim_near == 1))then 
        if(var_type == ori_type .OR. var_type == character_type)then
            if(LC_sideth == 1)then    
                u(2,:) = varOriL(1,l,:)
                u(3,:) = varOriL(2,l,:)
            else if(LC_sideth == 2)then
                u(2,:) = varOriL(l,nsp,:)
                u(3,:) = varOriL(l,nsp-1,:)
            else if(LC_sideth == 3)then
                u(2,:) = varOriL(nsp,l,:)
                u(3,:) = varOriL(nsp-1,l,:)
            else if(LC_sideth == 4)then
                u(2,:) = varOriL(l,1,:)
                u(3,:) = varOriL(l,2,:)
            end if     

            call get_nearCellInfo(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,varOriL,varOriR,u,uA,ui_all,dP1P2)

            call FaceFluxC2NNW2(index,nearCellIndex,LC_sideth,l,u,uA,ui_all,dP1P2,oriL(:),innerDisOri)
            
            oriL_CNNW2 = oriL
 
        end if

        ! �ӵ�Ԫ�Ȱѽ���ԭʼ������������֮��絥Ԫ��������һ��  
        if(LC_sideth==1)then        
            cellset(index).fluxG_innerfp(l,2,:) = innerDisOri
        elseif(LC_sideth==2)then
            cellset(index).fluxF_innerfp(l,nsp,:) = innerDisOri
        elseif(LC_sideth==3)then
            cellset(index).fluxG_innerfp(l,nsp,:) = innerDisOri
        elseif(LC_sideth==4)then
            cellset(index).fluxF_innerfp(l,2,:) = innerDisOri
        end if

        if( beta_dim_near == 1 .and. beta_dim == 0 )then
            oriL = smoother*oriL_CNNW2 + (1.0_prec - smoother)*oriL_CPR
        end if
            
    end if
    
end subroutine get_oriL_weight

subroutine get_oriR_weight(index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,oriR)

    !-----------------------------------------------------------------------------
    !
    !   ��������ⵥԪ���ͨ���㴦Riemannͨ����ֵ
    !   �������ǽṹ�ı������в����򣬸��ӳ������˼������Ӧ�Ĳ�߱�Ź�ϵ
    !         ����һ�ַ�������һ���м�������ڵ�Ԫ���С���ת��ʹ֮��Ӧ������1-3��2-4��Ӧ������
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none 
    
    integer :: index,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,m,j,k
    real(prec),dimension(:,:) :: u(3,4),ll(4,4),rr(4,4),chara_u(3,4),chara_con(3,4),ui_all(8,4)
    real(prec),dimension(:,:,:) :: varOriL(nsp,nsp,4),varOriR(nsp,nsp,4),varConL(nsp,nsp,4),varConR(nsp,nsp,4)
    real(prec),external :: LaI_nPs
    real(prec),dimension(:) :: coor_P1(2),innerDisOri(4),innerDisCon(4),chara_innerDisCon(4),conR(4),oriR(4),chara_conR(4),uA(4)
    real(prec) :: realDis,detJ_fp,dP1P2
    integer :: beta_dim,beta_dim_near
    real(prec),dimension(:) :: oriR_CPR(4),oriR_CNNW2(4)
    real(prec) :: smoother
    
    ! ȡ�����ҵ�Ԫ�ϵ�ԭʼ����   
    varOriL(:,:,:) = cellset(index).spvalue_ori(:,:,:)
    varOriR(:,:,:) = cellset(nearCellIndex).spvalue_ori(:,:,:)

    if(RC_sideth == 2 .or.RC_sideth == 4)then
        beta_dim = cellset(nearCellIndex).Beta_line(RC_l,1)
        smoother = cellset(nearCellIndex).Smooth_line_x(RC_l,1)
    elseif(RC_sideth == 1 .or. RC_sideth == 3)then
        beta_dim = cellset(nearCellIndex).Beta_line(RC_l,2)
        smoother = cellset(nearCellIndex).Smooth_line_y(RC_l,1)
    end if
    
    if(LC_sideth == 2 .or.LC_sideth == 4)then  
        beta_dim_near = cellset(index).Beta_line(l,1)        
    elseif(LC_sideth == 1 .or.LC_sideth == 3)then
        beta_dim_near  = cellset(index).Beta_line(l,2)              
    end if

    oriR_CPR   = 0.0_prec
    oriR_CNNW2 = 0.0_prec
    
    ! �⻬��Ԫ��Larange��ֵ
    if(beta_dim == 0 .or. (beta_dim == 0 .and. beta_dim_near == 1))then  
        if(RC_sideth == 1)then
            call LaI_nPs_arr(kesi_l,SPs,varOriR(:,RC_l,:),nsp,oriR)
        elseif(RC_sideth == 2)then
            call LaI_nPs_arr(kesi_r,SPs,varOriR(RC_l,:,:),nsp,oriR)
        else if(RC_sideth == 3)then
            call LaI_nPs_arr(kesi_r,SPs,varOriR(:,RC_l,:),nsp,oriR)
        else if(RC_sideth == 4)then
            call LaI_nPs_arr(kesi_l,SPs,varOriR(RC_l,:,:),nsp,oriR)
        end if 
        oriR_CPR   = oriR
    end if
    
    ! ���ⵥԪ��NNW��ֵ
    if(beta_dim == 1 .or. (beta_dim == 0 .and. beta_dim_near == 1))then
        if(var_type==ori_type .OR. var_type==character_type)then         
            if(RC_sideth == 1)then
                u(2,:) = varOriR(1,RC_l,:)
                u(3,:) = varOriR(2,RC_l,:)   
            else if(RC_sideth == 2)then
                u(2,:) = varOriR(RC_l,nsp,:)
                u(3,:) = varOriR(RC_l,nsp-1,:)               
            else if(RC_sideth == 3)then
                u(2,:) = varOriR(nsp,RC_l,:)
                u(3,:) = varOriR(nsp-1,RC_l,:)
            else if(RC_sideth == 4)then
                u(2,:) = varOriR(RC_l,1,:)
                u(3,:) = varOriR(RC_l,2,:)
            end if           
            
            call get_nearCellInfo(nearCellIndex,index,RC_sideth,LC_sideth,RC_l,l,varOriR,varOriL,u(:,:),uA,ui_all,dP1P2)
            
            call FaceFluxC2NNW2(nearCellIndex,index,RC_sideth,RC_l,u,uA,ui_all,dP1P2,oriR(:),innerDisOri)
            
            oriR_CNNW2   = oriR 
     
        end if        

        if(RC_sideth==1)then      
            cellset(nearCellIndex).fluxG_innerfp(RC_l,2,:) = innerDisOri
        elseif(RC_sideth==2)then
            cellset(nearCellIndex).fluxF_innerfp(RC_l,nsp,:) = innerDisOri
        elseif(RC_sideth==3)then
            cellset(nearCellIndex).fluxG_innerfp(RC_l,nsp,:) = innerDisOri
        elseif(RC_sideth==4)then
            cellset(nearCellIndex).fluxF_innerfp(RC_l,2,:) = innerDisOri
        end if

    end if
    
    if( beta_dim_near == 1 .and. beta_dim == 0)then
        oriR = smoother*oriR_CNNW2 + (1.0_prec - smoother)*oriR_CPR
    end if
    
end subroutine get_oriR_weight
    
subroutine get_wall_oriR(indexCell,LC_sideth,oriL,oriR)    

    !-----------------------------------------------------------------------------
    !
    !   ���������������Ĵ���
    !   ������**
    !   ���ߣ�gqShi 
    !   ��ʷ�����������Ĵ�������ȡ���ⵥԪ�ķ��������ǻ��ڵ��� //2021.12.17
    !
    !-----------------------------------------------------------------------------

    use global_var
    use parameter_setting
    use type_module
    implicit none 
    integer :: indexCell,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,m,j,k,vertexth1,vertexth2
    real(prec),dimension(:) :: oriL(4),oriR(4)
    real(prec) :: x1,y1,x2,y2
    real(prec) :: theta_x12,theta_VaX,theta_VaX12,theta_symVaX  ! ��1 2��x�� �ļнǣ���A���ٶ���x��ļнǣ��ԳƵ�A'��x��ļн�
    real(prec) :: r,u,v,p,r_sym,u_sym,v_sym,p_sym,VV
    vertexth1 = LC_sideth
    vertexth2 = LC_sideth+1
    if(vertexth1 == 4) vertexth2 = 1
    x1 = xy_coor(cellset(indexCell).nodes(vertexth1),1)
    y1 = xy_coor(cellset(indexCell).nodes(vertexth1),2)
    x2 = xy_coor(cellset(indexCell).nodes(vertexth2),1)
    y2 = xy_coor(cellset(indexCell).nodes(vertexth2),2)
    ! ����ԳƵ㴦�������ش�ֱ���淽���ٶȴ�С��ȣ������෴
    
    !write(*,*)theta_x12,theta_VaX,theta_VaX12
    r = oriL(1)
    u = oriL(2)
    v = oriL(3)
    p = oriL(4)
    !write(*,*)u,v
    VV = sqrt(u**2 + v**2)
    theta_x12 = atan((y2-y1),(x2-x1))
    theta_VaX = atan(v,u)
    theta_symVaX = 2.0_prec*theta_x12 - theta_VaX  !!!! -(theta_VaX - theta_x12) + theta_x12
    !write(*,*) theta_x12,theta_VaX,theta_symVaX
    r_sym = r
    u_sym = VV * cos(theta_symVaX)
    v_sym = VV * sin(theta_symVaX)
    p_sym = p
    !write(*,*)'sym',u_sym,v_sym
    oriR(1) = r_sym
    oriR(2) = u_sym
    oriR(3) = v_sym
    oriR(4) = p_sym
    !write(*,*) '1',oriR
    
    !
    !u = oriL(2)
    !v = oriL(3)
    !VV = sqrt(u**2 + v**2)
    !theta_x12 = atan((y2-y1),(x2-x1))
    !theta_VaX = atan(v,u)
    !theta_VaX12 = theta_x12 - theta_VaX
    !oriL(2) = VV*cos(theta_VaX12)*cos(theta_x12 )
    !oriL(3) = VV*cos(theta_VaX12)*sin(theta_x12 )
    !!write(*,*)theta_x12,theta_VaX,theta_VaX12
    !r = oriL(1)
    !u = oriL(2)
    !v = oriL(3)
    !p = oriL(4)
    !write(*,*)u,v
    !!write(*,*)'sym',u_sym,v_sym
    !oriR(1) = r
    !oriR(2) = u
    !oriR(3) = v
    !oriR(4) = p
end subroutine get_wall_oriR    
    
subroutine get_wall_oriR2(indexCell,LC_sideth,oriL,oriR)    
 
    !-----------------------------------------------------------------------------
    !
    !   ���������������Ĵ���
    !   ������**
    !   ���ߣ�gqShi 
    !   ��ʷ�����������Ĵ�������ȡ���ⵥԪ�ķ��������ǻ��ڵ��� //2021.12.17
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none 
    integer :: indexCell,nearCellIndex,LC_sideth,RC_sideth,l,RC_l,m,j,k,vertexth1,vertexth2
    real(prec),dimension(:) :: oriL(4),oriR(4)
    real(prec) :: x1,y1,x2,y2,xa,ya,x,y,n1,n2
    real(prec) :: theta_x12,theta_VaX,theta_VaX12,theta_symVaX  ! ��1 2��x�� �ļнǣ���A���ٶ���x��ļнǣ��ԳƵ�A'��x��ļн�
    real(prec) :: r,u,v,p,r_sym,u_sym,v_sym,p_sym,VV,c_b,s,VVbp_k
    real(prec),dimension(2) :: normal,VV_b,VVp_b
    !�����ⷨ������
    vertexth1 = LC_sideth
    vertexth2 = LC_sideth+1
    if(vertexth1 == 4) vertexth2 = 1
    x1 = xy_coor(cellset(indexCell).nodes(vertexth1),1)
    y1 = xy_coor(cellset(indexCell).nodes(vertexth1),2)
    x2 = xy_coor(cellset(indexCell).nodes(vertexth2),1)
    y2 = xy_coor(cellset(indexCell).nodes(vertexth2),2)
    do j = 1,4
        if(j .NE. vertexth1 .and. j .NE.vertexth2)then
            xa = xy_coor(cellset(indexCell).nodes(j),1)
            ya = xy_coor(cellset(indexCell).nodes(j),2)
            exit
        end if       
    end do   
    call sym_Point(x1,y1,x2,y2,xa,ya,x,y) 

    n1 = (x-xa)/sqrt((x-xa)**2 + (y-ya)**2)
    n2 = (y-ya)/sqrt((x-xa)**2 + (y-ya)**2)
    normal(1) = n1
    normal(2) = n2
    !!����������
    r = oriL(1)
    u = oriL(2)
    v = oriL(3)
    p = oriL(4)
    !write(*,*)r,u,v,p
    VVp_b(1) = u
    VVp_b(2) = v
    VVbp_k = dot_product(VVp_b,normal)!VVp_b(1)*normal(1) + VVp_b(2)*normal(2)
    VV_b = VVp_b - normal*VVbp_k
    !����
    u = VV_b(1)!VVp_b(1) - n1*VVbp_k
    v = VV_b(2)!VVp_b(2) - n2*VVbp_k
    !write(*,*)u,v
    !�����������
    c_b = sqrt(gamma*p/r) + 0.5_prec*gamma1*VVbp_k
    s = p/(r**gamma)
    r = (c_b**2/gamma/s)**(1.0_prec/gamma1)
    p = s*r**gamma
    !write(*,*)'0',normal
    !write(*,*)'1',VVbp_k,dot_product(VVp_b,normal)
    !write(*,*)c_b
    !write(*,*)indexCell,r
    !write(*,*)xa,ya
    !write(*,*)x,y
    !write(*,*)indexCell,n1,n2,n1**2+n2**2
    !stop
    oriL(1) = r
    oriL(2) = u
    oriL(3) = v
    oriL(4) = p
    
    oriR = oriL
end subroutine get_wall_oriR2        
!!---��һ�³��򣬴ӱ���Ԫ���ڵ�Ԫֱ�ӻ�ȡ�߽�ͨ���㴦��ֵ��֮�����ͬԭ��ʽ-----------------------------------------------------------------------------------------------------------------
!!PS �� ԭ�����ǻ�����ڵ�Ԫ�Ľ��ԭʼ����ֵ�ͽ���ͨ�����ʵ�ʾ���
subroutine get_nearCellInfo(index,nearCellIndex,Cell_sideth,nearCell_sideth,l,near_l,CellOri,nearCellOri,ui,uA,ui_all,dP1P2)

     
    !-----------------------------------------------------------------------------
    !
    !   �����������������Ϣ
    !   ����������NNW ��ֵģ��u(1)��uA������ռ����
    !         uA�Ĳ�ֵģ���������֣�һ����������Ȩ��������Χ6������Ȩ���������ҵ�Ԫ�������Ȼ��ƽ������4�㣩
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none 
    
    integer :: index,nearCellIndex,Cell_sideth,nearCell_sideth,l,near_l,m,i,j,near_l_t,l_t
    real(prec),dimension(:) :: uA(4),u1(4),u2(4),u3(4),u4(4),u5(4),u6(4),u1_2(4),u2_2(4),uA_1(4),uA_2(4)
    real(prec),dimension(:)   :: coor_C(2),coor_P(2),coor_P2(2),coor_Pfp1(2),coor_Pfp2(2)                   !----��׼��Ԫ��ʵ�ʵ�Ԫ
    real(prec),dimension(:,:) :: cell_vertex(4,2),cell_vertex_near(4,2)                                     !----ʵ�ʵ�Ԫ���㣬��׼��Ԫ��֪ 
    real(prec),dimension(:,:) :: ui(3,4),ui_all(8,4)
    real(prec),dimension(:,:,:) :: nearCellOri(nsp,nsp,4),CellOri(nsp,nsp,4)
    real(prec) :: realDis,kesifp,etafp,det_u1,det_u2,det_A
    real(prec) :: det_u1_2,det_u2_2
    real(prec) :: cc,d1,d2,d3,d4,d5,d6,dd1,dd2,dd3,dd4,dd5,dd6,dd,w1,w2,w3,w4,w5,w6,dP1P2                   !----������ʾ�ڱ���Ԫ��ż����ʾ�ڵ�Ԫ
    integer :: ikesi,ieta,ikesi2,ieta2,count_ui
    
    realDis = 0.0_prec
    ui_all = 0.0_prec
    count_ui = 1
    d1 = 0.0_prec;d2 = 0.0_prec;d3 = 0.0_prec;d4 = 0.0_prec;d5 = 0.0_prec;d6 = 0.0_prec
    cc = 1.0_prec !----����1������д���������滻
         
    !-----------------------------------------------------------------------------
    !
    !   �������ڵ�Ԫ
    !
    !-----------------------------------------------------------------------------
     
    if(nearCell_sideth == 1)then    !�ڵ�Ԫ
        ikesi = 1;            ieta  = near_l
        kesifp = SPs_local(ikesi,ieta,1); etafp = -1.0_prec
    elseif(nearCell_sideth == 2)then
        ikesi = near_l;            ieta  = nsp
        kesifp = 1.0_prec; etafp = SPs_local(ikesi,ieta,2)
    else if(nearCell_sideth == 3)then
        ikesi = nsp;            ieta  = near_l
        kesifp = SPs_local(ikesi,ieta,1); etafp = 1.0_prec
    else if(nearCell_sideth == 4)then
        ikesi = near_l;            ieta  = 1
        kesifp = -1.0_prec; etafp = SPs_local(ikesi,ieta,2) 
    end if
    det_u1 = cellset(nearCellIndex).det_J(ikesi,ieta)!u1����Jacobi
    
    ! ȡ�����ڵ�Ԫ��ģ���
    u1(:) = nearCellOri(ikesi,ieta,:)    
    do i =1,8
        ui_all(i,:) = u1(:)! �״θ�ֵ ���綨����ui�ķ�Χ����Ϊ�е�����ò������еĵ㣬����ȡ��ֵʱ�������⡣
    end do
    count_ui = count_ui+1
    
    ! �����1����������                            
    do i = 1,4
        cell_vertex_near(i,:) = xy_coor(cellset(nearCellIndex).nodes(i),:)
    end do
    coor_C = SPs_local(ikesi,ieta,:)
    call quad_C2Phy(cell_vertex_near,coor_C,coor_P)
    coor_C(1) = kesifp;    coor_C(2) = etafp
    call quad_C2Phy(cell_vertex_near,coor_C,coor_Pfp1)
    d1 = sqrt((coor_P(1)-coor_Pfp1(1))**2+(coor_P(2)-coor_Pfp1(2))**2)   
             
    !-----------------------------------------------------------------------------
    !
    !   ����������Ԫ
    !
    !-----------------------------------------------------------------------------
    
    ! ��2 (ikesi2��ieta2) �ǵ�2�ڵ�Ԫ�ı�ţ�(kesifp,etafp)ͨ���������
    if(Cell_sideth==1)then      !��1�ϵ�ͨ����
        ikesi2  = 1.0_prec;     ieta2 = l;
        kesifp = SPs_local(ikesi2,ieta2,1);    etafp = -1.0_prec;
    elseif(Cell_sideth==2)then  !��2�ϵ�ͨ����
        ikesi2 = l;     ieta2 = nsp;
        kesifp = 1.0_prec;     etafp = SPs_local(ikesi2,ieta2,2);
    elseif(Cell_sideth==3)then  !��3�ϵ�ͨ����
        ikesi2 = nsp;   ieta2 = l;
        kesifp = SPs_local(ikesi2,ieta2,1);    etafp = 1.0_prec;
    elseif(Cell_sideth==4)then  !��4�ϵ�ͨ����
        ikesi2 = l;     ieta2 = 1.0_prec;
        kesifp = -1.0_prec;     etafp = SPs_local(ikesi2,ieta2,2);
    end if

    u2(:) = CellOri(ikesi2,ieta2,:)    

    ui_all(count_ui,:) = u2(:);count_ui = count_ui+1
    do j = 1,4
        cell_vertex(j,:) = xy_coor(cellset(index).nodes(j),:)
    end do
    coor_C = SPs_local(ikesi2,ieta2,:)
    call quad_C2Phy(cell_vertex,coor_C,coor_P2)
    coor_C(1) = kesifp; coor_C(2) = etafp
    call quad_C2Phy(cell_vertex,coor_C,coor_Pfp2)    
    d2 = sqrt((coor_Pfp2(1)-coor_P2(1))**2 + (coor_Pfp2(2)-coor_P2(2))**2 )  
    
    ! 1,2֮����������
    dP1P2 = sqrt((coor_P(1)-coor_P2(1))**2 + (coor_P(2)-coor_P2(2))**2 ) 
    
             
    !-----------------------------------------------------------------------------
    !
    !   ������NNW��õ�Ԫ�߽�ͨ���㴦��ʼֵ�ķ���  TPs_weight_NNW = 0, MPs_weight_NNW = 1, TTPsL_ave_NNW = 2
    !
    !-----------------------------------------------------------------------------
      
    if(method_ori_fp == TPs_weight_NNW )then    
        
        ! ���㷴�����Ȩ
        dd1 = cc/d1
        dd2 = cc/d2 
        w1 = dd1/(dd1+dd2)
        w2 = dd2/(dd1+dd2)
        uA = w1*u1(:)+w2*u2(:)       
        
    elseif(method_ori_fp == MPs_weight_NNW)then
        
        ! ��㷴�����Ȩ���ñ���Ԫ��ΧһȦ��
        near_l_t = near_l-1
        if(near_l==1)then
            near_l_t = near_l+1
        elseif(near_l==nsp)then
            near_l_t = near_l-1
        end if
        if(nearCell_sideth == 1)then
            ikesi = 1;            ieta  = near_l_t
        elseif(nearCell_sideth == 2)then
            ikesi = near_l_t;            ieta  = nsp
        else if(nearCell_sideth == 3)then
            ikesi = nsp;            ieta  = near_l_t
        else if(nearCell_sideth == 4)then
            ikesi = near_l_t;            ieta  = 1
        end if  
        u3(:) = nearCellOri(ikesi,ieta,:)   
        ui_all(count_ui,:) = u3(:);count_ui = count_ui+1
        coor_C = SPs_local(ikesi,ieta,:)
        call quad_C2Phy(cell_vertex_near,coor_C,coor_P)
        d3 = sqrt((coor_P(1)-coor_Pfp1(1))**2+(coor_P(2)-coor_Pfp1(2))**2)   
        if(near_l/=1 .and. near_l/=nsp)then   
            near_l_t = near_l+1
            if(nearCell_sideth == 1)then
                ikesi = 1;            ieta  = near_l_t
            elseif(nearCell_sideth == 2)then
                ikesi = near_l_t;            ieta  = nsp
            else if(nearCell_sideth == 3)then
                ikesi = nsp;            ieta  = near_l_t
            else if(nearCell_sideth == 4)then
                ikesi = near_l_t;            ieta  = 1
            end if  
            u5(:) = nearCellOri(ikesi,ieta,:) 
            ui_all(count_ui,:) = u5(:);count_ui = count_ui+1
            coor_C = SPs_local(ikesi,ieta,:)
            call quad_C2Phy(cell_vertex_near,coor_C,coor_P)
            d5 = sqrt((coor_P(1)-coor_Pfp1(1))**2+(coor_P(2)-coor_Pfp1(2))**2)   
        end if
        l_t = l-1
        if(l==1)then
            l_t = l+1
        elseif(l==nsp)then
            l_t = l-1
        end if
        if(Cell_sideth==1)then      !��1�ϵ�ͨ����
            ikesi2  = 1;     ieta2 = l_t;
        elseif(Cell_sideth==2)then  !��2�ϵ�ͨ����
            ikesi2 = l_t;     ieta2 = nsp;
        elseif(Cell_sideth==3)then  !��3�ϵ�ͨ����
            ikesi2 = nsp;   ieta2 = l_t;
        elseif(Cell_sideth==4)then  !��4�ϵ�ͨ����
            ikesi2 = l_t;     ieta2 = 1;
        end if 
        u4(:) = CellOri(ikesi2,ieta2,:)
        ui_all(count_ui,:) = u4(:);count_ui = count_ui+1
        coor_C = SPs_local(ikesi2,ieta2,:)
        call quad_C2Phy(cell_vertex,coor_C,coor_P)
        d4 = sqrt((coor_P(1)-coor_Pfp2(1))**2+(coor_P(2)-coor_Pfp2(2))**2)   
        if(l/=1 .and. l/=nsp)then   
            l_t = l+1
            if(Cell_sideth==1)then      !��1�ϵ�ͨ����
                ikesi2  = 1;     ieta2 = l_t;
            elseif(Cell_sideth==2)then  !��2�ϵ�ͨ����
                ikesi2 = l_t;     ieta2 = nsp;
            elseif(Cell_sideth==3)then  !��3�ϵ�ͨ����
                ikesi2 = nsp;   ieta2 = l_t;
            elseif(Cell_sideth==4)then  !��4�ϵ�ͨ����
                ikesi2 = l_t;     ieta2 = 1;
            end if 
            u6(:) = CellOri(ikesi2,ieta2,:)   
            ui_all(count_ui,:) = u6(:);count_ui = count_ui+1
            coor_C = SPs_local(ikesi2,ieta2,:)
            call quad_C2Phy(cell_vertex,coor_C,coor_P)
            d6 = sqrt((coor_P(1)-coor_Pfp2(1))**2+(coor_P(2)-coor_Pfp2(2))**2)    
        end if
        
        dd1 = cc/d1;dd2 = cc/d2;dd3 = cc/d3;dd4 = cc/d4;dd5 = cc/d5;dd6 = cc/d6;
        if(l/=1 .and. l/=nsp.and.near_l/=1.and.near_l/=nsp)then
            dd = dd1+dd2+dd3+dd4+dd5+dd6
            w1 = dd1/dd;w2 = dd2/dd;w3 = dd3/dd;w4 = dd4/dd;w5 = dd5/dd;w6 = dd6/dd;
            uA = w1*u1(:)+w2*u2(:)+w3*u3(:)+w4*u4(:)+w5*u5(:)+w6*u6(:)
        else
            dd = dd1+dd2+dd3+dd4
            w1 = dd1/dd;w2 = dd2/dd;w3 = dd3/dd;w4 = dd4/dd;
            uA = w1*u1(:)+w2*u2(:)+w3*u3(:)+w4*u4(:)
        end if
    elseif(method_ori_fp == TTPsL_ave_NNW)then
        
        ! ����Ԫ�������Բ�ֵ��ƽ���� ÿ����Ԫ���������ֵ��ͨ���㴦��Ȼ��ƽ��
        if(nearCell_sideth == 1)then    !�ڵ�Ԫ
            ikesi = 1+1;            ieta  = near_l
        elseif(nearCell_sideth == 2)then
            ikesi = near_l;            ieta  = nsp-1
        else if(nearCell_sideth == 3)then
            ikesi = nsp-1;            ieta  = near_l
        else if(nearCell_sideth == 4)then
            ikesi = near_l;            ieta  = 1+1
        end if
        !ȡ�����ڵ�Ԫ��ģ���
        u1_2(:) = nearCellOri(ikesi,ieta,:) 
        
        if(Cell_sideth==1)then      !��1�ϵ�ͨ����
            ikesi2  = 1+1;     ieta2 = l;
        elseif(Cell_sideth==2)then  !��2�ϵ�ͨ����
            ikesi2 = l;     ieta2 = nsp-1;
        elseif(Cell_sideth==3)then  !��3�ϵ�ͨ����
            ikesi2 = nsp-1;   ieta2 = l;
        elseif(Cell_sideth==4)then  !��4�ϵ�ͨ����
            ikesi2 = l;     ieta2 = 1+1;
        end if
        u2_2(:) = CellOri(ikesi2,ieta2,:)    
    
        uA_1 = dis_sp_fp(1)/dis_sp_fp(2)*(u1(:) - u1_2(:)) +  u1(:) 
        uA_2 = dis_sp_fp(1)/dis_sp_fp(2)*(u2(:) - u2_2(:)) +  u2(:) 
        !uA_1 = dis_sp_fp(1)/dis_sp_fp(2)*(u2(:) - u1(:)) +  u2(:) 
        !uA_2 = dis_sp_fp(1)/dis_sp_fp(2)*(u2_2(:) - u1_2(:)) +  u2_2(:) 
        uA = (uA_1 + uA_2)*0.5_prec
        ui_all(count_ui,:) = u1_2(:);count_ui = count_ui+1
        ui_all(count_ui,:) = u2_2(:);count_ui = count_ui+1
    end if
    ui(1,:) = u1(:) !
    ui_all(count_ui,:) = ui(2,:);count_ui = count_ui+1
    ui_all(count_ui,:) = ui(3,:);
    
end subroutine get_nearCellInfo

subroutine FaceFluxC2NNW2(index,nearCellIndex,Cell_sideth,l,u,uA,ui_all,dP1P2,u_cell_L,u_cell_R)
     
    !-----------------------------------------------------------------------------
    !
    !   ��������Ԫ�߽紦��NNW��ֵ
    !   ������CNNW2��ʽ������ҵ���ġ�uA1�Ѿ����
    !   ���ߣ�gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    real(prec),dimension(:,:) :: u(3,4),u_keep(3,4),u_char(3,4),ui_all(8,4),ui_all_char(8,4)
    real(prec),dimension(2) :: coor_C,coor_P1,coor_P2,coor_Pfp
    real(prec),dimension(4) :: uA1,uB1,uA2,uB2,du1,du2,du,u_cell_L,u_cell_R,u2_phy,uA,limu ,theta0
    real(prec) :: d1,d2,d3,d4,dd1,dd2,dd3,dd4
    real(prec) :: cc,w1,w2,w3,w4,w5,w6,Fai,varMax,varMin,limA,limB   
    real(prec),external :: lim,vanLeerLimiter
    real(prec),dimension(:,:) :: cell_vertex(4,2)   !ʵ�ʵ�Ԫ����
    integer :: index,nearCellIndex,Cell_sideth,m,l,j,ikesi2,ieta2,k
    real(prec) :: kesifp,etafp,realDis,det_A,detJ_u2,dP1P2
    real(prec),dimension(:,:) :: ll(4,4),rr(4,4),ruvwp(4),ll2(4,4),rr2(4,4)
    
    ! ͨ��������֮��ľ���
    d2 = dis_sp_fp(1)   
    d3 = dis_sp_fp(2)
    d4 = dis_sp_fp(3)
    
    ! ������Ȩ
    cc = 1.0_prec
    dd2 = cc/d2
    dd3 = cc/d3
    dd4 = cc/d4
    
    w3 = dd3/(dd3+dd4)
    w4 = dd4/(dd3+dd4)    
    uA1 = uA
    
    ! ���ѡ������������ֵ������Ҫ���������任
    if(var_type == character_type)then  
        u_keep = u
        call proj_matrix(u(2,:),ll,rr,Cell_sideth)  
        do k = 1,3
            call Characteristic_projection(u(k,:),ll,u(k,:)) !----��k���� ����u������,�������������
        end do   
        call Characteristic_projection(uA1,ll,uA1)
    end if
    
    ! ͨ���㴦ֵ 
    uB1 = w3*u(2,:)+w4*u(3,:)
    
    ! ������Ȩ�͵�Ԫ����A,B�ֱ������ĵ�2���ݶ�
    w5 = dd2/(dd2+dd3)
    w6 = dd3/(dd2+dd3)
    du1 = (u(2,:)-uA1)/d2
    du2 = (uB1-u(2,:))/d3
    
    
    !-----------------------------------------------------------------------------
    !
    !   ������ԭʼ��CNNW2���Ʒ���
    !
    !-----------------------------------------------------------------------------
    
    ! 2�����ݶȣ����ּ��㷽��     
    du = w5*du1 + w6*du2    !----old����Ȩ
    !du=(uB1-uA1)/(d2+d3)   !----��Ԫ����ֵ�Ĳ�ֵ/��Ԫ����
    
    ! ���¼���uA,uB
    uA2 = u(2,:)-du*d2
    uB2 = u(2,:)+du*d3
    
    ! ȡ��ֵʱ��u1���ȡ���ڴ�����ռ�������ù��ı任������ڼ���ռ䣬��ʹ�ü���ռ��������Ҫ���Ǳ任
    if(var_type == character_type)then
        do k = 1,8
            call Characteristic_projection(ui_all(k,:),ll,ui_all(k,:))
        end do
    end if

    do m = 1,4
        varMax = maxVal(ui_all(:,m))
        varMin = minVal(ui_all(:,m))
        limA = lim(u(2,m),uA2(m),varMax,varMin)
        limB = lim(u(2,m),uB2(m),varMax,varMin)
        Fai = min(limA,limB) 
        
        if(method_subcell == method_NNW .and. limiter_switch == limiter_close)then
            ! �ر������������Բ�ֵ
            Fai = 1.0_prec
        elseif(method_subcell == method_Godunov)then        
            ! Godunovs' scheme
            Fai = 0.0_prec
        end if
        
        !��Ԫ��ֵ�ǽ�����ֵ����Ԫ��ֵ�ǽ�����ֵ
        u_cell_L(m) = u(2,m)-Fai*du(m)*d2
        u_cell_R(m) = u(2,m)+Fai*du(m)*d3      
        !u_cell_L(m) = u(2,m)-limA*du(m)*d2
        !u_cell_R(m) = u(2,m)+limB*du(m)*d3 
        
                 
        !-----------------------------------------------------------------------------
        !
        !   ��������һ��������
        !
        !-----------------------------------------------------------------------------

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
    end do
    
    ! ���ѡ������������ֵ������Ҫ���������任
    if(var_type == character_type)then  
        call Inverse_Characteristic_projection(u_cell_L(:),rr,u_cell_L(:))
        call Inverse_Characteristic_projection(u_cell_R(:),rr,u_cell_R(:))
        u = u_keep
    end if
    
    !if(method_ori_fp == MPs_weight_NNW .or. method_ori_fp == TTPsL_ave_NNW)then
    !    u_cell_L = uA
    !end if
    
end subroutine FaceFluxC2NNW2


