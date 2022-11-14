subroutine indicator_JST
    !Refs:
        !! Matthias Sonntag, Stephan Schmidt, Nicolas R. Gauger,
        !! Shape derivatives for the compressible Navier–Stokes equations in variational form,
        !! Journal of Computational and Applied Mathematics,
        !! Volume 296,
        !! 2016,
        !! Pages 334-351,
        !! ISSN 0377-0427,
        !! https://doi.org/10.1016/j.cam.2015.09.010.
    
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,Beta,Beta_line_x,Beta_line_y
    integer :: index,indexNearCell(4),sideth(4),near_kk(4)
    real(prec) :: I_JST,Threshold 
    integer,parameter :: nghn_s = 3
    real(prec),dimension(:,:) :: I_JST_nodes(nsp,nsp),Stencil_JST(2,-nghn_s:nghn_s),volume_nodes(nsp,nsp) 
    real(prec),dimension(:) :: width(nsp)
    !real(prec),dimension(:,:,:) :: Cell_kesi
    real(prec),dimension(nsp,nsp,4) :: varCell
    real(prec),dimension(4,nsp,nsp,4) :: varNearCell
    integer :: indexCellL,indexCellR,ikesiL,ietaL,ikesiR,ietaR,direc
    
    !阙值
    Threshold  = 0.01_prec
    do k = 1,nsp
        do l =1,nsp
            volume_nodes(k,l) = (FPs(k+1)-FPs(k))*(FPs(l+1)-FPs(l)) 
        end do
        width(k) = FPs(k+1)-FPs(k)
    end do  
    
    !!## 单元内每个节点都有一个indicator value
    do i = 1,ncells
        index = i
        varCell(:,:,:) = cellset(index).spvalue_ori(:,:,:)                  !Current cell
        !varCell(:,:,4) = cellset(index).spvalue_ori(:,:,1) * cellset(index).spvalue_ori(:,:,4)                  !Current cell rho p
        indexNearCell  = cellset(i).nearcells                               !Near cells
        do k = 1,4                                                          !Mark side index of the near cells
            varNearCell(k,:,:,:) = cellset(indexNearCell(k)).spvalue_ori(:,:,:)   !Near cells
            !varNearCell(k,:,:,4) = cellset(indexNearCell(k)).spvalue_ori(:,:,1)*cellset(indexNearCell(k)).spvalue_ori(:,:,4)!rho p
            do j = 1,4
                if(cellset(indexNearCell(k)).nearcells(j)==index)then
                    sideth(k) = j                                           !Sideth in near cells
                end if
            end do
        end do  
        
       if(detect_type == ByCell)then
            !!## 所有求解点
            direc = 0
            do k = 1,nsp
                do l = 1,nsp
                    !! First Step. Obtain the stencil of the detection                
                    call get_JST_stencil(k,l,nghn_s,varCell,varNearCell,sideth,Stencil_JST)
                    !! Second Step. Compute Ijst 
                    call compute_JST(direc,nghn_s,Stencil_JST,I_JST_nodes(k,l))                
                end do
            end do
         
            I_JST = sum(volume_nodes * abs(I_JST_nodes))/4.0_prec   
            if(I_JST > Threshold )then
                cellset(i).Beta_line = 1
                cellset(i).Beta = 1
            end if 
        
       elseif(detect_type == ByDim)then
            !kesi direction 
            direc = 1
            do k = 1,nsp
                do l = 1,nsp
                    !! First Step. Obtain the stencil of the detection                
                    call get_JST_stencil(k,l,nghn_s,varCell,varNearCell,sideth,Stencil_JST)
                    !! Second Step. Compute Ijst 
                    call compute_JST(direc,nghn_s,Stencil_JST,I_JST_nodes(k,l))                
                end do
            end do
            do k = 1,nsp
                I_JST = sum(width(:) * abs(I_JST_nodes(k,:)))/2.0_prec   
                if(I_JST > Threshold )then
                    cellset(i).Beta_line(k,1) = 1
                    cellset(i).Beta = 1
                end if                                                
            end do  
            
            !eta direction 
            direc = 2
            do k = 1,nsp
                do l = 1,nsp
                    !! First Step. Obtain the stencil of the detection                
                    call get_JST_stencil(k,l,nghn_s,varCell,varNearCell,sideth,Stencil_JST)
                    !! Second Step. Compute Ijst 
                    call compute_JST(direc,nghn_s,Stencil_JST,I_JST_nodes(k,l))                
                end do
               
            end do
            
            do k = 1,nsp
                 
                I_JST = sum(width(:) * abs(I_JST_nodes(:,k)))/2.0_prec   
                if(I_JST > Threshold )then
                    !write(*,*)varNearCell(3,:,k,1)!varCell(:,k,1)!I_JST_nodes(:,k),
                    cellset(i).Beta_line(k,2) = 1
                    cellset(i).Beta = 1
                end if         
            end do
        end if
    end do  
    !stop
    end subroutine indicator_JST
    
subroutine get_JST_stencil(k,l,nghn_s,varCell,varNearCell,sideth,Stencil_JST)
    use real_precision
    use parameter_setting,only : nsp
    implicit none
    real(prec),dimension(:,:) :: Stencil_JST(2,-nghn_s:nghn_s) 
    real(prec),dimension(nsp,nsp,4) :: varCell
    real(prec),dimension(4,nsp,nsp,4) :: varNearCell
    integer :: sideth(4),k,l,nghn_s,RC_k,RC_l,m
    integer,parameter :: varth = 4  ! p
    real(prec),dimension(:) :: stencil_temp(nsp*3)
    Stencil_JST = 0.0_prec
    stencil_temp = 0.0_prec
    !!kesi
    if(sideth(4) == 1)then
        call get_RC_l(4,1,k,RC_k)
        do m = 1,nsp
            stencil_temp(m) = varNearCell(4,nsp+1-m,RC_k,varth)
        end do   
    elseif(sideth(4) == 2)then
        call get_RC_l(4,2,k,RC_k)
        stencil_temp(1:nsp) = varNearCell(4,RC_k,:,varth)
    elseif(sideth(4) == 3)then
        call get_RC_l(4,3,k,RC_k)
        stencil_temp(1:nsp) = varNearCell(4,:,RC_k,varth)
    elseif(sideth(4) == 4)then
        call get_RC_l(4,4,k,RC_k)
        do m = 1,nsp
            stencil_temp(m) = varNearCell(4,RC_k,nsp+1-m,varth)
        end do   
    end if
    stencil_temp(nsp+1:nsp+nsp) = varCell(k,:,varth)
    if(sideth(2) == 1)then
        call get_RC_l(2,1,k,RC_k)
        stencil_temp(2*nsp+1:3*nsp) = varNearCell(2,:,RC_k,varth)
    elseif(sideth(2) == 2)then
        call get_RC_l(2,2,k,RC_k)
        do m = 1,nsp
            stencil_temp(2*nsp+m) = varNearCell(2,RC_k,nsp+1-m,varth)
        end do 
    elseif(sideth(2) == 3)then
        call get_RC_l(2,3,k,RC_k)
        do m = 1,nsp
            stencil_temp(2*nsp+m) = varNearCell(2,nsp+1-m,RC_k,varth)
        end do 
    elseif(sideth(2) == 4)then
        call get_RC_l(2,4,k,RC_k)
        stencil_temp(2*nsp+1:3*nsp) = varNearCell(2,RC_k,:,varth)
    end if
    Stencil_JST(1,-nghn_s:nghn_s) = stencil_temp(nsp+l-nghn_s:nsp+l+nghn_s)
    !!eta
    if(sideth(1) == 1)then
        call get_RC_l(1,1,l,RC_l)
        do m = 1,nsp
            stencil_temp(m) = varNearCell(1,nsp+1-m,RC_l,varth)
        end do 
        !stencil_temp(1:nsp) = varNearCell(1,:,RC_l,varth)
    elseif(sideth(1) == 2)then
        call get_RC_l(1,2,l,RC_l)
        stencil_temp(1:nsp) = varNearCell(1,RC_l,:,varth)
    elseif(sideth(1) == 3)then
        call get_RC_l(1,3,l,RC_l)
        stencil_temp(1:nsp) = varNearCell(1,:,RC_l,varth)
        !write(*,*)stencil_temp(1:nsp)
    elseif(sideth(1) == 4)then
        call get_RC_l(1,4,l,RC_l)
        do m = 1,nsp
            stencil_temp(m) = varNearCell(1,RC_l,nsp+1-m,varth)
        end do 
        !stencil_temp(1:nsp) = varNearCell(1,RC_l,:,varth)
    end if
    stencil_temp(nsp+1:nsp+nsp) = varCell(:,l,varth)
    if(sideth(3) == 1)then
        call get_RC_l(3,1,l,RC_l)
        stencil_temp(2*nsp+1:3*nsp) = varNearCell(3,:,RC_l,varth)
    elseif(sideth(3) == 2)then
        call get_RC_l(3,2,l,RC_l)
        stencil_temp(2*nsp+1:3*nsp) = varNearCell(3,RC_l,:,varth)
    elseif(sideth(3) == 3)then
        call get_RC_l(3,3,l,RC_l)
        do m = 1,nsp
            stencil_temp(2*nsp+m) = varNearCell(3,nsp+1-m,RC_l,varth)
        end do 
        !stencil_temp(2*nsp+1:3*nsp) = varNearCell(3,:,RC_l,varth)
    elseif(sideth(3) == 4)then
        call get_RC_l(3,4,l,RC_l)
        do m = 1,nsp
            stencil_temp(2*nsp+m) = varNearCell(3,RC_l,nsp+1-m,varth)
        end do 
        !stencil_temp(2*nsp+1:3*nsp) = varNearCell(3,RC_l,:,varth)
    end if
        
    Stencil_JST(2,-nghn_s:nghn_s) = stencil_temp(nsp+k-nghn_s:nsp+k+nghn_s)
    
end subroutine get_JST_stencil
    
subroutine compute_JST(direc,nghn_s,Stencil_JST,I_JST_node)
    use real_precision
    implicit none
    integer :: nghn_s,direc
    real(prec),dimension(:,:) :: Stencil_JST(2,-nghn_s:nghn_s),Stencil_JST_temp(2,1:nghn_s+nghn_s)
    real(prec) :: I_JST_node,min_var,max_var
    
    
    Stencil_JST_temp(:,1:nghn_s) = Stencil_JST(:,-nghn_s:-1)
    Stencil_JST_temp(:,nghn_s+1:2*nghn_s) = Stencil_JST(:,1:nghn_s)
    if(direc == 0)then
        min_var = minVal(Stencil_JST_temp)
        max_var = maxVal(Stencil_JST_temp)        
    elseif(direc == 1)then!kesi
        min_var = minVal(Stencil_JST_temp(1,:))
        max_var = maxVal(Stencil_JST_temp(1,:))
    elseif(direc == 2)then!eta
        min_var = minVal(Stencil_JST_temp(2,:))
        max_var = maxVal(Stencil_JST_temp(2,:)) 
        !write(*,*)min_var,max_var
    end if
    I_JST_node = (min_var - 2.0_prec*Stencil_JST(1,0) + max_var)/(min_var + 2.0_prec*Stencil_JST(1,0) + max_var)
    
end subroutine compute_JST