
!!!---!modal energy-----------------------------------------------------------------------------------------------------------------
subroutine indicator_MDH!modal energy
    !Refs:
    !   [3]	Hennemann S, Rueda-Ramírez A M, Hindenlang F J, et al. A provably entropy stable subcell shock capturing approach for high order split form DG for the compressible Euler equations [J]. Journal of Computational Physics, 2021, 426(109935).    use parameter_setting
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,index,Beta,Beta_line_x,Beta_line_y
    real(prec) :: sum_M1,sum_M2,E_Le ,T_Le,a,np,c,onp
    real(prec),dimension(:) :: M_Le(nsp)  

    real(prec),dimension(nsp,nsp,4) :: varOriCell           
    
    !T_Le = 0.005_prec*10**(-1.8*(nsp)**0.25)
    a = 0.5_prec
    c = -1.8_prec
    np = nsp-1
    !onp = 1.0_prec/np
    onp = 0.25_prec
    T_Le = a*10**(c*(np+1)**onp)
    do i  = 1,ncells
        varOriCell(:,:,1) = cellset(i).spvalue_ori(:,:,1)*cellset(i).spvalue_ori(:,:,4)
        
        do m = 1,1      
            do k = 1,nsp
                !x direction :varOriCell(k,:,m)
                do l = 1,nsp
                    M_Le(l) = dot_product(K_La_Le(l,:),varOriCell(k,:,m))                      
                end do

                sum_M1 = 0.0_prec
                sum_M2 = 0.0_prec
                do j = 1,nsp
                    sum_M1 = sum_M1+M_Le(j)**2
                end do
                do j = 1,nsp-1
                    sum_M2 = sum_M2+M_Le(j)**2
                end do
                E_Le = max(M_Le(nsp)**2/sum_M1,M_Le(nsp-1)**2/sum_M2)
                
                if(E_Le > T_Le)then                
                    cellset(i).Beta_line(k,1) = 1
                    cellset(i).Beta = 1     
                end if                     
                
                !y direction :varOriCell(:,k,m)
                do l = 1,nsp
                    M_Le(l) = dot_product(K_La_Le(l,:),varOriCell(:,k,m))                      
                end do

                sum_M1 = 0.0_prec
                sum_M2 = 0.0_prec
                do j = 1,nsp
                    sum_M1 = sum_M1+M_Le(j)**2
                end do
                do j = 1,nsp-1
                    sum_M2 = sum_M2+M_Le(j)**2
                end do
                E_Le = max(M_Le(nsp)**2/sum_M1,M_Le(nsp-1)**2/sum_M2)
                if(E_Le > T_Le)then                
                    cellset(i).Beta_line(k,2) = 1
                    cellset(i).Beta = 1    
                end if           
            end do              
        end do
    end do                                                                                

end subroutine indicator_MDH

!!!---!modal energy-----------------------------------------------------------------------------------------------------------------
subroutine indicator_MDH2!modal energy
    !Refs:
    !   [3]	Hennemann S, Rueda-Ramírez A M, Hindenlang F J, et al. A provably entropy stable subcell shock capturing approach for high order split form DG for the compressible Euler equations [J]. Journal of Computational Physics, 2021, 426(109935).    use parameter_setting
    !   对ME做了一点改进：引进边界值，构造nsp+1阶多项式，在单元界面之间有大的间断时也能捕捉到。边界值的rho*p 采取Roe Average
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    
    integer :: i,j,k,l,m,index,Beta,Beta_line_x,Beta_line_y
    real(prec) :: sum_M1,sum_M2,E_Le ,T_Le,a,np,c,onp,a2,T_Le2
    real(prec),dimension(:) :: M_Le(nsp+2),ruvpL(4),ruvpR(4),RoeAverage_sides(4,nsp)  
    real(prec),dimension(nsp,nsp+2,4) :: varOriCell_kesi,varOriCell_eta
    integer :: indexCellL,indexCellR,indexNearCell,sideth(4),near_kk(4),ikesiL,ietaL,ikesiR,ietaR

    varOriCell_kesi = 0.0_prec
    varOriCell_eta = 0.0_prec
    a = 0.5_prec
    c = -1.8_prec
    np = nsp+1
    onp = 0.25_prec
    T_Le = a*10**(c*(np+1)**onp)

    do i  = 1,ncells       
        do j = 1,4  !4个相邻单元
            indexNearCell  = cellset(i).nearcells(j)!相邻单元的索引
            do k = 1,4
                if(i==cellset(indexNearCell).nearcells(k))then
                    sideth(j) = k   !记录本单元的侧边是相邻单元 的第几侧边
                end if
            end do        
        end do 
        do k = 1,nsp
            !内部点赋值
            varOriCell_kesi(k,2:nsp+1,1) = cellset(i).spvalue_ori(k,:,1)*cellset(i).spvalue_ori(k,:,4)
            varOriCell_eta(k,2:nsp+1,1)  = cellset(i).spvalue_ori(:,k,1)*cellset(i).spvalue_ori(:,k,4)
        end do
        !补充边界点
        do j = 1,4
            do k = 1,nsp
                !邻单元k的编号          
                if(sideth(j)*j==2.or.sideth(j)*j==12.or.sideth(j)==j)then
                    near_kk(j) = nsp+1-k    !判断邻单元是第几条
                else 
                    near_kk(j) = k
                end if
                !不同情况
                if(j==1)then    
                    ikesiL = 1;     ietaL = k;
                elseif(j==2)then
                    ikesiL = k;    ietaL = nsp
                elseif(j==3)then
                    ikesiL = nsp;   ietaL = k;
                elseif(j==4)then
                    ikesiL = k;     ietaL = 1
                end if
                if(sideth(j)==1)then    
                    ikesiR = 1;          ietaR = near_kk(j);
                elseif(sideth(j)==2)then
                    ikesiR = near_kk(j);    ietaR = nsp
                elseif(sideth(j)==3)then
                    ikesiR = nsp;   ietaR = near_kk(j);
                elseif(sideth(j)==4)then
                    ikesiR = near_kk(j);     ietaR = 1
                end if
            
                indexCellL = i
                indexCellR = cellset(i).nearcells(j)
                ruvpL = cellset(indexCellL).spvalue_ori(ikesiL,ietaL,:)
                ruvpR = cellset(indexCellR).spvalue_ori(ikesiR,ietaR,:)
                call RoeAverage(ruvpL,ruvpR,RoeAverage_sides(j,k))     
                !write(*,*)RoeAverage_sides(j,k)
            end do
        end do
        varOriCell_kesi(:,1,1)      = RoeAverage_sides(4,:)
        varOriCell_kesi(:,nsp+2,1)  = RoeAverage_sides(2,:)
        varOriCell_eta(:,1,1)       = RoeAverage_sides(1,:)
        varOriCell_eta(:,nsp+2,1)   = RoeAverage_sides(3,:)
        !write(*,*)varOriCell_eta(:,:,1)
        !stop
        do m = 1,1           
            do k = 1,nsp
                !x direction :varOriCell(k,:,m)
                do l = 1,nsp+2
                    M_Le(l) = dot_product(K_La_Le_2(l,:),varOriCell_kesi(k,:,m))                      
                end do
                !if(i==11)then
                !    write(*,"(2I4 ,7F20.10)")i,k,M_Le
                !    stop
                !end if
                sum_M1 = 0.0_prec
                sum_M2 = 0.0_prec
                do j = 1,nsp+2
                    sum_M1 = sum_M1+M_Le(j)**2
                end do
                do j = 1,nsp+1
                    sum_M2 = sum_M2+M_Le(j)**2
                end do
                E_Le = max(M_Le(nsp+2)**2/sum_M1,M_Le(nsp+1)**2/sum_M2)
                !write(*,"(2I4 ,2F20.10)")i,k,M_Le(nsp+2)**2/sum_M1,M_Le(nsp+1)**2/sum_M2
                if(E_Le > T_Le)then                
                    !write(*,*) i,k,E_Le,T_Le
                    !cellset(i).Beta_line(k,1) = 1
                    !cellset(i).Beta = 1     
                    if(detect_type == ByDim )then
                        cellset(i).Beta_line(k,1)=1
                        cellset(i).Beta = 1
                    elseif(detect_type == ByCell)then
                        cellset(i).Beta_line = 1
                        cellset(i).Beta = 1
                    end if
                end if                     
                
                !y direction :varOriCell(:,k,m)
                do l = 1,nsp+2
                    M_Le(l) = dot_product(K_La_Le_2(l,:),varOriCell_eta(k,:,m))                      
                end do
                !write(*,"(2I4 ,5F20.10)")i,k,M_Le
                sum_M1 = 0.0_prec
                sum_M2 = 0.0_prec
                do j = 1,nsp+2
                    sum_M1 = sum_M1+M_Le(j)**2
                end do
                do j = 1,nsp+1
                    sum_M2 = sum_M2+M_Le(j)**2
                end do
                E_Le = max(M_Le(nsp+2)**2/sum_M1,M_Le(nsp+1)**2/sum_M2)
                if(E_Le > T_Le)then                
                    !write(*,*) i,k,E_Le,T_Le
                    if(detect_type == ByDim )then
                        cellset(i).Beta_line(k,2)=1
                        cellset(i).Beta = 1!在分维侦测中，此句仅为标记整个单元位置
                    elseif(detect_type == ByCell)then
                        cellset(i).Beta_line = 1
                        cellset(i).Beta = 1
                    end if
                end if           
            end do              
        end do
    end do                                                                                

end subroutine indicator_MDH2
subroutine indicator_MDH3!modal energy
    !Refs:
    !   [3]	Hennemann S, Rueda-Ramírez A M, Hindenlang F J, et al. A provably entropy stable subcell shock capturing approach for high order split form DG for the compressible Euler equations [J]. Journal of Computational Physics, 2021, 426(109935).    use parameter_setting
    !   对ME做了一点改进：引进边界值，构造nsp+1阶多项式，在单元界面之间有大的间断时也能捕捉到。边界值的rho*p 采取Roe Average
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,index,Beta,Beta_line_x,Beta_line_y
    real(prec) :: sum_M1,sum_M2,E_Le ,T_Le,T_down,a,np,c,onp,a2,T_Le2,T_h,T_m,T_s,T_s_l,T_s_h
    real(prec),dimension(:) :: M_Le(nsp+2),ruvpL(4),ruvpR(4),RoeAverage_sides(4,nsp)  
    real(prec),dimension(nsp,nsp+2,4) :: varOriCell_kesi,varOriCell_eta
    integer :: indexCellL,indexCellR,indexNearCell,sideth(4),near_kk(4),ikesiL,ietaL,ikesiR,ietaR
    integer :: Beta_line_down,Beta_down

    varOriCell_kesi = 0.0_prec
    varOriCell_eta = 0.0_prec
 
    a = 0.5_prec
    c = -1.8_prec
    np = nsp+1
    onp = 0.25_prec
    T_Le = log(a*10**(c*(np+1)**onp))!-4.0_prec  双马赫分维计算时采取这个阈值
    !d_smooth = 4.0_prec  !
    
    !权值过渡
    T_s = T_Le - 2.0_prec*d_smooth
    T_s_l =  T_s - d_smooth
    T_s_h =  T_s + d_smooth
    
    !上下阙值
    T_down = T_Le - 4.0_prec
    !write(*,*) T_s_h,T_Le,-3.0_prec*log(5.0_prec)
    !stop
    !T_Le = T_m
    do i  = 1,ncells       
        do j = 1,4  !4个相邻单元
            indexNearCell  = cellset(i).nearcells(j)!相邻单元的索引
            do k = 1,4
                if(i==cellset(indexNearCell).nearcells(k))then
                    sideth(j) = k   !记录本单元的侧边是相邻单元 的第几侧边
                end if
            end do        
        end do 
        do k = 1,nsp
            !内部点赋值
            varOriCell_kesi(k,2:nsp+1,1) = cellset(i).spvalue_ori(k,:,1)*cellset(i).spvalue_ori(k,:,4)
            varOriCell_eta(k,2:nsp+1,1)  = cellset(i).spvalue_ori(:,k,1)*cellset(i).spvalue_ori(:,k,4)
        end do
        !补充边界点
        do j = 1,4
            do k = 1,nsp
                !邻单元k的编号          
                if(sideth(j)*j==2.or.sideth(j)*j==12.or.sideth(j)==j)then
                    near_kk(j) = nsp+1-k    !判断邻单元是第几条
                else 
                    near_kk(j) = k
                end if
                !不同情况
                if(j==1)then    
                    ikesiL = 1;     ietaL = k;
                elseif(j==2)then
                    ikesiL = k;    ietaL = nsp
                elseif(j==3)then
                    ikesiL = nsp;   ietaL = k;
                elseif(j==4)then
                    ikesiL = k;     ietaL = 1
                end if
                if(sideth(j)==1)then    
                    ikesiR = 1;          ietaR = near_kk(j);
                elseif(sideth(j)==2)then
                    ikesiR = near_kk(j);    ietaR = nsp
                elseif(sideth(j)==3)then
                    ikesiR = nsp;   ietaR = near_kk(j);
                elseif(sideth(j)==4)then
                    ikesiR = near_kk(j);     ietaR = 1
                end if
            
                indexCellL = i
                indexCellR = cellset(i).nearcells(j)
                ruvpL = cellset(indexCellL).spvalue_ori(ikesiL,ietaL,:)
                ruvpR = cellset(indexCellR).spvalue_ori(ikesiR,ietaR,:)
                call RoeAverage(ruvpL,ruvpR,RoeAverage_sides(j,k))     
                !write(*,*)RoeAverage_sides(j,k)
            end do
        end do
        varOriCell_kesi(:,1,1)      = RoeAverage_sides(4,:)
        varOriCell_kesi(:,nsp+2,1)  = RoeAverage_sides(2,:)
        varOriCell_eta(:,1,1)       = RoeAverage_sides(1,:)
        varOriCell_eta(:,nsp+2,1)   = RoeAverage_sides(3,:)
        !write(*,*)varOriCell_eta(:,:,1)
        !stop
        do m = 1,1           
            do k = 1,nsp
        !! ## x direction :varOriCell(k,:,m)
                do l = 1,nsp+2
                    M_Le(l) = dot_product(K_La_Le_2(l,:),varOriCell_kesi(k,:,m))                      
                end do

                sum_M1 = 0.0_prec
                sum_M2 = 0.0_prec
                do j = 1,nsp+2
                    sum_M1 = sum_M1+M_Le(j)**2
                end do
                do j = 1,nsp+1
                    sum_M2 = sum_M2+M_Le(j)**2
                end do
                E_Le = max(M_Le(nsp+2)**2/sum_M1,M_Le(nsp+1)**2/sum_M2)
                !write(*,"(2I4 ,2F20.10)")i,k,M_Le(nsp+2)**2/sum_M1,M_Le(nsp+1)**2/sum_M2
                E_Le = log(E_Le)
                if(E_Le > T_Le)then                    
                    if(detect_type == ByDim )then
                        cellset(i).Beta_line(k,1)=1
                        cellset(i).Beta = 1
                    elseif(detect_type == ByCell)then
                        cellset(i).Beta_line = 1
                        cellset(i).Beta = 1
                    end if
                end if         
                !! smoother
                if(E_Le < T_s_l)then
                    cellset(i).Smooth_line_x(k,:) = 0.0_prec
                elseif(E_Le .GE. T_s_l .and. E_Le .LE. T_s_h)then
                    cellset(i).Smooth_line_x(k,:) = 0.5_prec*(1.0_prec+sin(pi*(E_Le - T_s )/d_smooth))
                elseif(E_Le .GT. T_s_h .and. E_Le .LE. T_Le)then                    
                    cellset(i).Smooth_line_x(k,:) = 1.0_prec
                    !write(*,*)E_Le/T_h,log(E_Le/T_h),min(1.0_prec,E_Le/T_h)   
                elseif(E_Le .GT. T_Le)then                    
                    !间断。不起作用
                    cellset(i).Smooth_line_x(k,:) = 0.0_prec
                end if
                !write(*,*)i,k,cellset(i).Smooth_line_x(k,1)
                !! up down
                
                if(switch_updown == open_updown)then
                    if(E_Le > T_down)then                    
                        Beta_line_down = 1
                        Beta_down = 1
                    else
                        Beta_line_down = 0
                        Beta_down = 0
                    end if 
                    if(cellset(i).Beta_old == 0 .and. cellset(i).Beta == 1)then
                        cellset(i).Beta_line(k,1)=1
                        cellset(i).Beta = 1
                    elseif(cellset(i).Beta_old == 1 .and. Beta_down == 0)then
                        cellset(i).Beta_line(k,1)=0
                        cellset(i).Beta = 0
                    elseif(cellset(i).Beta_old == 1 .and. Beta_down == 1)then
                        cellset(i).Beta_line(k,1)=1
                        cellset(i).Beta = 1
                    endif          
                    
                end if
                
                
    ! ! ##   y direction :varOriCell(:,k,m)
                do l = 1,nsp+2
                    M_Le(l) = dot_product(K_La_Le_2(l,:),varOriCell_eta(k,:,m))                      
                end do
                !write(*,"(2I4 ,5F20.10)")i,k,M_Le
                sum_M1 = 0.0_prec
                sum_M2 = 0.0_prec
                do j = 1,nsp+2
                    sum_M1 = sum_M1+M_Le(j)**2
                end do
                do j = 1,nsp+1
                    sum_M2 = sum_M2+M_Le(j)**2
                end do
                E_Le = max(M_Le(nsp+2)**2/sum_M1,M_Le(nsp+1)**2/sum_M2)
                E_Le = log(E_Le)
                if(E_Le > T_Le)then                
                    !write(*,*) i,k,E_Le,T_Le
                    if(detect_type == ByDim )then
                        cellset(i).Beta_line(k,2)=1
                        cellset(i).Beta = 1!在分维侦测中，此句仅为标记整个单元位置
                    elseif(detect_type == ByCell)then
                        cellset(i).Beta_line = 1
                        cellset(i).Beta = 1
                    end if
                end if 
                
                if(E_Le < T_s_l)then
                    cellset(i).Smooth_line_y(k,:) = 0.0_prec
                elseif(E_Le .GE. T_s_l .and. E_Le .LE. T_s_h)then
                    cellset(i).Smooth_line_y(k,:) = 0.5_prec*(1.0_prec+sin(pi*(E_Le - T_s )/(2.0_prec*d_smooth)))
                elseif(E_Le .GT. T_s_h .and. E_Le .LE. T_Le)then                    
                    cellset(i).Smooth_line_y(k,:) = 1.0_prec        
                elseif(E_Le .GT. T_Le)then                    
                    !间断。不起作用
                    cellset(i).Smooth_line_y(k,:) = 0.0_prec
                end if
                
                if(switch_updown == open_updown)then
                    !只考虑了整个单元侦测，如分维需要beta_line_old
                    if(E_Le > T_down)then                    
                        Beta_line_down = 1
                        Beta_down = 1
                    else
                        Beta_line_down = 0
                        Beta_down = 0
                    end if 
                    if(cellset(i).Beta_old == 0 .and. cellset(i).Beta == 1)then
                        cellset(i).Beta_line(k,2)=1
                        cellset(i).Beta = 1
                    elseif(cellset(i).Beta_old == 1 .and. Beta_down == 0)then
                        cellset(i).Beta_line(k,2)=0
                        cellset(i).Beta = 0
                    elseif(cellset(i).Beta_old == 1 .and. Beta_down == 1)then
                        cellset(i).Beta_line(k,2)=1
                        cellset(i).Beta = 1
                    endif                       
                end if
            end do              
        end do
    end do                                                                                

end subroutine indicator_MDH3
    
subroutine RoeAverage(ruvpL,ruvpR,rhoPRoe)
    use real_precision
    use global_var
    implicit none
    real(prec) :: ruvpL(4),ruvpR(4),rhoL,uL,vL,pL,EL,rhoR,uR,vR,pR,ER,rhoRoe,uRoe,vRoe,pRoe,hRoe,U2_Roe,cRoe,orhoL,orhoR,rhoPRoe
    real(prec) :: srL,srR,srLR,osrLR,hL,hR
        
    !Left and right variables
    rhoL = ruvpL(1);        uL = ruvpL(2);        vL = ruvpL(3);        pL = ruvpL(4);
    rhoR = ruvpR(1);        uR = ruvpR(2);        vR = ruvpR(3);        pR = ruvpR(4);
    !write(*,"('L',4F20.10)")ruvpL
    !write(*,"('R',4F20.10)")ruvpR
        
    orhoL = 1.0_prec/rhoL
    orhoR = 1.0_prec/rhoR
        
    EL = pL/(gamma - 1.0_prec) + 0.5_prec * rhoL*(uL**2 + vL**2)
    ER = pR/(gamma - 1.0_prec) + 0.5_prec * rhoR*(uR**2 + vR**2)
        
    !// Left and right pressures
    !pL = (gamma - 1.0_prec) * (EL - 0.5_prec * rhoL*(uL**2 + vL**2))
    !pR = (gamma - 1.0_prec) * (ER - 0.5_prec * rhoR*(uR**2 + vR**2))
        
    !// Left and right enthalpy
    hL = (EL + pL) * orhoL
    hR = (ER + pR) * orhoR

    !// Square root of rhoL and rhoR.
    srL  = sqrt(rhoL)
    srR  = sqrt(rhoR)
    srLR = srL + srR
    osrLR = 1.0_prec/srLR
        
    !// Velocity, enthalpy and sound speed Roe averages (equation 11.60).
    rhoRoe = sqrt(rhoL*rhoR)
    uRoe   = (srL * uL + srR * uR) * osrLR
    vRoe   = (srL * vL + srR * vR) * osrLR
    hRoe   = (srL * hL + srR * hR) * osrLR
    U2_Roe   = (uRoe * uRoe + vRoe * vRoe)
    !cRoe   = sqrt((gamma - 1.0)*(hRoe - 0.5 * URoe))
    !ocRoe  = 1.0/cRoe      
    pRoe = (hRoe*rhoRoe - 0.5_prec*rhoRoe*U2_Roe)*(gamma - 1.0_prec)/gamma
        
    rhoPRoe = rhoRoe*pRoe
    !write(*,*)rhoPRoe,rhoRoe,pRoe
    !stop
end subroutine RoeAverage
    