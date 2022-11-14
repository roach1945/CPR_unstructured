
subroutine indicator_ATV
    !##       Troubled-cell indicator based on the average total variation of the solution.
    !Refs:[1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29.

    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,near_k,ikesi1,ikesi2,ieta1,ieta2,ikesi,ieta
    integer :: sideth(4),near_kk(4),sidethL,sidethR,indexNearCell,indexNearCell1, indexNearCell2, indexNearCell3, indexNearCell4 
    integer :: index,indexCellL,indexCellR,Beta,Beta_line_x,Beta_line_y
    real(prec) :: aveCell,varFluxL,varFluxR,varmodL,varmodR,a1,a2,a3,dh2
    real(prec),dimension(4) :: aveNearCell,TV,ATV,Max_TV,Max_v,Min_v
    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
    real(prec),dimension(nsp,nsp,4) :: varOri,varOriL,varOriR
    real(prec),dimension(4,nsp,nsp,4) :: varOriNearCell
    real(prec) :: deta_L,deta_R,deta_CL,deta_CR,aa,temp,theta,Vtemp_x(4),Vtemp_y(4)
    
    theta = 0.5_prec
    !-------------------------------------------------------
    call solve_ATV(ATV)
    !write(*,*) ATV
    
    !-------------------------------------------------------
    !!本单元 -- 物理域原始变量
    do i = 1,ncells
        Vtemp_x = 0.0_prec
        Vtemp_y = 0.0_prec
        varOri(:,:,:) = cellset(i).spvalue_ori(:,:,:)    
        !邻单元
        do j = 1,4  !4个相邻单元
            indexNearCell  = cellset(i).nearcells(j)!相邻单元的索引
            do k = 1,4
                if(i==cellset(indexNearCell).nearcells(k))then
                    sideth(j) = k   !记录本单元的侧边是相邻单元 的第几侧边
                end if
            end do       
            varOriNearCell(j,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
        end do 
        !write(*,*)i,sideth
        !----x direction --------------------------------------------------------------------
        beta = 0
        do m = 1,1
            !write(*,*)'12'
            do k =1,nsp
                if(sideth(4)==1)then    !x方向，靠近侧边4的第一个点，需要和邻单元相对的点做差。取出该点的当地坐标
                    ikesi = 1;          ieta = k;
                elseif(sideth(4)==2)then
                    ikesi = k;    ieta = nsp
                elseif(sideth(4)==3)then
                    ikesi = nsp;   ieta = nsp+1-k;
                elseif(sideth(4)==4)then
                    ikesi = nsp+1-k;     ieta = 1
                end if  
                !write(*,*)'13',ikesi,ieta
                
                temp =  abs(varOri(k,1,m)-varOriNearCell(4,ikesi,ieta,m))
                Vtemp_x(m) = max(Vtemp_x(m),temp)                
                
                do l = 2,nsp
                    temp =  abs(varOri(k,l,m)-varOri(k,l-1,m))
                    Vtemp_x(m) = max(Vtemp_x(m),temp)
                    if(temp>theta*ATV(m))then    !MV : theta*Max_TV(m);  ATV : theta*ATV(m)
                        cellset(i).Beta_line(k,1) = 1
                        cellset(i).Beta = 1
                    end if
                    !write(*,*)'2',temp,theta*Max_TV(m)
                end do 
                !x方向，靠近侧边2的第一个点，需要和邻单元相对的点做差。取出该点的当地坐标
                if(sideth(2)==1)then    
                    ikesi = 1;          ieta = nsp+1-k;
                elseif(sideth(2)==2)then
                    ikesi = nsp+1-k;    ieta = nsp
                elseif(sideth(2)==3)then
                    ikesi = nsp;   ieta = k;
                elseif(sideth(2)==4)then
                    ikesi = k;     ieta = 1
                end if  
                temp =  abs(varOri(k,nsp,m)-varOriNearCell(2,ikesi,ieta,m))
                Vtemp_x(m) = max(Vtemp_x(m),temp)

            !--------------------------------------------------------------------------------------------------       
    
            !----y direction ----------------------------------------------------------------------------------
                if(sideth(1)==1)then    !x方向，靠近侧边1的第一个点，需要和邻单元相对的点做差。取出该点的当地坐标
                    ikesi = 1;          ieta = nsp+1-k;
                elseif(sideth(1)==2)then
                    ikesi = nsp+1-k;    ieta = nsp
                elseif(sideth(1)==3)then
                    ikesi = nsp;   ieta = k;
                elseif(sideth(1)==4)then
                    ikesi = k;     ieta = 1
                end if  
                temp =  abs(varOri(1,k,m)-varOriNearCell(1,ikesi,ieta,m))
                Vtemp_y(m) = max(Vtemp_y(m),temp)
        
                do l = 2,nsp
                    temp =  abs(varOri(l,k,m)-varOri(l-1,k,m))
                    Vtemp_y(m) = max(Vtemp_y(m),temp)
                end do 
                !y方向，靠近侧边3的第一个点，需要和邻单元相对的点做差。取出该点的当地坐标
                if(sideth(3)==1)then    
                    ikesi = 1;          ieta = k;
                elseif(sideth(3)==2)then
                    ikesi = k;    ieta = nsp
                elseif(sideth(3)==3)then
                    ikesi = nsp;   ieta = nsp+1-k;
                elseif(sideth(3)==4)then
                    ikesi = nsp+1-k;     ieta = 1
                end if  
                temp =  abs(varOri(nsp,k,m)-varOriNearCell(3,ikesi,ieta,m))
                Vtemp_y(m) = max(Vtemp_y(m),temp)
                
                if(Vtemp_x(m)>theta*ATV(m))then
                    cellset(i).Beta_line(k,1) = 1
                    cellset(i).Beta = 1
                    if(detect_type == ByCell )then
                        cellset(i).Beta_line = 1
                        exit   
                    end if 
                end if
                if(Vtemp_y(m)>theta*ATV(m))then
                    cellset(i).Beta_line(k,2) = 1
                    cellset(i).Beta = 1
                    if(detect_type == ByCell )then
                        cellset(i).Beta_line = 1
                        exit   
                    end if
                end if
            end do     
        end do        
    end do

end subroutine indicator_ATV

subroutine solve_ATV(ATV)
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,near_k
    integer :: sidethL,sidethR
    integer :: index,indexCellL,indexCellR
    real(prec) :: temp,ikesi1,ikesi2,ieta1,ieta2
    real(prec),dimension(4) :: aveNearCell,TV,ATV,Max_TV
    real(prec),dimension(nsp,nsp,4) :: varOri,varOriL,varOriR

    
    ! 先计算TV（the total variation of the solution）.相邻求解点之间都需要计算。
    TV = 0.0_prec
    !首先计算一个单元内部求解点的TV
    do i = 1,ncells               
        varOri(:,:,:) = cellset(i).spvalue_ori(:,:,:)!!依次取出单元 -- 物理域原始变量             
        do m = 1,1  !侦测变量
            do k = 1,nsp
                do l = 1,nsp-1
                    temp = abs(varOri(k,l+1,m) - varOri(k,l,m)) !x direction
                    TV(m) = TV(m) + temp
                end do
            end do      
        end do
    end do        
    !然后根据边计算边的相邻单元内紧接边的相对求解点
    do i = 1,nsides        
        indexCellL = sideset(i).nearcells(1)
        indexCellR = sideset(i).nearcells(2)
        if(indexCellL == 0 .OR. indexCellR == 0)cycle
        varOriL(:,:,:) = cellset(indexCellL).spvalue_ori(:,:,:)
        varOriR(:,:,:) = cellset(indexCellR).spvalue_ori(:,:,:)
        do j = 1,4
            if(i == cellset(indexCellL).sides(j)) sidethL = j            !计算属于邻单元的第几侧边
            if(i == cellset(indexCellR).sides(j)) sidethR = j
        end do
                
        do m = 1,1
            do k = 1,nsp
                if((sidethL*sidethR==2).or.(sidethL*sidethR==12).or.(sidethL==sidethR))then!若是满足此条件，则求解点编号相反
                    near_k = nsp+1-k
                else
                    near_k = k
                end if
                if(sidethL==1)then
                    ikesi1 = 1;     ieta1 = k;
                elseif(sidethL==2)then
                    ikesi1 = k;     ieta1 = nsp
                elseif(sidethL==3)then
                    ikesi1 = nsp;   ieta1 = k;
                elseif(sidethL==4)then
                    ikesi1 = k;     ieta1 = 1
                end if
                if(sidethR==1)then
                    ikesi2 = 1;     ieta2 = near_k
                elseif(sidethR==2)then
                    ikesi2 = near_k;ieta2 = nsp
                elseif(sidethR==3)then
                    ikesi2 = nsp;   ieta2 = near_k
                elseif(sidethR==4)then
                    ikesi2 = near_k;ieta2 = 1
                end if
                temp = abs(varOriL(ikesi1,ieta1,m)-varOriR(ikesi2,ieta2,m))
                TV(m) = TV(m) + temp
            end do
        end do
    end do
    ATV = TV/(ncells*nsp*(nsp-1)*2+(nsides-nbdsides)*nsp)
end subroutine solve_ATV