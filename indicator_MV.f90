subroutine indicator_MV
    !Refs:[1]	Li G, Qiu J. Hybrid weighted essentially non-oscillatory schemes with different indicators [J]. Journal of Computational Physics, 2010, 229(21): 8105-29.
    !   实际上，对ATV有了小的改动。原文中是采取Average Total Variation,对问题依赖严重。在此将平均总变差改为最大变差，
    !   但是如此，也会有问题，那就是光滑问题必有问题单元，且在强激波存在的情况下，可能对弱间断侦测不灵敏。
    !   留待验证方案：1、采取全局最大变量值*theta作为判断依据；2、采取（全局最大变量值-最小值）*theta作为判断依据。
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,m,near_k,ikesi1,ikesi2,ieta1,ieta2,ikesi,ieta
    integer :: sideth(4),near_kk(4),sidethL,sidethR,indexNearCell
    integer :: index,indexCellL,indexCellR,Beta
    real(prec),dimension(4) :: Max_TV,Max_v,Min_v
    real(prec),dimension(nsp,nsp,4) :: varOriCell,varOriCellL,varOriCellR
    real(prec),dimension(4,nsp,nsp,4) :: varOriCellNearCell
    real(prec) :: temp,theta,Vtemp_x(4),Vtemp_y(4)
    
    theta = 0.5_prec
    ! 先计算TV（the total variation of the solution）.相邻求解点之间都需要计算。
    Max_TV = 0.0_prec
    Max_v = 0.0_prec
    Min_v = 0.0_prec
    !首先计算一个单元内部求解点的TV
    do i = 1,ncells               
        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)!!依次取出单元 -- 物理域原始变量             
        do m = 1,1  !侦测变量
            Max_v = MaxVal(varOriCell(:,:,m))
            Min_v = MinVal(varOriCell(:,:,m))
            do k = 1,nsp
                do l = 1,nsp-1
                    temp = abs(varOriCell(k,l+1,m) - varOriCell(k,l,m)) !x direction
                    Max_TV(m) = max(Max_TV(m),temp)
                    temp = abs(varOriCell(l+1,k,m) - varOriCell(l,k,m)) !y
                    Max_TV(m) = max(Max_TV(m),temp)
                end do
            end do      
        end do
    end do        
    
    !然后根据边计算边的相邻单元内紧接边的相对求解点
    do i = 1,nsides        
        indexCellL = sideset(i).nearcells(1)
        indexCellR = sideset(i).nearcells(2)
        if(indexCellL == 0 .OR. indexCellR == 0)cycle
        varOriCellL(:,:,:) = cellset(indexCellL).spvalue_ori(:,:,:)
        varOriCellR(:,:,:) = cellset(indexCellR).spvalue_ori(:,:,:)
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
                temp = abs(varOriCellL(ikesi1,ieta1,m)-varOriCellR(ikesi2,ieta2,m))               
                Max_TV(m) = max(Max_TV(m),temp)
            end do
        end do
    end do

    !-------------------------------------------------------
    !!本单元 -- 物理域原始变量
    do i = 1,ncells
        Vtemp_x = 0.0_prec
        Vtemp_y = 0.0_prec
        varOriCell(:,:,:) = cellset(i).spvalue_ori(:,:,:)    
        !邻单元
        do j = 1,4  !4个相邻单元
            indexNearCell  = cellset(i).nearcells(j)!相邻单元的索引
            do k = 1,4
                if(i==cellset(indexNearCell).nearcells(k))then
                    sideth(j) = k   !记录本单元的侧边是相邻单元 的第几侧边
                end if
            end do 
            varOriCellNearCell(j,:,:,:) = cellset(indexNearCell).spvalue_ori(:,:,:)        
        end do 
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
   
                temp =  abs(varOriCell(k,1,m)-varOriCellNearCell(4,ikesi,ieta,m))
                Vtemp_x(m) = max(Vtemp_x(m),temp)                
                
                do l = 2,nsp
                    temp =  abs(varOriCell(k,l,m)-varOriCell(k,l-1,m))
                    Vtemp_x(m) = max(Vtemp_x(m),temp)
                    if(temp>theta*Max_TV(m))then    !MV : theta*Max_TV(m);  
                        cellset(i).Beta_line(k,1) = 1
                        cellset(i).Beta = 1
                    end if
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
                temp =  abs(varOriCell(k,nsp,m)-varOriCellNearCell(2,ikesi,ieta,m))
                Vtemp_x(m) = max(Vtemp_x(m),temp)

            !--------------------------------------------------------------------------------------------------       
    
            !----y direction ----------------------------------------------------------------------------------
                if(sideth(1)==1)then    !y方向，靠近侧边1的第一个点，需要和邻单元相对的点做差。取出该点的当地坐标
                    ikesi = 1;          ieta = nsp+1-k;
                elseif(sideth(1)==2)then
                    ikesi = nsp+1-k;    ieta = nsp
                elseif(sideth(1)==3)then
                    ikesi = nsp;   ieta = k;
                elseif(sideth(1)==4)then
                    ikesi = k;     ieta = 1
                end if  
                temp =  abs(varOriCell(1,k,m)-varOriCellNearCell(1,ikesi,ieta,m))
                Vtemp_y(m) = max(Vtemp_y(m),temp)
        
                do l = 2,nsp
                    temp =  abs(varOriCell(l,k,m)-varOriCell(l-1,k,m))
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
                temp =  abs(varOriCell(nsp,k,m)-varOriCellNearCell(3,ikesi,ieta,m))
                Vtemp_y(m) = max(Vtemp_y(m),temp)
                
                !标记问题单元
                if(Vtemp_x(m)>theta*Max_TV(m))then
                    cellset(i).Beta_line(k,1) = 1
                    cellset(i).Beta = 1
                    if(detect_type == ByCell )then
                        cellset(i).Beta_line = 1
                        exit   
                    end if 
                end if
                if(Vtemp_y(m)>theta*Max_TV(m))then
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
    

end subroutine indicator_MV