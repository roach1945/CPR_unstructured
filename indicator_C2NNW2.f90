subroutine indicator_C2NNW2
    
    ! 好像没有实现完全，废弃
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none

    integer :: i,j,k,l,Beta
    integer :: index,indexNearCell(4),sideth(4),near_jL,near_jR
    real(prec),dimension(3,4) :: varOri         !(sps, var num)
    real(prec) :: dis(4)
    varOri = 0.0_prec                               !initialization
    do i = 1,ncells                                 !Traverse all cells
        index = i
        indexNearCell  = cellset(i).nearcells
        do k = 1,4                                  !mark side index of the near cells
            do j = 1,4
                if(cellset(indexNearCell(k)).nearcells(j)==index)then
                    sideth(k) = j
                end if
            end do
        end do     

        loop1:do j = 1,nsp                                
            !Traverse all rows         4 -- L, 2  -- R
            if((sideth(4)==3).or.(sideth(4)==4))then
                near_jL = nsp+1-j
            else
                near_jL = j
            end if
            !外部点--左
            if(sideth(4) == 1)then
                varOri(1,:) = cellset(indexNearCell(4)).spvalue_ori(1,near_jL,:)
            elseif(sideth(4) == 2 )then
                varOri(1,:) = cellset(indexNearCell(4)).spvalue_ori(near_jL,nsp,:)
            elseif(sideth(4) == 3)then
                varOri(1,:) = cellset(indexNearCell(4)).spvalue_ori(nsp,near_jL,:)
            elseif(sideth(4) == 4 )then
                varOri(1,:) = cellset(indexNearCell(4)).spvalue_ori(near_jL,1,:)
            end if

            varOri(2:3,:) = cellset(index).spvalue_ori(j,1:2,:)
            !write(*,*)i,j,indexNearCell(4)
            !write(*,*)varOri(:,1)
            call indicator_C2NNW2_fun(dis,varOri,Beta)
            if(Beta == 1)then
                if(detect_type == ByDim )then
                    !cellset(index).Beta_line(j,1)=1
                    cellset(index).Beta = 1
                elseif(detect_type == ByCell)then
                    cellset(index).Beta_line = 1
                    cellset(index).Beta = 1                        
                    exit loop1    
                end if       
            end if
            !内部点
            do l = 2,nsp-1
                !!原始变量
                dis(1) = dis_sp_fp(2*(l-1))
                dis(2) = dis_sp_fp(2*(l-1)+1)  
                dis(3) = dis_sp_fp(2*l)
                dis(4) = dis_sp_fp(2*l+1)
                !write(*,*) dis
                varOri = cellset(i).spvalue_ori(j,l-1:l+1,:)
                call indicator_C2NNW2_fun(dis,varOri,Beta)
                if(Beta == 1)then
                    exit
                end if
            end do
            !write(*,*)i,Beta
            cellset(index).Beta_line(j,1)=1
            if(Beta == 1)then          
                if(detect_type == ByDim )then
                    !cellset(index).Beta_line(j,1)=1
                    cellset(index).Beta = 1
                elseif(detect_type == ByCell)then
                    cellset(index).Beta_line = 1
                    cellset(index).Beta = 1                        
                    exit loop1     
                end if       
            end if      
            
            if((sideth(2)==1).or.(sideth(2)==2))then
                near_jR = nsp+1-j
            else
                near_jR = j
            end if
            !外部点--右
            if(sideth(2) == 1)then
                varOri(3,:) = cellset(indexNearCell(2)).spvalue_ori(1,near_jL,:)
            elseif(sideth(2) == 2 )then
                varOri(3,:) = cellset(indexNearCell(2)).spvalue_ori(near_jL,nsp,:)
            elseif(sideth(2) == 3)then
                varOri(3,:) = cellset(indexNearCell(2)).spvalue_ori(nsp,near_jL,:)
            elseif(sideth(2) == 4 )then
                varOri(3,:) = cellset(indexNearCell(2)).spvalue_ori(near_jL,1,:)
            end if

            varOri(1:2,:) = cellset(index).spvalue_ori(j,nsp-1:nsp,:)
            call indicator_C2NNW2_fun(dis,varOri,Beta)
            if(Beta == 1)then
                if(detect_type == ByDim )then
                    !cellset(index).Beta_line(j,1)=1
                    cellset(index).Beta = 1
                elseif(detect_type == ByCell)then
                    cellset(index).Beta_line = 1
                    cellset(index).Beta = 1                        
                    exit loop1    
                end if       
            end if
        end do loop1                     
        !write(*,*) i,Beta
    end do  
    
end subroutine indicator_C2NNW2
subroutine indicator_C2NNW2_fun(dis,u,Beta)
    use global_var
    use parameter_setting
    use type_module
    implicit none
    real(prec),dimension(:,:) :: u(3,4)
    real(prec),dimension(2) :: coor_C,coor_P1,coor_P2,coor_Pfp
    real(prec),dimension(4) :: uA1,uB1,uA2,uB2,du1,du2,du,u_cell_L(4),u_cell_R(4)
    real(prec) :: d1,d2,d3,d4,dd1,dd2,dd3,dd4,dis(4)
    real(prec) :: cc,w1,w2,w3,w4,w5,w6,Fai,varMax,varMin,limA,limB  
    real(prec),external :: lim
    integer :: index,nearCellIndex,Cell_sideth,m,l,j,Beta
    real(prec) :: lim_temp,VarC,Var,temp
    !write(*,*) 'u',u
    !通量点与解点之间的距离
    d1 = dis(1)
    d2 = dis(2)
    d3 = dis(3)
    d4 = dis(4) 
    !write(*,*) d1,d2,d3,d4
    !反距离权
    cc = 1.0_prec!常数1，方便写代码做的替换
    dd1 = cc/d1
    dd2 = cc/d2    
    dd3 = cc/d3
    dd4 = cc/d4
    
    w1 = dd1/(dd1+dd2)
    w2 = dd2/(dd1+dd2)
    w3 = dd3/(dd3+dd4)
    w4 = dd4/(dd3+dd4)
    
    !通量点处值
    uA1 = w1*u(1,:)+w2*u(2,:)
    uB1 = w3*u(2,:)+w4*u(3,:)
    
    !单元导数
    w5 = dd2/(dd2+dd3)
    w6 = dd3/(dd2+dd3)
    du1 = (u(2,:)-uA1)/d2
    du2 = (uB1-u(2,:))/d3
    du = w5*du1 + w6*du2

    !重新计算uA,uB
    uA2 = u(2,:)-du*d2
    uB2 = u(2,:)+du*d3
    
    Beta = 0
    do m = 1,1      
        varMax = maxVal(u(:,m))
        varMin = minVal(u(:,m))
        call indicator_C2NNW2_fun2(u(2,m),uA2(m),varMax,varMin,Beta)
        if(Beta == 1)then
            exit
        end if
        call indicator_C2NNW2_fun2(u(2,m),uB2(m),varMax,varMin,Beta)
        if(Beta == 1)then
            exit
        end if
    end do
    !write(*,*) Beta
end subroutine indicator_C2NNW2_fun

subroutine indicator_C2NNW2_fun2(Var,VarC,varMax,varMin,Beta)
    use real_precision
    implicit none
    real(prec) :: Var,VarC,varMax,varMin,temp,delta
    integer :: Beta
    temp = 0.0_prec
    delta = 1.0e-6
    if(VarC > Var)then
        temp = (varMax-Var)/(VarC-Var+delta)
    elseif(VarC < Var)then
        temp = (varMin-Var)/(VarC-Var+delta)
    endif
    !write(*,*) temp
    if(temp > 1.0_prec)then
        Beta = 1
    else
        Beta = 0
    endif
    
end subroutine indicator_C2NNW2_fun2
