subroutine indicator_TVB
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none

    integer :: i,j,k,l,Beta
    integer :: index,indexNearCell(4),sideth(4),near_iL,near_iR
    real(prec),dimension(3,nsp,4) :: varOri         !(cell index, sp num, var num)
    real(prec) :: dh
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
        !write(*,*) indexNearCell
        !write(*,*) sideth
        do j = 1,nsp                                
            !Traverse all rows         4 -- L, 2  -- R
            if((sideth(4)==3).or.(sideth(4)==4))then
                near_iL = nsp+1-j
            else
                near_iL = j
            end if
            
            if((sideth(2)==1).or.(sideth(2)==2))then
                near_iR = nsp+1-j
            else
                near_iR = j
            end if

            !write(*,*) '000',near_iL,near_iR
            varOri(1,:,:) = cellset(indexNearCell(4)).spvalue_ori(near_iL,:,:)
            varOri(2,:,:) = cellset(index).spvalue_ori(j,:,:)
            varOri(3,:,:) = cellset(indexNearCell(2)).spvalue_ori(near_iR,:,:)
            
            !!!----
            dh = (2.0_prec*cellset(index).MJacobi(j,1,1))
            !write(*,*)dh
            !dh = 2.0_prec
            call indicator_TVB_fun(varOri,dh,Beta)
            cellset(index).Beta_line(j,1) = Beta    
            !write(*,*)i,j,Beta
            if(Beta == 1)then
                if(detect_type == ByDim )then
                    !cellset(index).Beta_line(j,1)=1
                    cellset(index).Beta = 1
                elseif(detect_type == ByCell)then
                    cellset(index).Beta_line = 1
                    cellset(index).Beta = 1                        
                    exit    
                end if       
            end if
            !
            !!Traverse all Column s         1 -- L, 3  -- R
            if((sideth(1)==1).or.(sideth(1)==2))then
                near_iL = nsp+1-j
            else
                near_iL = j
            end if
            
            if((sideth(3)==3).or.(sideth(3)==4))then
                near_iR = nsp+1-j
            else
                near_iR = j
            end if
            
            !write(*,*) '000',near_iL,near_iR
            varOri(1,:,:) = cellset(indexNearCell(1)).spvalue_ori(:,near_iL,:)
            varOri(2,:,:) = cellset(index).spvalue_ori(:,j,:)
            varOri(3,:,:) = cellset(indexNearCell(3)).spvalue_ori(:,near_iR,:)
            
            !!!----
            dh = (2.0_prec*cellset(index).MJacobi(1,j,4))
            !write(*,*)dh
            !dh = 2.0_prec
            call indicator_TVB_fun(varOri,dh,Beta)
            cellset(index).Beta_line(j,2) = Beta       
            !write(*,*)i,j,Beta
            if(Beta == 1)then
                if(detect_type == ByDim )then
                    !cellset(index).Beta_line(j,2)=1
                    cellset(index).Beta = 1
                elseif(detect_type == ByCell)then
                    cellset(index).Beta_line = 1
                    cellset(index).Beta = 1                        
                    exit    
                end if       
            end if
        end do           
    end do  
    
end subroutine indicator_TVB

subroutine indicator_TVB_fun(var,dh,Beta)
    use real_precision
    use parameter_setting
    use global_var
    implicit none
    real(prec),external :: Gauss_integral_SPs,Gauss_double_integral,LaI_nPs,m_func,TVB_minmod
    real(prec) :: var(3,nsp,4)        !(cell index, sp num, var num)
    real(prec) :: dh
    real(prec) :: ave,aveL,aveR,faceL,faceR,deta_L,deta_R,deta_CL,deta_CR,varmodL,varmodR
    integer :: m,Beta
    Beta = 0
    !var(:,:,1) = var(:,:,1)*var(:,:,4) !rho*p
    do m =1,1       !r u v p
        ave  = Gauss_integral_SPs(var(2,:,m))/2.0_prec
        aveL = Gauss_integral_SPs(var(1,:,m))/2.0_prec
        aveR = Gauss_integral_SPs(var(3,:,m))/2.0_prec
        !write(*,*)ave,aveL,aveR
        faceL = LaI_nPs(kesi_l,SPs,var(2,:,m),nsp)
        faceR = LaI_nPs(kesi_r,SPs,var(2,:,m),nsp)
        
        deta_L = ave - aveL
        deta_R = aveR - ave          
        deta_CL = ave - faceL
        deta_CR = faceR-ave

        varmodL = TVB_minmod(deta_CL,deta_R,deta_L,dh)      
        varmodR = TVB_minmod(deta_CR,deta_R,deta_L,dh)
            
        if(abs(varmodL - deta_CL)>1.0e-15 .or. abs(varmodR - deta_CR) >1.0e-15)then
            !write(*,*)deta_CL,deta_CR
           Beta = 1
        end if     
    end do
    
end subroutine indicator_TVB_fun


function TVB_minmod(a1,a2,a3,dh)
    use real_precision
    use parameter_setting,only: TVB_M
    implicit none
    real(prec),external :: m_func
    real(prec) :: TVB_minmod,a1,a2,a3,dh,dh2
    dh2 = dh**2
    
    if(abs(a1) <= TVB_M*dh2)then
        TVB_minmod = a1      
    elseif(abs(a1) > TVB_M*dh2)then
        TVB_minmod = m_func(a1,a2,a3)
        !write(*,*) a1,dh2
    endif

end function TVB_minmod

function m_func(a1,a2,a3)
    !m(a1,a2,a3,...,ak)
    use real_precision 
    implicit none
    real(prec),external::min_abs
    real(prec)::a1,a2,a3,m_func
    
    !相同符号
    if((a1*a2>0.0_prec).and.(a1*a3>0.0_prec).and.(a1*a2*a3>0.0_prec))then
        m_func = min_abs(a1,a2,a3)
    elseif((a1*a2>0.0_prec).and.(a1*a3>0.0_prec).and.(a1*a2*a3<0.0_prec))then
        m_func = -1.0_prec*min_abs(a1,a2,a3)
    else
        m_func = 0.0_prec
    end if

    !if(sign(1.0_prec,a1) == sign(1.0_prec,a2) .and. sign(1.0_prec,a1) == sign(1.0_prec,a3)  .and. sign(1.0_prec,a1)>0 )then
    !    m_func = min_abs(a1,a2,a3)
    !elseif(sign(1.0_prec,a1) == sign(1.0_prec,a2) .and. sign(1.0_prec,a1) == sign(1.0_prec,a3)  .and. sign(1.0_prec,a1) < 0 )then
    !    m_func = -1.0_prec*min_abs(a1,a2,a3)
    !else
    !    m_func = 0.0_prec
    !end if
    return
    
end function

function min_abs(a,b,c)
    !取绝对值最小值
    use real_precision
    implicit none
    real(prec)::a,b,c,min_abs
    if(abs(a)<=abs(b).and.abs(a)<=abs(c))then
        min_abs=abs(a)
    elseif(abs(b)<=abs(a).and.abs(b)<=abs(c))then
        min_abs=abs(b)
    else
        min_abs=abs(c)
    end if
    return
end function
!!!--------------------------------------------------------------------------------------------------------------------
