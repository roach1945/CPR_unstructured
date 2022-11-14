!!!--------------------------------------------------------------------------------------------------------------------
!!!---KXRCF-----------------------------------------------------------------------------------------------------------------
subroutine indicator_KXRCF
    !Refs:
    !   [2]	Krivodonova L, Xin J, Remacle J F, et al. Shock detection and limiting with discontinuous Galerkin methods for hyperbolic conservation laws [J]. Applied Numerical Mathematics, 2004, 48(3-4): 323-38.
    use real_precision
    use type_module
    use parameter_setting
    use global_var
    implicit none
    integer :: i,j,k,l,Beta
    integer :: index,indexNearCell(4),sideth(4),near_jL,near_jR
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
                near_jL = nsp+1-j
            else
                near_jL = j
            end if
            
            if((sideth(2)==1).or.(sideth(2)==2))then
                near_jR = nsp+1-j
            else
                near_jR = j
            end if

            !----kesi direction-------------------
            varOri(1,:,:) = cellset(indexNearCell(4)).spvalue_ori(near_jL,:,:)
            varOri(2,:,:) = cellset(index).spvalue_ori(j,:,:)
            varOri(3,:,:) = cellset(indexNearCell(2)).spvalue_ori(near_jR,:,:)
            
            !!!----
            dh = (2.0_prec*cellset(index).MJacobi(j,1,1))
            call indicator_KXRCF_fun(varOri,sideth(4),dh,Beta)
            !write(*,*) Beta
            cellset(index).Beta_line(j,1) = Beta    
            !write(*,*)i,j,Beta
            if(Beta == 1)then
                if(detect_type == ByDim )then       
                    cellset(index).Beta = 1
                elseif(detect_type == ByCell)then
                    cellset(index).Beta_line = 1
                    cellset(index).Beta = 1                        
                    exit    
                end if       
            end if
            
            if(solver_dim == dim_1D)cycle       !1维情况 eta方向不做侦测
            
            !!Traverse all Column s         1 -- L, 3  -- R
            if((sideth(1)==1).or.(sideth(1)==2))then
                near_jL = nsp+1-j
            else
                near_jL = j
            end if
            
            if((sideth(3)==3).or.(sideth(3)==4))then
                near_jR = nsp+1-j
            else
                near_jR = j
            end if

            !----eta direction----------------------------------
            varOri(1,:,:) = cellset(indexNearCell(1)).spvalue_ori(:,near_jL,:)
            varOri(2,:,:) = cellset(index).spvalue_ori(:,j,:)
            varOri(3,:,:) = cellset(indexNearCell(3)).spvalue_ori(:,near_jR,:)
            
            !!!----
            dh = (2.0_prec*cellset(index).MJacobi(1,j,4))
            !write(*,*)dh
            !dh = 2.0_prec
            call indicator_KXRCF_fun(varOri,sideth(1),dh,Beta)
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

end subroutine indicator_KXRCF

!!!--------------------------------------------------------------------------------------------------------------------
subroutine indicator_KXRCF_fun(var,sideth,dh,Beta)
    use real_precision
    use parameter_setting
    use global_var
    implicit none
    real(prec),external :: Gauss_integral_SPs,LaI_nPs
    real(prec) :: var(3,nsp,4)        !(cell index, sp num, var num)
    real(prec) :: dh,temp
    real(prec) :: ave,currentFace,neighborFace,detaFace,varmodL,varmodR
    integer :: m,Beta,sideth
    Beta = 0
    do m = 1,1       !r u v p
        ave  = Gauss_integral_SPs(var(2,:,m))/2.0_prec
        !write(*,*) dh
        !flow into
        currentFace = LaI_nPs(kesi_l,SPs,var(2,:,m),nsp)
        if(sideth == 4 .or. sideth == 1)then
            neighborFace = LaI_nPs(kesi_l,SPs,var(1,:,m),nsp)
        elseif(sideth == 2 .or. sideth == 3)then
            neighborFace = LaI_nPs(kesi_r,SPs,var(1,:,m),nsp)
        end if
        
        detaFace = currentFace - neighborFace
        temp = abs(detaFace)/(dh**(nsp*0.5_prec) * abs(ave))  !p+1 = nsp p is the degree of polynomial
        !write(*,*)temp
        if(temp > 1.0_prec)then
           Beta = 1
        end if     
    end do
    
end subroutine indicator_KXRCF_fun