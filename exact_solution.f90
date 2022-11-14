subroutine exact_solution(Time)
    use global_var
    use type_module
    use parameter_setting,only:nsp
    implicit none
    integer :: i,j,k,l
    real(prec) :: Time
    !write(*,*) Time
    
    select case(case_comp)
    case(equEntropy_case)
        do i =1,ncells           
            do j = 1,nsp      
                do k = 1,nsp
                    call equEntropy_init(cellset(i).sp_coor(j,k,:),Time,cellset(i).spvalue_ori_exa(j,k,1:4))
                    !write(*,"(4F10.5)")cellset(i).sp_coor(j,k,:)
                end do
            end do
        end do    
    case(SodShockTube_case)
        call shockTube_exa(Time)
    case(LaxShockTube_case)
        call shockTube_exa(Time) 
    case(SinWave_2D_case)
        do i =1,ncells            
            do j = 1,nsp      
                do k = 1,nsp
                    call Euler_1D_order_init(cellset(i).sp_coor(j,k,:),Time,cellset(i).spvalue_ori_exa(j,k,1:4))
                end do
            end do
        end do 
    case default
        write(*,*)'no exact_solution'
        !stop        
    end select
    
end subroutine exact_solution



subroutine shockTube_exa(endT)
!!---Sod 激波管准确解---------------------------------------
    use parameter_setting
    use type_module
    use global_var
    implicit none
    real(prec),external::P_c_star,f_p_star
    real(prec)::error,delt_p_star,F1,F2,dF_dp_star,p_star_new,delt_u
    real(prec)::c_L,c_R,f_L,f_R,A_L,A_R,x
    real(prec)::c_a,c_L_star,c_R_star,c_star,temp,r_star_L,r_star_R,z_L,z_R,z_L_head,z_L_tail,z_R_head,z_R_tail
    real(prec)::u_star,p_star,eps,endT,xc0
    integer :: i,j,k,l
    real(prec) :: rho_L,u_L,v_L,p_L  !初始时刻，分隔面左侧参数
    real(prec) :: rho_R,u_R,v_R,p_R  !初始时刻，分隔面右侧参数
    real(prec),dimension(:,:,:) :: oriExact(nsp,nsp,4)
    if(case_comp == SodShockTube_case)then
        rho_L = 1.0_prec
        u_L = 0.0_prec
        v_L = 0.0_prec
        p_L = 1.0_prec   !初始时刻，分隔面左侧参数,sod激波管
        rho_R = 0.125_prec
        u_R = 0.0_prec
        v_R = 0.0_prec
        p_R = 0.1_prec !初始时刻，分隔面右侧参数,sod激波管
    else if(case_comp == LaxShockTube_case)then
        rho_L = 0.445_prec
        u_L = 0.698_prec
        v_L = 0.0_prec
        p_L = 3.528_prec   !初始时刻，分隔面左侧参数,lax激波管
        rho_R = 0.5_prec
        u_R = 0.0_prec
        v_R = 0.0_prec
        p_R = 0.571_prec !初始时刻，分隔面右侧参数,lax激波管  
    end if
    !xc0 = 0.5_prec*(xl+xr)
    eps = 1.0e-6
    delt_p_star = 1.0e-4
    delt_u = u_L-u_R
    p_star = min(p_L,p_R)
    error=1.0_prec
    do while (error>eps)
        !牛顿迭代求解中心区P*，用以判断两侧是激波/膨胀波
        F1 = P_c_star(p_star+delt_p_star,p_L,rho_L,p_R,rho_R)
        F2 = P_c_star(p_star,p_L,rho_L,p_R,rho_R)
        dF_dp_star = (F1-F2)/delt_p_star                  !导数
        p_star_new = p_star-(F2-delt_u)/dF_dp_star        !推进
        error = abs(p_star_new-p_star)                    !误差
        p_star = p_star_new
    end do

    !---------------------------------------
    !方程中的重复公式
    !---------------------------------------
    c_L = sqrt(gamma*p_L/rho_L) !左侧c
    c_R = sqrt(gamma*p_R/rho_R) !右侧c
    f_L = f_p_star(p_star,p_L,rho_L)
    f_R = f_p_star(p_star,p_R,rho_R)
    u_star = (u_L+u_R-f_L+f_R)/2.0_prec
    A_L = sqrt((gamma+1.0_prec)/(2.0_prec*gamma)*p_star/p_L+(gamma-1.0_prec)/(2.0_prec*gamma))
    A_R = sqrt((gamma+1.0_prec)/(2.0_prec*gamma)*p_star/p_R+(gamma-1.0_prec)/(2.0_prec*gamma))

    !-----------------------------------------
    !两侧不同情况
    !-----------------------------------------
    do i =1,ncells
        
        do k = 1,nsp
            do l = 1,nsp
                select case(dire_shock)
                case(direct_x)
                    temp = cellset(i).sp_coor(k,l,1) - (xl+xr)*0.5_prec     !x方向流动
                case(direct_y)
                    temp = cellset(i).sp_coor(k,l,2) - (yl+yr)*0.5_prec    !y方向流动
                end select
                if (p_star>p_L .and. p_star>p_R)then              !两侧均为激波      
                    z_L = u_L-A_L*c_L                             !左行激波速度
                    r_star_L = rho_L*c_L*A_L/(u_star-u_L+A_L*c_L)   !左行激波后密度
                    z_R = u_R+A_R*c_R                             !右行激波速度
                    r_star_R = rho_R*c_R*A_R/(u_R-u_star+A_R*c_R)   !右行激波后密度        
        
                    if (temp<(z_L*endT))then                      !左行激波波前参数
                        oriExact(k,l,1) = rho_L
                        oriExact(k,l,2) = u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_L
                    else if (temp<(u_star*endT))then                !左行激波波后参数
                        oriExact(k,l,1) = r_star_L
                        oriExact(k,l,2) = u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_star
                    else if (temp<(z_R*endT))then                   !右行激波波后参数
                        oriExact(k,l,1) = r_star_R
                        oriExact(k,l,2) = u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_star                    
                    else                              !右行激波波前参数
                        oriExact(k,l,1) = rho_R
                        oriExact(k,l,2) = u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_R                     
                    end if
                else if (p_star<p_L .and. p_star<p_R .and. p_star<=0)then     !两侧均为膨胀波，但中心区是真空
                    z_L_head=u_L-c_L                               !左行膨胀波波头速度
                    z_L_tail=u_star-c_star                   !左行膨胀波波尾速度
                    z_R_head=u_R+c_R                               !右行膨胀波波头速度
                    z_R_tail=u_star+c_star                   !右行膨胀波波尾速度

                    if (temp<z_L_head*endT)then               !左行膨胀波波前参数
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp <z_L_tail*endT)then              !左行膨胀波波内参数
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(u_L-temp/endT)+2.0_prec/(gamma+1.0_prec)*c_L
                        oriExact(k,l,2)=temp/endT+c_a
                        oriExact(k,l,4)=p_L*(c_a/c_L)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                    else if (temp<z_R_tail*endT) then               !膨胀波波后参数，由于是真空，所以没有接触间断
                        oriExact(k,l,1)=0.0_prec
                        oriExact(k,l,2)=0.0_prec
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=0.0_prec
                    else if (temp<z_R_head*endT)then               !右行膨胀波波内参数
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*((temp/endT-u_R)+2.0_prec/(gamma+1.0_prec)*c_R)
                        oriExact(k,l,2)=temp/endT-c_a
                        oriExact(k,l,4)=p_R*(c_a/c_R)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                    else                            !右行膨胀波波前参数
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
    
                else if (p_star<p_L .and. p_star<p_R .and. p_star>0)then         !两侧均为膨胀波，但中心区不是真空
                    c_L_star=c_L+(gamma-1.0_prec)/2.0_prec*(u_L-u_star)           !左行膨胀波波后速度
                    z_L_head=u_L-c_L                                !左行膨胀波波头速度
                    z_L_tail=u_star-c_L_star                        !左行膨胀波波尾速度
                    c_R_star=c_R+(gamma-1.0_prec)/2.0_prec*(u_R-u_star)           !右行膨胀波波后速度
                    z_R_head=u_R+c_R                                !右行膨胀波波头速度
                    z_R_tail=u_star+c_R_star                        !右行膨胀波波尾速度

                    if (temp<z_L_head*endT)then                        !左行膨胀波波前参数
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp <z_L_tail*endT)then                     !左行膨胀波波内参数
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(u_L-temp/endT)+2.0_prec/(gamma+1.0_prec)*c_L   !计算膨胀波内部声速
                        oriExact(k,l,2)=temp/endT+c_a
                        oriExact(k,l,4)=p_L*(c_a/c_L)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                    else if (temp<u_star*endT)then                        !左行膨胀波波后参数
                        oriExact(k,l,1)=r_star_L
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_tail*endT)then                     !右行膨胀波波后参数，由于不是真空，所以存在接触间断
                        oriExact(k,l,1)=r_star_R
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_head*endT)then                      !右行膨胀波波内参数
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(temp/endT-u_R)+2.0_prec/(gamma+1.0_prec)*c_R
                        oriExact(k,l,2)=temp/endT-c_a
                        oriExact(k,l,4)=p_R*(c_a/c_R)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1) = gamma*oriExact(k,l,4)/(c_a**2)
                        oriExact(k,l,3) = 0.0_prec
                    else                                   !右行膨胀波波前参数
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
                else if (p_star<p_L .and. p_star>p_R)then		        !左侧为膨胀波，右侧为激波
                    c_L_star=c_L+(gamma-1.0_prec)/2.0_prec*(u_L-u_star)       !左行膨胀波波后速度
                    r_star_L=rho_L*(p_star/p_L)**(1.0_prec/gamma)        !左行膨胀波波后密度
                    z_L_head=u_L-c_L                            !左行膨胀波波头速度
                    z_L_tail=u_star-c_L_star                    !左行膨胀波波尾速度
                    z_R=u_R+A_R*c_R                             !右行激波速度
                    r_star_R=rho_R*c_R*A_R/(u_R-u_star+A_R*c_R)   !右行激波波后密度

                    if (temp<z_L_head*endT)then                    !左行膨胀波波前参数
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp<z_L_tail*endT)then                  !左行膨胀波波内参数
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(u_L-temp/endT)+2.0_prec/(gamma+1.0_prec)*c_L
                        oriExact(k,l,2)=temp/endT+c_a
                        oriExact(k,l,4)=p_L*(c_a/c_L)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                        oriExact(k,l,3) = 0.0_prec
                    else if (temp<u_star*endT)then                    !左行膨胀波波后参数
                        oriExact(k,l,1)=r_star_L
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R*endT)then                       !右行激波波后参数
                        oriExact(k,l,1)=r_star_R
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else                               !右行激波波前参数
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
                else                                             !左侧为激波，右侧为膨胀波
                    z_L=u_L-A_L*c_L                              !左行激波速度
                    r_star_L=rho_L*c_L*A_L/(u_star-u_L+A_L*c_L)      !左行激波波后密度
                    c_R_star=c_R+(gamma-1.0_prec)/2.0_prec*(u_R-u_star)        !右行膨胀波波后声速
                    z_R_head=u_R+c_R                             !右行膨胀波波头速度
                    z_R_tail=u_star+c_R_star                     !右行膨胀波波尾速度
                    r_star_R=gamma*p_star/(c_R_star*c_R_star)    !右行膨胀波波后密度
                    if (temp<z_L*endT)then                          !左行激波波前参数
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp<u_star*endT)then                     !左行激波波后参数
                        oriExact(k,l,1)=r_star_L
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_tail*endT)then                   !右行膨胀波波后参数
                        oriExact(k,l,1)=r_star_R
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_head*endT)then                   !右行膨胀波波内参数
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(temp/endT-u_R)+2.0_prec/(gamma+1.0_prec)*c_R
                        oriExact(k,l,2)=temp/endT-c_a
                        oriExact(k,l,4)=p_R*(c_a/c_R)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                        oriExact(k,l,3) = 0.0_prec
                    else                                !右行膨胀波波前参数
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
                end if
            end do
        end do
        select case(dire_shock)
        case(direct_x)
            cellset(i).spvalue_ori_exa = oriExact
        case(direct_y)
            cellset(i).spvalue_ori_exa(:,:,1) = oriExact(:,:,1) 
            cellset(i).spvalue_ori_exa(:,:,2) = oriExact(:,:,3) 
            cellset(i).spvalue_ori_exa(:,:,3) = oriExact(:,:,2) 
            cellset(i).spvalue_ori_exa(:,:,4) = oriExact(:,:,4) 
        end select
        
    end do
end subroutine shockTube_exa

function f_p_star(p_star,P_1,r_1)
    ! P*
    use real_precision
    use global_var,only:gamma
    implicit none
    real(prec)::c_1,p_star,P_1,r_1,f_p_star
    c_1=sqrt(gamma*P_1/r_1)
    if (p_star<0.0_prec)then
        write(*,*) 'Error pressure <0!!!'
    else if (p_star<1.0e-10 .or. p_star==0.0_prec)then
        f_p_star=-2.0_prec*c_1/(gamma-1.0_prec)
    else if(p_star>P_1)then
        f_p_star=(p_star-P_1)/(r_1*c_1*sqrt((gamma+1.0_prec)/(2.0_prec*gamma)*p_star/P_1+(gamma-1.0_prec)/(2.0_prec*gamma)))
    else
        f_p_star=2.0_prec*c_1/(gamma-1.0_prec)*((p_star/P_1)**((gamma-1.0_prec)/(2.0_prec*gamma))-1.0_prec)
    end if

    return
end function
function P_c_star(p_star,P_1,r_1,P_2,r_2)
    !中心区P*
    use real_precision
    real(prec),external::f_p_star
    real(prec)::P_c_star,p_star,P_1,r_1,P_2,r_2
    P_c_star = f_p_star(p_star,P_1,r_1)+f_p_star(p_star,P_2,r_2)
    return
end function