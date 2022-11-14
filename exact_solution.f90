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
!!---Sod ������׼ȷ��---------------------------------------
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
    real(prec) :: rho_L,u_L,v_L,p_L  !��ʼʱ�̣��ָ���������
    real(prec) :: rho_R,u_R,v_R,p_R  !��ʼʱ�̣��ָ����Ҳ����
    real(prec),dimension(:,:,:) :: oriExact(nsp,nsp,4)
    if(case_comp == SodShockTube_case)then
        rho_L = 1.0_prec
        u_L = 0.0_prec
        v_L = 0.0_prec
        p_L = 1.0_prec   !��ʼʱ�̣��ָ���������,sod������
        rho_R = 0.125_prec
        u_R = 0.0_prec
        v_R = 0.0_prec
        p_R = 0.1_prec !��ʼʱ�̣��ָ����Ҳ����,sod������
    else if(case_comp == LaxShockTube_case)then
        rho_L = 0.445_prec
        u_L = 0.698_prec
        v_L = 0.0_prec
        p_L = 3.528_prec   !��ʼʱ�̣��ָ���������,lax������
        rho_R = 0.5_prec
        u_R = 0.0_prec
        v_R = 0.0_prec
        p_R = 0.571_prec !��ʼʱ�̣��ָ����Ҳ����,lax������  
    end if
    !xc0 = 0.5_prec*(xl+xr)
    eps = 1.0e-6
    delt_p_star = 1.0e-4
    delt_u = u_L-u_R
    p_star = min(p_L,p_R)
    error=1.0_prec
    do while (error>eps)
        !ţ�ٵ������������P*�������ж������Ǽ���/���Ͳ�
        F1 = P_c_star(p_star+delt_p_star,p_L,rho_L,p_R,rho_R)
        F2 = P_c_star(p_star,p_L,rho_L,p_R,rho_R)
        dF_dp_star = (F1-F2)/delt_p_star                  !����
        p_star_new = p_star-(F2-delt_u)/dF_dp_star        !�ƽ�
        error = abs(p_star_new-p_star)                    !���
        p_star = p_star_new
    end do

    !---------------------------------------
    !�����е��ظ���ʽ
    !---------------------------------------
    c_L = sqrt(gamma*p_L/rho_L) !���c
    c_R = sqrt(gamma*p_R/rho_R) !�Ҳ�c
    f_L = f_p_star(p_star,p_L,rho_L)
    f_R = f_p_star(p_star,p_R,rho_R)
    u_star = (u_L+u_R-f_L+f_R)/2.0_prec
    A_L = sqrt((gamma+1.0_prec)/(2.0_prec*gamma)*p_star/p_L+(gamma-1.0_prec)/(2.0_prec*gamma))
    A_R = sqrt((gamma+1.0_prec)/(2.0_prec*gamma)*p_star/p_R+(gamma-1.0_prec)/(2.0_prec*gamma))

    !-----------------------------------------
    !���಻ͬ���
    !-----------------------------------------
    do i =1,ncells
        
        do k = 1,nsp
            do l = 1,nsp
                select case(dire_shock)
                case(direct_x)
                    temp = cellset(i).sp_coor(k,l,1) - (xl+xr)*0.5_prec     !x��������
                case(direct_y)
                    temp = cellset(i).sp_coor(k,l,2) - (yl+yr)*0.5_prec    !y��������
                end select
                if (p_star>p_L .and. p_star>p_R)then              !�����Ϊ����      
                    z_L = u_L-A_L*c_L                             !���м����ٶ�
                    r_star_L = rho_L*c_L*A_L/(u_star-u_L+A_L*c_L)   !���м������ܶ�
                    z_R = u_R+A_R*c_R                             !���м����ٶ�
                    r_star_R = rho_R*c_R*A_R/(u_R-u_star+A_R*c_R)   !���м������ܶ�        
        
                    if (temp<(z_L*endT))then                      !���м�����ǰ����
                        oriExact(k,l,1) = rho_L
                        oriExact(k,l,2) = u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_L
                    else if (temp<(u_star*endT))then                !���м����������
                        oriExact(k,l,1) = r_star_L
                        oriExact(k,l,2) = u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_star
                    else if (temp<(z_R*endT))then                   !���м����������
                        oriExact(k,l,1) = r_star_R
                        oriExact(k,l,2) = u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_star                    
                    else                              !���м�����ǰ����
                        oriExact(k,l,1) = rho_R
                        oriExact(k,l,2) = u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4) = p_R                     
                    end if
                else if (p_star<p_L .and. p_star<p_R .and. p_star<=0)then     !�����Ϊ���Ͳ����������������
                    z_L_head=u_L-c_L                               !�������Ͳ���ͷ�ٶ�
                    z_L_tail=u_star-c_star                   !�������Ͳ���β�ٶ�
                    z_R_head=u_R+c_R                               !�������Ͳ���ͷ�ٶ�
                    z_R_tail=u_star+c_star                   !�������Ͳ���β�ٶ�

                    if (temp<z_L_head*endT)then               !�������Ͳ���ǰ����
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp <z_L_tail*endT)then              !�������Ͳ����ڲ���
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(u_L-temp/endT)+2.0_prec/(gamma+1.0_prec)*c_L
                        oriExact(k,l,2)=temp/endT+c_a
                        oriExact(k,l,4)=p_L*(c_a/c_L)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                    else if (temp<z_R_tail*endT) then               !���Ͳ������������������գ�����û�нӴ����
                        oriExact(k,l,1)=0.0_prec
                        oriExact(k,l,2)=0.0_prec
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=0.0_prec
                    else if (temp<z_R_head*endT)then               !�������Ͳ����ڲ���
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*((temp/endT-u_R)+2.0_prec/(gamma+1.0_prec)*c_R)
                        oriExact(k,l,2)=temp/endT-c_a
                        oriExact(k,l,4)=p_R*(c_a/c_R)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                    else                            !�������Ͳ���ǰ����
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
    
                else if (p_star<p_L .and. p_star<p_R .and. p_star>0)then         !�����Ϊ���Ͳ������������������
                    c_L_star=c_L+(gamma-1.0_prec)/2.0_prec*(u_L-u_star)           !�������Ͳ������ٶ�
                    z_L_head=u_L-c_L                                !�������Ͳ���ͷ�ٶ�
                    z_L_tail=u_star-c_L_star                        !�������Ͳ���β�ٶ�
                    c_R_star=c_R+(gamma-1.0_prec)/2.0_prec*(u_R-u_star)           !�������Ͳ������ٶ�
                    z_R_head=u_R+c_R                                !�������Ͳ���ͷ�ٶ�
                    z_R_tail=u_star+c_R_star                        !�������Ͳ���β�ٶ�

                    if (temp<z_L_head*endT)then                        !�������Ͳ���ǰ����
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp <z_L_tail*endT)then                     !�������Ͳ����ڲ���
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(u_L-temp/endT)+2.0_prec/(gamma+1.0_prec)*c_L   !�������Ͳ��ڲ�����
                        oriExact(k,l,2)=temp/endT+c_a
                        oriExact(k,l,4)=p_L*(c_a/c_L)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                    else if (temp<u_star*endT)then                        !�������Ͳ��������
                        oriExact(k,l,1)=r_star_L
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_tail*endT)then                     !�������Ͳ�������������ڲ�����գ����Դ��ڽӴ����
                        oriExact(k,l,1)=r_star_R
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_head*endT)then                      !�������Ͳ����ڲ���
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(temp/endT-u_R)+2.0_prec/(gamma+1.0_prec)*c_R
                        oriExact(k,l,2)=temp/endT-c_a
                        oriExact(k,l,4)=p_R*(c_a/c_R)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1) = gamma*oriExact(k,l,4)/(c_a**2)
                        oriExact(k,l,3) = 0.0_prec
                    else                                   !�������Ͳ���ǰ����
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
                else if (p_star<p_L .and. p_star>p_R)then		        !���Ϊ���Ͳ����Ҳ�Ϊ����
                    c_L_star=c_L+(gamma-1.0_prec)/2.0_prec*(u_L-u_star)       !�������Ͳ������ٶ�
                    r_star_L=rho_L*(p_star/p_L)**(1.0_prec/gamma)        !�������Ͳ������ܶ�
                    z_L_head=u_L-c_L                            !�������Ͳ���ͷ�ٶ�
                    z_L_tail=u_star-c_L_star                    !�������Ͳ���β�ٶ�
                    z_R=u_R+A_R*c_R                             !���м����ٶ�
                    r_star_R=rho_R*c_R*A_R/(u_R-u_star+A_R*c_R)   !���м��������ܶ�

                    if (temp<z_L_head*endT)then                    !�������Ͳ���ǰ����
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp<z_L_tail*endT)then                  !�������Ͳ����ڲ���
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(u_L-temp/endT)+2.0_prec/(gamma+1.0_prec)*c_L
                        oriExact(k,l,2)=temp/endT+c_a
                        oriExact(k,l,4)=p_L*(c_a/c_L)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                        oriExact(k,l,3) = 0.0_prec
                    else if (temp<u_star*endT)then                    !�������Ͳ��������
                        oriExact(k,l,1)=r_star_L
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R*endT)then                       !���м����������
                        oriExact(k,l,1)=r_star_R
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else                               !���м�����ǰ����
                        oriExact(k,l,1)=rho_R
                        oriExact(k,l,2)=u_R
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_R
                    end if
                else                                             !���Ϊ�������Ҳ�Ϊ���Ͳ�
                    z_L=u_L-A_L*c_L                              !���м����ٶ�
                    r_star_L=rho_L*c_L*A_L/(u_star-u_L+A_L*c_L)      !���м��������ܶ�
                    c_R_star=c_R+(gamma-1.0_prec)/2.0_prec*(u_R-u_star)        !�������Ͳ���������
                    z_R_head=u_R+c_R                             !�������Ͳ���ͷ�ٶ�
                    z_R_tail=u_star+c_R_star                     !�������Ͳ���β�ٶ�
                    r_star_R=gamma*p_star/(c_R_star*c_R_star)    !�������Ͳ������ܶ�
                    if (temp<z_L*endT)then                          !���м�����ǰ����
                        oriExact(k,l,1)=rho_L
                        oriExact(k,l,2)=u_L
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_L
                    else if (temp<u_star*endT)then                     !���м����������
                        oriExact(k,l,1)=r_star_L
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_tail*endT)then                   !�������Ͳ��������
                        oriExact(k,l,1)=r_star_R
                        oriExact(k,l,2)=u_star
                        oriExact(k,l,3) = 0.0_prec
                        oriExact(k,l,4)=p_star
                    else if (temp<z_R_head*endT)then                   !�������Ͳ����ڲ���
                        c_a=(gamma-1.0_prec)/(gamma+1.0_prec)*(temp/endT-u_R)+2.0_prec/(gamma+1.0_prec)*c_R
                        oriExact(k,l,2)=temp/endT-c_a
                        oriExact(k,l,4)=p_R*(c_a/c_R)**(2.0_prec*gamma/(gamma-1.0_prec))
                        oriExact(k,l,1)=gamma*oriExact(k,l,4)/(c_a*c_a)
                        oriExact(k,l,3) = 0.0_prec
                    else                                !�������Ͳ���ǰ����
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
    !������P*
    use real_precision
    real(prec),external::f_p_star
    real(prec)::P_c_star,p_star,P_1,r_1,P_2,r_2
    P_c_star = f_p_star(p_star,P_1,r_1)+f_p_star(p_star,P_2,r_2)
    return
end function