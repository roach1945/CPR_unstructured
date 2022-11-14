subroutine set_IC 
     
    !-----------------------------------------------------------------------------
    !
    !   方法：初始流场设置
    !   作者：gqShi 
    !
    !-----------------------------------------------------------------------------
 
    use global_var
    use parameter_setting
    use type_module
    implicit none
    
    integer i,j,k,l
    
    select case(case_comp)
    case(equEntropy_case)
        do i =1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call equEntropy_init(cellset(i).sp_coor(j,k,:),0.0,cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do  
    case(SodShockTube_case)
        do i =1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call shockWave_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do  
    case(LaxShockTube_case)
        do i =1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call shockWave_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do     
    case(ShuOsher_case)
        do i =1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call ShuOsher_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(Riemann2D_case)
        do i =1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call Riemann2D_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(DoubleMach_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call DoubleMach_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(VortexShock_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call VortexShock_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(CompositeVortexShock_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call CompositeVortexShock_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(steadyShock_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call steadyShock_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(HyperCylinder_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call HyperCylinder_case_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(WingFlow_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call WingFlow_case_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(SinWave_2D_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call Euler_1D_order_init(cellset(i).sp_coor(j,k,:),0.0,cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(Doublerarefaction_1D_case)     !一维稀疏波问题
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    !call Euler_1D_order_init(cellset(i).sp_coor(j,k,:),0.0,cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case(test_case)
        do i = 1,ncells+nbdsides
            allocate(cellset(i).spvalue_ori(nsp,nsp,4))
            do j = 1,nsp      
                do k = 1,nsp
                    call test_init(cellset(i).sp_coor(j,k,:),cellset(i).spvalue_ori(j,k,1:4))
                end do
            end do
        end do
    case default
        write(*,*)'non. Occur in subroutine set_IC '
        stop
    end select
      
end subroutine set_IC 

subroutine equEntropy_init(coor,endT,slt)

    use global_var
    use parameter_setting
    implicit none
    
    real(prec) :: xx0,yy0,r0,u0,v0,p0,xc0,yc0,conser(1:4),slt(1:4),endT
    real(prec) :: xcore,ycore,EKC,RB2,RB,temexp,TE,tt,disx,disy,xc1,yc1
    real(prec),dimension(:) :: coor(2)
    
    xx0 = coor(1)
    yy0 = coor(2)
    EKC = 5.0_prec
    !初始时刻，向右均匀场
    r0 = 1.0_prec
    u0 = 1.0_prec
    v0 = 0.0_prec
    p0 = 1.0_prec
    xlong = xr - xl
    ylong = yr - yl
    !涡的位置
    xc0 = 0.5_prec*(xr+xl)
    yc0 = 0.5_prec*(yr+yl)

    disx=(u0*endT/xlong-int(u0*endT/xlong))*xlong
    disy=(v0*endT/ylong-int(v0*endT/ylong))*ylong

    !真实的涡心位置
    if(disx<=0.5_prec*xlong-xc0)then
        xc1=xc0 + disx
    else
        xc1=xc0 + disx -xlong
    end if

    if(disy<=0.5_prec*xlong-yc0)then
        yc1=yc0 + disy
    else
        yc1=yc0 + disy -xlong
    end if

    !计算所用的涡心位置
    if(xc1>=0.0_prec)then
        if (xx0>=(xc1-0.5_prec*xlong)) then
            xcore = xc1
        else
            xcore = xc1 - xlong
        end if
    else
        if (xx0<=(xc1+0.5_prec*xlong)) then
            xcore = xc1
        else
            xcore = xc1 + xlong
        end if
    end if


    if(yc1>=0.0_prec)then
        if (yy0>=(yc1-0.5_prec*ylong)) then
            ycore = yc1
        else
            ycore = yc1 - ylong
        end if
    else
        if (yy0<=(yc1+0.5_prec*ylong)) then
            ycore = yc1
        else
            ycore = yc1 + ylong
        end if
    end if


    RB2 =  (xx0-xcore)*(xx0-xcore) + (yy0-ycore)*(yy0-ycore)
    RB  = sqrt(RB2)
    temexp = exp(0.5_prec*(1.0_prec-RB2))

    TE = -(gamma-1.0_prec)*EKC*EKC*temexp*temexp/(8.0_prec*gamma*pi*pi)
    tt = 1.0_prec + TE
    slt(1)  = tt**(1.0_prec/(gamma-1.0_prec))     
    slt(2)  = u0 -EKC*temexp*(yy0-ycore)/(2.0_prec*pi)
    slt(3)  = v0 +EKC*temexp*(xx0-xcore)/(2.0_prec*pi)
    slt(4) = slt(1)**gamma
    
    !slt(1) = r0
    !slt(2) = u0
    !slt(3) = v0
    !slt(4) = p0
end subroutine equEntropy_init

subroutine shockWave_init(coor,ruvp)

    !激波管初始场
    use global_var
    use parameter_setting
    implicit none
    
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
    real(prec) :: rho_L,u_L,v_L,p_L  
    real(prec) :: rho_R,u_R,v_R,p_R  
    

    if(case_comp == SodShockTube_case)then
        rho_L = 1.0_prec
        u_L = 0.0_prec
        v_L = 0.0_prec
        p_L = 1.0_prec      
        
        rho_R = 0.125_prec
        u_R = 0.0_prec
        v_R = 0.0_prec
        p_R = 0.1_prec      
    else if(case_comp == LaxShockTube_case)then
        rho_L = 0.445_prec
        u_L = 0.698_prec
        v_L = 0.0_prec
        p_L = 3.528_prec    
        
        rho_R = 0.5_prec
        u_R = 0.0_prec
        v_R = 0.0_prec
        p_R = 0.571_prec    
              
    end if
    
    x = coor(1)
    y = coor(2)
    select case(dire_shock)
    case(direct_x)
        !x 方向
        if(x .LT. 0.5_prec*(xl+xr))then
            ruvp(1)=rho_L
            ruvp(2)=u_L
            ruvp(3)=v_L
            ruvp(4)=p_L
        elseif(x .GE. 0.5_prec*(xl+xr))then
            ruvp(1)=rho_R
            ruvp(2)=u_R
            ruvp(3)=v_R
            ruvp(4)=p_R    
        endif
    case(direct_y)
        !y 方向
        if(y .LT. 0.5_prec*(yl+yr))then
            ruvp(1)=rho_L
            ruvp(2)=v_L
            ruvp(3)=u_L
            ruvp(4)=p_L
        elseif(y .GE. 0.5_prec*(yl+yr))then
            ruvp(1)=rho_R
            ruvp(2)=v_R
            ruvp(3)=u_R
            ruvp(4)=p_R    
        endif
    end select
    
    !ruvp(1)=1.0_prec
    !ruvp(2)=1.0_prec
    !ruvp(3)=0.0_prec
    !ruvp(4)=1.0_prec
end subroutine shockWave_init

subroutine ShuOsher_init(coor,ruvp)

    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
    real(prec) :: rho_L,u_L,v_L,p_L
    real(prec) :: rho_R,u_R,v_R,p_R  
    
    
    rho_L = 3.857143_prec
    u_L = 2.629369_prec
    v_L = 0.0_prec
    p_L = 10.33333_prec   
        
    
    u_R = 0.0_prec
    v_R = 0.0_prec
    p_R = 1.0_prec 

    
    x = coor(1)
    y = coor(2)
    select case(dire_shock)
    case(direct_x)
        !x 方向
        if(x .LT. xl+0.1_prec*(xr-xl))then
            ruvp(1)=rho_L
            ruvp(2)=u_L
            ruvp(3)=v_L
            ruvp(4)=p_L
        elseif(x .GE. xl+0.1_prec*(xr-xl))then
            ruvp(1)=1.0_prec + 0.2_prec*sin(50.0_prec/xlong*x)
            ruvp(2)=u_R
            ruvp(3)=v_R
            ruvp(4)=p_R    
        endif
    case(direct_y)
        !y 方向
        if(y .LT.  0.1_prec)then
            ruvp(1)=rho_L
            ruvp(2)=v_L
            ruvp(3)=u_L
            ruvp(4)=p_L
        elseif(y .GE. 0.1_prec)then
            ruvp(1)=1.0_prec + 0.2_prec*sin(50.0_prec*y)
            ruvp(2)=v_R
            ruvp(3)=u_R
            ruvp(4)=p_R    
        endif
    end select
end subroutine ShuOsher_init
    
subroutine Euler_1D_order_init(coor,time,ruvp)

    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
    real(prec) :: rho_L,u_L,v_L,p_L  
    real(prec) :: rho_R,u_R,v_R,p_R 
    real(prec) :: time
          
    u_R = 2.0_prec
    v_R = 0.0_prec
    p_R = 1.0_prec 

    
    x = coor(1)
    y = coor(2)
    select case(dire_shock)
    case(direct_x)
        
        !x 方向
        ruvp(1)=1.0_prec + 0.2_prec*sin(2.0_prec*pi*(x-time))
        ruvp(2)=u_R
        ruvp(3)=v_R
        ruvp(4)=p_R    
    case(direct_y)
        
        !y 方向
        ruvp(1)=1.0_prec + 0.2_prec*sin(2.0_prec*pi*(y-time))
        ruvp(2)=u_R
        ruvp(3)=v_R
        ruvp(4)=p_R  
    end select
end subroutine Euler_1D_order_init

    
subroutine Doublerarefaction_1D_init(coor,ruvp)

    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
    real(prec) :: rho_L,u_L,v_L,p_L
    real(prec) :: rho_R,u_R,v_R,p_R  
    
    rho_L = 7.0_prec
    u_L = -1.0_prec
    v_L = 0.0_prec
    p_L = 0.2_prec   
        
    rho_R = 7.0_prec
    u_R = 1.0_prec
    v_R = 0.0_prec
    p_R = 0.2_prec 




end subroutine Doublerarefaction_1D_init




subroutine Riemann2D_init(coor,ruvp)

    use global_var
    use parameter_setting
    implicit none
    
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
 
end subroutine Riemann2D_init

subroutine DoubleMach_init(coor,ruvp)
    !双马赫反射初始场
    !初始马赫数 M = 10
    !亚比alpha1 = p2/p1 = 116.5
    !温比alpha2 = T2/T1 = 20.388
    !波前波后速度比 / 密度比 alpha3 = V2/V1 = r2/r1 = 5.7143
    !U V根据角度和波前波后速比计算
    !   波前：静止。波上激波前速度无量纲化V1 波后V2 = V1*alpha3 
    !   波后：固定坐标系 V2' = V1 - V2 = V1(1-alpha3)
    !         u2 = V2'sin(θ) v2 = V2'cos(θ)
    !Refs:《流体力学》下册（吴望一著）
    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
        
    x = coor(1)
    y = coor(2)
    
    if(y < sqrt(3.0_prec)*(x - 1.0_prec/6.0_prec))then
        r = 1.4_prec
        u = 0.0_prec
        v = 0.0_prec
        p = 1.0_prec
    elseif(y >= sqrt(3.0_prec)*(x - 1.0_prec/6.0_prec))then
        r = 8.0_prec
        u = 7.145_prec
        v = -4.125_prec
        p = 116.5_prec   
    endif
    !Refs: [1]Guo J, Zhu H, Yan Z-G, et al. High-Order Hybrid WCNS-CPR Scheme for Shock Capturing of Conservation Laws [J]. 
    !         International Journal of Aerospace Engineering, 2020, 2020: 1-13.
    ruvp(1) = r
    ruvp(2) = u
    ruvp(3) = v
    ruvp(4) = p
end subroutine DoubleMach_init
    
subroutine VortexShock_init(coor,ruvp)
    !涡-激波相互作用初始场
    !初始马赫数 M = 1.1
    !Refs:《计算流体力学——典型算法与算例》（高歌等著）p80
    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r0,u0,v0,p0,r,u,v,p,ruvp(1:4),Ma,xcore,ycore,epsilon,RR0,alpha,T0,RB2,RB,temexp,Te,theta,tao,tt,ratio
    real(prec),dimension(:) :: coor(2)
        
    x = coor(1)
    y = coor(2)
    !激波强度
    Ma = 1.1_prec
    Ms = Ma
    !激波前参数
    r0 = 1.0_prec
    u0 = Ma*sqrt(gamma)
    v0 = 0.0_prec
    p0 = 1.0_prec
           
    if(x < 2.0_prec)then
        !激波前添加涡
        !涡的位置
        xcore = 0.25_prec
        ycore = 0.5_prec
        epsilon = 0.3_prec
        RR0 = 0.05_prec
        alpha = 0.204_prec
        T0 = p0/r0
        
        RB2 =  (x-xcore)**2 + (y-ycore)**2
        RB  = sqrt(RB2)
        tao = RB/RR0
        temexp = exp(alpha*(1.0_prec-tao**2))
        Te = -(gamma-1.0_prec)*epsilon*epsilon*temexp*temexp/(4.0_prec*alpha*gamma)
        theta = atan2((y-ycore),(x-xcore))
        
        r = (T0 + TE)**(1.0_prec/(gamma-1.0_prec))
        p = r**gamma
        u = u0 + tao*epsilon*temexp*sin(theta)
        v = v0 - tao*epsilon*temexp*cos(theta)
        
    else
        !激波后参数
        ratio = ((gamma+1.0_prec)*Ma**2/(2.0_prec+(gamma-1.0_prec)*Ma**2))
        r = r0*ratio
        u = u0/ratio
        v = 0.0_prec
        p = p0*(1.0_prec + 2.0_prec*gamma*(Ma**2-1.0_prec)/(gamma+1.0_prec))   

    endif

    ruvp(1) = r
    ruvp(2) = u
    ruvp(3) = v
    ruvp(4) = p
end subroutine VortexShock_init
    
subroutine CompositeVortexShock_init(coor,ruvp)
    !涡-激波相互作用初始场
    !初始马赫数 M = 1.1
    !Refs:《计算流体力学——典型算法与算例》（高歌等著）p80
    !   Audrey Rault, Guillaume Chiavassa, and Rosa Donat. Shock–vortex interactions at high mach numbers. Journal of Scientific Computing, 19(1-3):347–371, 2003.
    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r0,u0,v0,p0,r,u,v,p,ruvp(1:4),xcore,ycore,epsilon,RR,alpha,T0,RB2,RB,temexp,Te,theta,tao,tt,ratio,Vm,constant1,constant2,constant22
    real(prec) :: a,b,ra,Vtheta,C1,C2
    real(prec),dimension(:) :: coor(2)
        
    x = coor(1)
    y = coor(2)
    !激波强度
    
    !Mv = 1.7_prec
    !Ms = 1.1_prec
    
    a = 0.075_prec
    b = 0.175_prec
    RR = 1.0_prec
    Vm = Mv*sqrt(gamma)
    
    !激波前参数
    r0 = 1.0_prec
    u0 = Ms*sqrt(gamma)
    v0 = 0.0_prec
    p0 = 1.0_prec
           
    if(x < 0.5_prec)then

        !涡的位置
        xcore = 0.25_prec
        ycore = 0.5_prec
      
        T0 = p0/r0  
        RB2 =  (x-xcore)**2 + (y-ycore)**2
        RB  = sqrt(RB2)
        constant1 = (gamma-1.0_prec)/(RR*gamma)
        constant2 = Vm*(a/(a**2-b**2))
        constant22 = constant2**2
        if(RB<=a)then
            Vtheta = Vm*RB/a
            C1 = constant1*constant22*0.5_prec*(a**2-4.0_prec*b**2*dlog(a)-b**4/(a**2)+4.0_prec*b**2*dlog(b))-0.5_prec*constant1*Vm**2
            Te = constant1*Vtheta**2*0.5_prec + C1          
        elseif(RB>a .and. RB < b)then       
            
            Vtheta = constant2*(RB-b**2/RB)            
            C2 = 2.0_prec*b**2*dlog(b)*constant1*constant22
            Te = constant1*constant22*0.5_prec*(RB**2 - 4.0_prec*(b**2)*dlog(RB)-(b**4/RB**2))+C2           
        elseif(RB>=b)then
            Vtheta = 0.0_prec
            Te = 0.0_prec
        end if

        theta = atan2((y-ycore),(x-xcore))    
        r = (T0+Te)**(1.0_prec/(gamma-1.0_prec))  
        p = r**gamma
        u = u0 - Vtheta*sin(theta)
        v = v0 + Vtheta*cos(theta)       
        !u = u0 - Vtheta*(y-ycore)/RB
        !v = v0 + Vtheta*(x-xcore)/RB
    else
        
        !激波后参数
        ratio = ((gamma+1.0_prec)*Ms**2/(2.0_prec+(gamma-1.0_prec)*Ms**2))
                
        r = r0*ratio
        u = u0/ratio
        v = 0.0_prec
        p = p0*(1.0_prec + 2.0_prec*gamma*(Ms**2-1.0_prec)/(gamma+1.0_prec))   
    
    endif

    ruvp(1) = r
    ruvp(2) = u
    ruvp(3) = v
    ruvp(4) = p
    
end subroutine CompositeVortexShock_init
    
subroutine steadyShock_init(coor,ruvp)

    ! 稳态激波初始场

    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r0,u0,v0,p0,r,u,v,p,ruvp(1:4),ratio
    real(prec),dimension(:) :: coor(2)
        
    x = coor(1)
    y = coor(2)

    select case(dire_shock)
    case(direct_x)
        ! x 方向
        ! 激波前参数
        r0 = 1.0_prec
        u0 = Ms*sqrt(gamma)
        v0 = 0.0_prec
        p0 = 1.0_prec
        if(x < 0.5_prec)then   
            ! 测试定激波        
            r = r0
            p = p0
            u = u0
            v = v0
        else
            ! 激波后参数
            ratio = (gamma+1.0_prec)*Ms**2/(2.0_prec+(gamma-1.0_prec)*Ms**2)         
            r = r0*ratio
            u = u0/ratio
            v = v0
            p = p0*(1.0_prec + 2.0_prec*gamma*(Ms**2-1.0_prec)/(gamma+1.0_prec))   
        endif   
    case(direct_y)
        ! y 方向
        ! 激波前参数
        r0 = 1.0_prec
        u0 = 0.0_prec
        v0 = Ms*sqrt(gamma)
        p0 = 1.0_prec
        if(y < 0.5_prec)then   
            ! 测试定激波        
            r = r0
            p = p0
            u = u0
            v = v0
        else
            ! 激波后参数
            ratio = (gamma+1.0_prec)*Ms**2/(2.0_prec+(gamma-1.0_prec)*Ms**2)         
            r = r0*ratio
            u = u0
            v = v0/ratio
            p = p0*(1.0_prec + 2.0_prec*gamma*(Ms**2-1.0_prec)/(gamma+1.0_prec))   
        endif  
    end select       

    ruvp(1) = r
    ruvp(2) = u
    ruvp(3) = v
    ruvp(4) = p
end subroutine steadyShock_init
    
subroutine HyperCylinder_case_init(coor,ruvp)

    ! 圆柱绕流

    use global_var
    use parameter_setting
    implicit none
    
    real(prec) :: x,y,r0,u0,v0,p0,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
        
    x = coor(1)
    y = coor(2)
    
    r0 = 1.4_prec
    u0 = Ms
    v0 = 0.0_prec
    p0 = 1.0_prec
           
    ruvp(1) = r0
    ruvp(2) = u0
    ruvp(3) = v0
    ruvp(4) = p0
end subroutine HyperCylinder_case_init

subroutine WingFlow_case_init(coor,ruvp)

    ! 机翼绕流

    use global_var
    use parameter_setting
    implicit none
    
    real(prec) :: x,y,r0,u0,v0,p0,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)

    x = coor(1)
    y = coor(2)
    !激波强度 

    r0 = 1.4_prec
    u0 = Ms*cos(attackAngle*PI/180.0_prec )
    v0 = Ms*sin(attackAngle*PI/180.0_prec )
    p0 = 1.0_prec
    
    ruvp(1) = r0
    ruvp(2) = u0
    ruvp(3) = v0
    ruvp(4) = p0

end subroutine WingFlow_case_init    
    
subroutine test_init(coor,ruvp)

    ! 测试专用

    use global_var
    use parameter_setting
    implicit none
    real(prec) :: x,y,r,u,v,p,ruvp(1:4)
    real(prec),dimension(:) :: coor(2)
    real(prec) :: rho_L,u_L,v_L,p_L 
    real(prec) :: rho_R,u_R,v_R,p_R  
    
    
    !rho_L = 3.857143_prec
    u_L = 1.0_prec
    v_L = 0.0_prec
    p_L = 1.0_prec    !初始时刻，分隔面左侧参数
        
    
    u_R = 1.0_prec
    v_R = 0.0_prec
    p_R = 1.0_prec !初始时刻，分隔面右侧参数

    x = coor(1)
    y = coor(2)
    ruvp(1)=1.0_prec + 0.2_prec*sin(10.0_prec*pi*x)
    ruvp(2)=u_L
    ruvp(3)=v_L
    ruvp(4)=p_L
    !select case(dire_shock)
    !case(direct_x)
    !    !x 方向
    !    if(x .LT. 0.5_prec)then
    !        ruvp(1)=1.0_prec + 0.2_prec*sin(10.0_prec*pi*x)
    !        ruvp(2)=u_L
    !        ruvp(3)=v_L
    !        ruvp(4)=p_L
    !    elseif(x .GE. 0.5_prec)then
    !        ruvp(1)=1.0_prec + 0.2_prec*sin(10.0_prec*pi*x)
    !        ruvp(2)=u_R
    !        ruvp(3)=v_R
    !        ruvp(4)=p_R    
    !    endif
    !case(direct_y)
    !    !y 方向
    !    if(y .LT.0.1_prec)then
    !        ruvp(1)=1.0_prec + 0.2_prec*sin(10.0_prec*pi*y)
    !        ruvp(2)=v_L
    !        ruvp(3)=u_L
    !        ruvp(4)=p_L
    !    elseif(y .GE. 0.1_prec)then
    !        ruvp(1)=1.0_prec + 0.2_prec*sin(10.0_prec*pi*y)
    !        ruvp(2)=v_R
    !        ruvp(3)=u_R
    !        ruvp(4)=p_R    
    !    endif
    !end select
end subroutine test_init