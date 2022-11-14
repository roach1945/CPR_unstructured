!----整体---------------------------------------------------------------------------------------------------
subroutine ori_to_con
    !原始变量装换到守恒变量
    use global_var
    use type_module
    implicit none
    integer i
    real(prec),dimension(:,:,:) :: r(nsp,nsp),u(nsp,nsp),v(nsp,nsp),p(nsp,nsp)
    
    do i =1,ncells+nbdsides
        r = cellset(i).spvalue_ori(:,:,1)
        u = cellset(i).spvalue_ori(:,:,2)
        v = cellset(i).spvalue_ori(:,:,3)
        p = cellset(i).spvalue_ori(:,:,4)
        
        cellset(i).spvalue_con(:,:,1) = r                                           !rho
        cellset(i).spvalue_con(:,:,2) = r*u                                         !rho*u
        cellset(i).spvalue_con(:,:,3) = r*v                                         !rho*v
        cellset(i).spvalue_con(:,:,4) = p/(gamma-1.0_prec)+0.5_prec*r*(u**2+v**2)   !E
        
    end do
    
end subroutine ori_to_con

subroutine ori_to_flu
    !原始变量转换到对流项
    use global_var
    use type_module
    implicit none
    integer i
    real(prec),dimension(:,:,:) :: r(nsp,nsp),u(nsp,nsp),v(nsp,nsp),p(nsp,nsp)
    real(prec),dimension(:,:,:) :: con1(nsp,nsp),con2(nsp,nsp),con3(nsp,nsp),con4(nsp,nsp)
    
    do i =1,ncells+nbdsides
        r = cellset(i).spvalue_ori(:,:,1)
        u = cellset(i).spvalue_ori(:,:,2)
        v = cellset(i).spvalue_ori(:,:,3)
        p = cellset(i).spvalue_ori(:,:,4)
        con1 = cellset(i).spvalue_con(:,:,1)                                    !rho
        con2 = cellset(i).spvalue_con(:,:,2)                                    !rho*u
        con3 = cellset(i).spvalue_con(:,:,3)                                    !rho*v
        con4 = cellset(i).spvalue_con(:,:,4)                                    !E
        !F(ru,ru^2+p,ruv,u(E+p))
        cellset(i).spvalue_fluF(:,:,1) = con2                                   !r*u
        cellset(i).spvalue_fluF(:,:,2) = con2*u+p                               !r*u^2+p
        cellset(i).spvalue_fluF(:,:,3) = con2*v                                 !r*u*v
        cellset(i).spvalue_fluF(:,:,4) = u*(con4+p)                             !u(E+p)
        !G(rv,ruv,rv^2+p,v(E+p))
        cellset(i).spvalue_fluG(:,:,1) = con3                                   !r*v
        cellset(i).spvalue_fluG(:,:,2) = con3*u                                 !r*u*v
        cellset(i).spvalue_fluG(:,:,3) = con3*v+p                               !r*v^2+p
        cellset(i).spvalue_fluG(:,:,4) = v*(con4+p)                             !v(E+p)
    end do
    
    
end subroutine ori_to_flu


subroutine con_to_ori
    !原始变量装换到守恒变量
    use global_var
    use type_module
    implicit none
    integer i
    real(prec),dimension(:,:,:) :: con1(nsp,nsp),con2(nsp,nsp),con3(nsp,nsp),con4(nsp,nsp)
    
    do i =1,ncells
        con1 = cellset(i).spvalue_con(:,:,1)                                    !rho
        con2 = cellset(i).spvalue_con(:,:,2)                                    !rho*u
        con3 = cellset(i).spvalue_con(:,:,3)                                    !rho*v
        con4 = cellset(i).spvalue_con(:,:,4)                                    !E
        
        cellset(i).spvalue_ori(:,:,1) = con1                                    !rho
        cellset(i).spvalue_ori(:,:,2) = con2/con1                               !u
        cellset(i).spvalue_ori(:,:,3) = con3/con1                               !v
        cellset(i).spvalue_ori(:,:,4) = (gamma-1.0_prec)*(con4-0.5_prec*con1*(cellset(i).spvalue_ori(:,:,2)**2+cellset(i).spvalue_ori(:,:,3)**2))!p

    end do   

end subroutine con_to_ori

subroutine phy_to_com
    !求解计算域下各变量的值，conservation,F,G
    use global_var
    use type_module
    implicit none
    integer i,k,l   
    real(prec),dimension(:,:)   :: detJ(nsp,nsp) 
    real(prec) :: kesix,kesiy,etax,etay
    
    do i = 1,ncells+nbdsides
        detJ  = cellset(i).det_J        
        do k = 1,nsp
            do l =1,nsp           
                kesix = cellset(i).Mdirect(k,l,1)
                kesiy = cellset(i).Mdirect(k,l,2)
                etax  = cellset(i).Mdirect(k,l,3)
                etay  = cellset(i).Mdirect(k,l,4)
                cellset(i).spvalue_con_loc(k,l,:)  = detJ(k,l)*cellset(i).spvalue_con(k,l,:)
                cellset(i).spvalue_fluF_loc(k,l,:) = detJ(k,l)*(kesix*cellset(i).spvalue_fluF(k,l,:) + kesiy*cellset(i).spvalue_fluG(k,l,:))
                cellset(i).spvalue_fluG_loc(k,l,:) = detJ(k,l)*(etax *cellset(i).spvalue_fluF(k,l,:) + etay *cellset(i).spvalue_fluG(k,l,:))  
            end do
        end do  
    end do
    
end subroutine phy_to_com


!----局部---------------------------------------------------------------------------------------------------
subroutine Func_ori_to_flu(ori,f,g)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p,f(4),g(4),ori(4)
    real(prec) ::va
    r=ori(1)
    u=ori(2)
    v=ori(3)
    p=ori(4)
    va=u*u+v*v
    f(1)=r*u
    f(2)=r*u*u+p
    f(3)=r*u*v
    f(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*u

    g(1)=r*v
    g(2)=r*u*v
    g(3)=r*v*v+p
    g(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*v

end subroutine Func_ori_to_flu
subroutine Func_ori_to_fluxF(ori,f)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p,f(4),g(4),ori(4)
    real(prec) ::va
    r=ori(1)
    u=ori(2)
    v=ori(3)
    p=ori(4)
    va=u*u+v*v
    f(1)=r*u
    f(2)=r*u*u+p
    f(3)=r*u*v
    f(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*u

end subroutine Func_ori_to_fluxF
subroutine Func_ori_to_fluxG(ori,g)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p,f(4),g(4),ori(4)
    real(prec) ::va
    r=ori(1)
    u=ori(2)
    v=ori(3)
    p=ori(4)
    va=u*u+v*v

    g(1)=r*v
    g(2)=r*u*v
    g(3)=r*v*v+p
    g(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*v
end subroutine Func_ori_to_fluxG

subroutine Func_ori_to_con(ori,conser)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p
    real(prec) ::ori(4),conser(4)
    r=ori(1)
    u=ori(2)
    v=ori(3)
    p=ori(4)
    conser(1)=r
    conser(2)=r*u
    conser(3)=r*v
    conser(4)=p/gamma1+0.5_prec*r*(u*u+v*v) 

end subroutine Func_ori_to_con

subroutine Func_con_to_ori(conser,ori)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p
    real(prec) ::ori(4),conser(4)
    
    r=conser(1)
    u=conser(2)/r
    v=conser(3)/r
    p=(conser(4)-0.5_prec*r*(u*u+v*v))*gamma1 
    ori(1)=r
    ori(2)=u
    ori(3)=v
    ori(4)=p
    
end subroutine Func_con_to_ori
subroutine Func_con_to_fluxF(conser,fluxF)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    
    real(prec) ::ori(4),conser(4),fluxF(4),fluxG(4)
    
    call Func_con_to_ori(conser,ori)
    
    call Func_ori_to_fluxF(ori,fluxF)
    
end subroutine Func_con_to_fluxF
    
subroutine Func_con_to_fluxG(conser,fluxG)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    
    real(prec) ::ori(4),conser(4),fluxF(4),fluxG(4)
    
    call Func_con_to_ori(conser,ori)
    
    call Func_ori_to_fluxG(ori,fluxG)
    
end subroutine Func_con_to_fluxG
subroutine Func_PhyTranCom(index,detJ,ori,oriCom)
    !获得计算域上的原始变量
    use real_precision
    implicit none
    integer :: index
    real(prec) ::ori(4),oriCom(4),detJ

    oriCom(1) = ori(1)!*detJ        !注释掉，采取原始变量。计算空间的原始变量
    oriCom(2) = ori(2)
    oriCom(3) = ori(3)
    oriCom(4) = ori(4)!*detJ
       
end subroutine Func_PhyTranCom
subroutine Func_ComTranPhy(index,detJ,Com,Phy)
    !获得物理域上的原始变量
    use real_precision
    use parameter_setting
    implicit none
    integer :: index
    real(prec) ::ori(4),oriCom(4),Com(4),detJ,Phy(4)

    oriCom = Com
    Phy(1) = oriCom(1)!/detJ
    Phy(2) = oriCom(2)
    Phy(3) = oriCom(3)
    Phy(4) = oriCom(4)!/detJ

end subroutine Func_ComTranPhy

subroutine proj_matrix(ori,ll,rr,director)   
    !特征投影左特征向量矩阵，右特征向量矩阵
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none  
    real(prec) ::ori(1:4),ll(1:4,1:4),rr(1:4,1:4)
    real(prec) :: r,u,v,p,H,c,c2,Vector
    integer :: director !0 -- F;1 -- G
    
    r = ori(1)
    u = ori(2)
    v = ori(3)
    p = ori(4)
    
    c2 = gamma*p/r    !c^2
    c = sqrt(c2)      !c
    Vector = (u**2.0_prec + v**2.0_prec) 
    H =  0.5_prec*Vector + c2/gamma1
    if(mod(director,2) == 0)then   
    !A(U)
        ll(1,1)=  H + c*(u-c)/gamma1
        ll(1,2)= -(u + c/gamma1)
        ll(1,3)= -v
        ll(1,4)=  1.0_prec
     
        ll(2,1)= -2.0_prec*H + 4.0_prec*c2/gamma1
        ll(2,2)=  2.0_prec*u
        ll(2,3)=  2.0_prec*v
        ll(2,4)= -2.0_prec
     
        ll(3,1)= -2.0_prec*v*c2/gamma1
        ll(3,2)=  0.0_prec
        ll(3,3)=  2.0_prec*c2/gamma1
        ll(3,4)=  0.0_prec     
     
        ll(4,1)=  H - c*(u+c)/gamma1
        ll(4,2)= -u + c/gamma1
        ll(4,3)= -v
        ll(4,4)=  1.0_prec

        !!
        rr(1,1)=  1.0_prec
        rr(1,2)=  1.0_prec
        rr(1,3)=  0.0_prec
        rr(1,4)=  1.0_prec
     
        rr(2,1)=  u - c
        rr(2,2)=  u
        rr(2,3)=  0.0_prec
        rr(2,4)=  u + c
     
        rr(3,1)=  v
        rr(3,2)=  v
        rr(3,3)=  1.0_prec
        rr(3,4)=  v  
     
        rr(4,1)=  H - c*u
        rr(4,2)=  0.5_prec*Vector
        rr(4,3)=  v
        rr(4,4)=  H + c*u
        
        rr = gamma1/(2.0_prec*c2)*rr
    elseif(mod(director,2) == 1)then
    !B(U)
        ll(1,1)=  H + c*(v-c)/gamma1
        ll(1,2)= -u
        ll(1,3)= -(v + c/gamma1)
        ll(1,4)=  1.0_prec
        
        ll(2,1)= -2.0_prec*u*c2/gamma1
        ll(2,2)=  2.0_prec*c2/gamma1
        ll(2,3)=  0.0_prec
        ll(2,4)=  0.0_prec       
        
        ll(3,1)= -2.0_prec*H + 4.0_prec*c2/gamma1
        ll(3,2)=  2.0_prec*u
        ll(3,3)=  2.0_prec*v
        ll(3,4)= -2.0_prec 
     
        ll(4,1)=  H - c*(v+c)/gamma1
        ll(4,2)= -u 
        ll(4,3)= -v + c/gamma1
        ll(4,4)=  1.0_prec
     
        rr(1,1)=  1.0_prec
        rr(1,2)=  0.0_prec
        rr(1,3)=  1.0_prec
        rr(1,4)=  1.0_prec
     
        rr(2,1)=  u
        rr(2,2)=  1.0_prec
        rr(2,3)=  u
        rr(2,4)=  u
     
        rr(3,1)=  v - c
        rr(3,2)=  0.0_prec
        rr(3,3)=  v
        rr(3,4)=  v + c 
     
        rr(4,1)=  H - c*v
        rr(4,2)=  u
        rr(4,3)=  0.5_prec*Vector
        rr(4,4)=  H + c*v    
        
        rr = gamma1/(2.0_prec*c2)*rr

    end if
    !write(*,*) ll
end subroutine    

subroutine proj_matrix2(ruvp,ll,rr,director)   
    !特征投影左特征向量矩阵，右特征向量矩阵
    !贾振岭对照组
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none  
    
    integer :: director !0 -- F;1 -- G
    real(prec) :: ruvp(1:4),ll(1:4,1:4),rr(1:4,1:4),c,c2,r,u,v,p
    real(prec) :: nx,ny,q2,qn,lx,ly,qL,H
    integer :: i
    
    r = ruvp(1)
    u = ruvp(2)
    v = ruvp(3)
    p = ruvp(4)

    if(mod(director,2) == 0)then   
    !A(U)      
        nx = 1.0_prec
        ny = 0.0_prec
    elseif(mod(director,2) == 1)then
    !B(U) 
        nx = 0.0_prec
        ny = 1.0_prec
    end if
    
        lx = -ny
        ly = nx
        qn = u*nx + v*ny
        qL = u*lx + v*ly
        q2 = u**2 + v**2
        
        c2=gamma*p/r
        c=sqrt(c2)       
        H =  0.5_prec*q2 + c2/gamma1
        ll(1,1)=  0.5_prec*(gamma1/2.0_prec/c2*q2 + qn/c)
        ll(1,2)= -0.5_prec*(gamma1/c2*u + nx/c)
        ll(1,3)= -0.5_prec*(gamma1/c2*v + ny/c)
        ll(1,4)=  gamma1/2.0_prec/c2
     
        ll(2,1)=  1.0_prec - gamma1/2.0_prec/c2*q2
        ll(2,2)=  gamma1/c2*u
        ll(2,3)=  gamma1/c2*v
        ll(2,4)= -gamma1/c2
     
        ll(3,1)= 0.5_prec*(gamma1/2.0_prec/c2*q2 - qn/c)
        ll(3,2)= -0.5_prec*(gamma1/c2*u - nx/c)
        ll(3,3)= -0.5_prec*(gamma1/c2*v - ny/c)
        ll(3,4)= gamma1/2.0_prec/c2
     
        ll(4,1)= -qL
        ll(4,2)= lx
        ll(4,3)= ly
        ll(4,4)= 0.0_prec
     
     
        rr(1,1)=  1.0_prec
        rr(1,2)=  1.0_prec
        rr(1,3)=  1.0_prec
        rr(1,4)=  0.0_prec
     
        rr(2,1)=  u - c*nx
        rr(2,2)=  u
        rr(2,3)=  u + c*nx
        rr(2,4)=  lx
     
        rr(3,1)=  v - c*ny
        rr(3,2)=  v
        rr(3,3)=  v + c*ny
        rr(3,4)=  ly
     
        rr(4,1)=  H - qn*c
        rr(4,2)=  q2/2.0_prec
        rr(4,3)=  H + qn*c
        rr(4,4)=  qL
    !write(*,*) ll
end subroutine 
!subroutine proj_matrix_x(ruvwp,ll,rr)   
!    use real_precision
!     use global_variables,only: gamma,gamma1
!     implicit none  
!     real(prec) ::ruvwp(1:5),ll(1:5,1:5),rr(1:5,1:5),cc,c2,ae,alfa
!     integer :: i
!     c2=gamma*ruvwp(5)/ruvwp(1)
!     cc=sqrt(c2)
!     ae = gamma1
!     alfa=0.5_prec*(ruvwp(2)**2.0_prec + ruvwp(3)**2.0_prec + ruvwp(4)**2.0_prec) 
!     
!     ll(1,1)=  ae*alfa-c2
!     ll(1,2)= -ae*ruvwp(2)
!     ll(1,3)= -ae*ruvwp(3)
!     ll(1,4)= -ae*ruvwp(4)
!     ll(1,5)=  ae
!     
!     ll(2,1)= -ruvwp(4)
!     ll(2,2)=  0.0_prec
!     ll(2,3)=  0.0_prec
!     ll(2,4)=  1.0_prec
!     ll(2,5)=  0.0_prec
!     
!     ll(3,1)= ruvwp(3)
!     ll(3,2)= 0.0_prec
!     ll(3,3)= -1.0_prec
!     ll(3,4)= 0.0_prec
!     ll(3,5)= 0.0_prec
!     
!     ll(4,1)= ae*alfa-cc*ruvwp(2)
!     ll(4,2)= cc-ae*ruvwp(2)
!     ll(4,3)= -ae*ruvwp(3)
!     ll(4,4)= -ae*ruvwp(4)
!     ll(4,5)=  ae
!     
!     ll(5,1)= ae*alfa+cc*ruvwp(2)
!     ll(5,2)=-cc-ae*ruvwp(2)
!     ll(5,3)=-ae*ruvwp(3)
!     ll(5,4)=-ae*ruvwp(4)
!     ll(5,5)= ae
!     
!     
!      rr(1,1) = -2.0_prec
!      rr(2,1) = -2.0_prec*ruvwp(2)
!      rr(3,1) = -2.0_prec*ruvwp(3)
!      rr(4,1) = -2.0_prec*ruvwp(4)
!      rr(5,1) = -2.0_prec*alfa
!                
!      rr(1,2) = 0.0_prec
!      rr(2,2) = 0.0_prec
!      rr(3,2) = 0.0_prec
!      rr(4,2) = 2.0_prec*c2
!      rr(5,2) = 2.0_prec*ruvwp(4)*c2 
!                
!      rr(1,3) = 0.0_prec
!      rr(2,3) = 0.0_prec
!      rr(3,3) =-2.0_prec*c2
!      rr(4,3) = 0.0_prec
!      rr(5,3) =-2.0_prec*ruvwp(3)*c2 
!                
!      rr(1,4) = 1.0_prec
!      rr(2,4) = ruvwp(2)+cc
!      rr(3,4) = ruvwp(3)
!      rr(4,4) = ruvwp(4)
!      rr(5,4) = alfa+cc*ruvwp(2)+c2/ae
!                
!      rr(1,5) = 1.0_prec
!      rr(2,5) = ruvwp(2)-cc
!      rr(3,5) = ruvwp(3)
!      rr(4,5) = ruvwp(4)
!      rr(5,5) = alfa-cc*ruvwp(2)+c2/ae
!      rr=rr*0.5_prec/c2
!end subroutine proj_matrix_x  
! 
!subroutine proj_matrix_y(ruvwp,ll,rr)   
!       use real_precision
!     use global_variables,only: gamma,gamma1
!     implicit none  
!     real(prec) ::ruvwp(1:5),ll(1:5,1:5),rr(1:5,1:5),cc,c2,ae,alfa
!     integer :: i
!     c2=gamma*ruvwp(5)/ruvwp(1)
!     cc=sqrt(c2)
!     ae = gamma1
!     alfa=0.5_prec*(ruvwp(2)**2.0_prec + ruvwp(3)**2.0_prec + ruvwp(4)**2.0_prec) 
!     ll=0.0_prec
!     rr=0.0_prec
!
!    ll(3,1)= ruvwp(4)
!    ll(3,2)= 0.0_prec
!    ll(3,3)= 0.0_prec
!    ll(3,4)= -1.0_prec
!    ll(3,5)= 0.0_prec   
!    
!    ll(1,1)=  ae*alfa-c2
!    ll(1,2)= -ae*ruvwp(2)
!    ll(1,3)= -ae*ruvwp(3)
!    ll(1,4)= -ae*ruvwp(4)
!    ll(1,5)=  ae
!    
!    ll(2,1)=  ruvwp(2)
!    ll(2,2)=  -1.0_prec
!    ll(2,3)=  0.0_prec
!    ll(2,4)=  0.0_prec
!    ll(2,5)=  0.0_prec
!    
!    
!    
!    ll(4,1)= ae*alfa-cc*ruvwp(3)
!    ll(4,2)= -ae*ruvwp(2)
!    ll(4,3)= cc-ae*ruvwp(3)
!    ll(4,4)= -ae*ruvwp(4)
!    ll(4,5)=  ae
!    
!    ll(5,1)= ae*alfa+cc*ruvwp(3)
!    ll(5,2)=-ae*ruvwp(2)
!    ll(5,3)=-cc-ae*ruvwp(3)
!    ll(5,4)=-ae*ruvwp(4)
!    ll(5,5)= ae
!    
!    
!                 
!     rr(1,3) = 0.0_prec
!     rr(2,3) = 0.0_prec
!     rr(3,3) = 0.0_prec
!     rr(4,3) = -2.0_prec*c2
!     rr(5,3) = -2.0_prec*ruvwp(4)*c2 
!    
!     rr(1,1) = -2.0_prec
!     rr(2,1) = -2.0_prec*ruvwp(2)
!     rr(3,1) = -2.0_prec*ruvwp(3)
!     rr(4,1) = -2.0_prec*ruvwp(4)
!     rr(5,1) = -2.0_prec*alfa
!               
!     rr(1,2) = 0.0_prec
!     rr(2,2) = -2.0_prec*c2
!     rr(3,2) = 0.0_prec
!     rr(4,2) = 0.0_prec
!     rr(5,2) = -2.0_prec*ruvwp(2)*c2 
!    
!     rr(1,4) = 1.0_prec
!     rr(2,4) = ruvwp(2)
!     rr(3,4) = ruvwp(3)+cc
!     rr(4,4) = ruvwp(4)
!     rr(5,4) = alfa+cc*ruvwp(3)+c2/ae
!               
!     rr(1,5) = 1.0_prec
!     rr(2,5) = ruvwp(2)
!     rr(3,5) = ruvwp(3)-cc
!     rr(4,5) = ruvwp(4)
!     rr(5,5) = alfa-cc*ruvwp(3)+c2/ae
!     rr=rr*0.5_prec/c2
!end subroutine proj_matrix_y     
subroutine Characteristic_projection(ori,ll,chara_var)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::conser(1:4),chara_var(1:4),ori(1:4),ll(1:4,1:4)
    integer :: i,j
    
    call Func_ori_to_con(ori,conser) 

    do i=1,4
        chara_var(i)=0.0_prec
        do j=1,4
            chara_var(i)=chara_var(i) + ll(i,j)*conser(j)          
        end do
    end do
    
end subroutine Characteristic_projection       

subroutine Inverse_Characteristic_projection(interface_chara_var,rr,ori)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) :: interface_chara_var(1:4),conser(1:4),ori(1:4)
    real(prec) :: rr(1:4,1:4)
    integer :: i,j,k
    do i=1,4
        conser(i)=0.0_prec
        do j=1,4 
            conser(i)=conser(i) + rr(i,j)*interface_chara_var(j)
        end do
    end do
     
    call Func_con_to_ori(conser,ori)

end subroutine Inverse_Characteristic_projection

 