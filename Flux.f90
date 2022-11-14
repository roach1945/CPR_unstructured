!各种迎风通量求解方法，Roe, Godunov, L-F
subroutine getRiemannFlux(ori1,ori2,norm,fupw)
    use real_precision 
    use parameter_setting
    implicit none
    real(prec) :: fupw(4),ori1(4),ori2(4),norm(2)
    select case(Riemann_flux)
    case(LLF_flux)
        call laxf_flux(ori1,ori2,norm,fupw)
    case(Roe_flux)
        call Roe_FDS(ori1,ori2,norm,fupw)
    end select
    
    
end subroutine getRiemannFLux
!---Roe Flux-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function Roe_Flux(ul,ur,fl,fr)
    !一维情况
    use real_precision
    implicit none
    real(prec) :: Roe_Flux,ul,ur,fl,fr,a
    if(ur .EQ. ul)then
        a = 0.0_prec!避免ur=ul，得到无穷大数
    else
        a = abs((fr-fl)/(ur-ul))
    end if
    Roe_Flux = 0.5_prec*((fl+fr)-a*(ur-ul))
    return
end function Roe_Flux 

subroutine Roe_FDS(ql,qr,norm,ff)
    use real_precision
    use global_var,only: gamma,gamma1
    use parameter_setting,only:efix
    implicit none
    integer :: i,m,nl
    real(prec) :: ff(1:4),df(1:4),ql(1:4),qr(1:4),qroe(1:4),dq(1:4),norm(2)
    real(prec) :: l1,l3,l4,rm,um,vm,wm,cm,hm,l12,l32,l42
    real(prec) :: nx,ny,nz,nt,cta,cgm,cgm1,rrl,rrl_1,ccgm,nn
    real(prec) :: c2,c,dcta,cta1,gama2,deti,det2,rodctac
    real(prec) :: rodcta,dp_c2,a1,a2,a3,a4,a5,a6,a7
    real(prec) :: r1,u1,v1,p1,r2,u2,v2,p2,n1,n2,flr(4),f1(4),g1(4),f2(4),g2(4),conser1(4),conser2(4)
    nx=norm(1)
    ny=norm(2)
    nn = sqrt(nx*nx+ny*ny)
    nx=norm(1)/nn
    ny=norm(2)/nn

    r1=ql(1)
    u1=ql(2)   
    v1=ql(3)
    p1=ql(4)
   
    r2=qr(1)
    u2=qr(2)   
    v2=qr(3) 
    p2=qr(4)
    call flux_fgh(r1,u1,v1,p1,f1,g1) 
    call flux_fgh(r2,u2,v2,p2,f2,g2)  
    do i=1,4
    flr(i)=(f1(i)+f2(i))*nx+(g1(i)+g2(i))*ny
    end do

!   the difference of variables between right and left sides
    do m=1,4
        dq(m) = qr(m) - ql(m)
    enddo

!   to calculate H at left and right sides
    gama2 = gamma/gamma1
    ql(4) = gama2*ql(4)/ql(1) + 0.5_prec*(ql(2)*ql(2)+ql(3)*ql(3))
    qr(4) = gama2*qr(4)/qr(1) + 0.5_prec*(qr(2)*qr(2)+qr(3)*qr(3))

!   to calculate density of Roe average
    qroe(1)=sqrt(ql(1)*qr(1))

!   to calculate velocity and H of Roe average
    rrl = sqrt(qr(1)/ql(1))
    rrl_1 = rrl + 1.0_prec
    do m=2,4
        qroe(m)=(ql(m)+qr(m)*rrl)/rrl_1
    enddo

    rm = qroe(1)
    um = qroe(2)
    vm = qroe(3)
    hm = qroe(4)

!   to calculate the speed of sound
    c2 = (hm-0.5_prec*(um*um+vm*vm))*gamma1
    if(c2.LT.0.0_prec)then
        write(*,*)'sqrt(c2) is a math error in subroutine Roe_FDS'
        stop
    endif
    c = sqrt(c2)

    cta = um*nx+vm*ny+nt

    cgm = max(sqrt(nx*nx + ny*ny ),1.0e-30_prec)
    cgm1= 1.0_prec/cgm

    nx = nx*cgm1
    ny = ny*cgm1

    ccgm = c*cgm

    l1 = abs(cta-ccgm)
    l3 = abs(cta     )
    l4 = abs(cta+ccgm)

!   Harten's entropy modification
    l12   = l1*l1
    l32   = l3*l3
    l42   = l4*l4

    deti  = efix*cgm  
    det2  = deti*deti  

    l1 = sqrt(l12+det2)
    l3 = sqrt(l32+det2)
    l4 = sqrt(l42+det2)

    dcta= dq(2)*nx+dq(3)*ny

    cta1= cta*cgm1

    rodcta = qroe(1)*dcta
    dp_c2  = dq(4)/c2

    rodctac = rodcta/c

    a1 = l3*(dq(1)-dp_c2)
    a2 = l4*(dp_c2+rodctac)*0.5_prec
    a3 = l1*(dp_c2-rodctac)*0.5_prec
    a4 = a1+a2+a3
    a5 = c*(a2-a3)
    a6 = l3*(rm*dq(2)-nx*rodcta)
    a7 = l3*(rm*dq(3)-ny*rodcta)


    df(1) = a4
    df(2) = um*a4+nx*a5+a6
    df(3) = vm*a4+ny*a5+a7
    df(4) = hm*a4+cta1*a5+um*a6+vm*a7-c2*a1/gamma1

    do m=1,4
        ff(m) = 0.5_prec*(flr(m)-df(m))
    enddo
    ff(:) = ff(:)*nn
end subroutine Roe_FDS
    

subroutine flux_fgh(r,u,v,p,f,g)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p,f(4),g(4)
    real(prec) ::va
    va=u*u+v*v
    f(1)=r*u
    f(2)=r*u*u+p
    f(3)=r*u*v
    f(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*u

    g(1)=r*v
    g(2)=r*u*v
    g(3)=r*v*v+p
    g(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*v

end subroutine flux_fgh

subroutine get_conser_var(ruvwp,conser)
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec) ::r,u,v,p
    real(prec) ::ruvwp(4),conser(4)
    r=ruvwp(1)
    u=ruvwp(2)
    v=ruvwp(3)
    p=ruvwp(4)
    conser(1)=r
    conser(2)=r*u
    conser(3)=r*v
    conser(4)=p/gamma1+0.5_prec*r*(u*u+v*v) 

end subroutine get_conser_var

subroutine get_origi_var(conser,ruvwp)
    use real_precision
    use global_var,only: gamma,gamma1
    real(prec) ::ruvwp(4),conser(4)
    ruvwp(1)=conser(1)
    ruvwp(2)=conser(2)/conser(1)
    ruvwp(3)=conser(3)/conser(1)
    ruvwp(4)=gamma1*(conser(4)-0.5_prec*(ruvwp(2)*ruvwp(2)+ruvwp(3)*ruvwp(3))*ruvwp(1))
end subroutine get_origi_var
!---Roe Flux End-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!---Lax-Frid Flux-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

subroutine laxf_flux(ori1,ori2,norm,fupw)
    !二维
    use real_precision
    use global_var,only: gamma,gamma1
    implicit none
    real(prec),intrinsic:: abs
    real(prec),intrinsic:: sqrt
    real(prec),intrinsic::max      

    real(prec) :: r1,u1,v1,p1,r2,u2,v2,p2,n1,n2,n3,fupw(4),conser1(4),conser2(4),ori1(4),ori2(4),norm(2)
    real(prec) :: hh1,hh2,lhh,nn
    real(prec) :: lu,lv,lc,a0,f1(4),g1(4),f2(4),g2(4),lc1,lc2
    real(prec) :: c,cl,cr
    integer :: j
    
    !write(*,*)'ori',ori1
    call Func_ori_to_con(ori1,conser1) 
    call Func_ori_to_con(ori2,conser2) 

    r1=ori1(1)
    u1=ori1(2)   
    v1=ori1(3)
    p1=ori1(4)
   
    r2=ori2(1)
    u2=ori2(2)   
    v2=ori2(3)
    p2=ori2(4)   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !if(r1<0 .or.r2<0)then
    !    write(*,*)r1,r2
    !end if
    n1=norm(1)
    n2=norm(2)

    nn = sqrt(n1*n1+n2*n2)
    n1=norm(1)/nn
    n2=norm(2)/nn
    !write(*,*)'1', r1,r2
    lu=(u1*sqrt(r1)+u2*sqrt(r2))/(sqrt(r1)+sqrt(r2))
    lv=(v1*sqrt(r1)+v2*sqrt(r2))/(sqrt(r1)+sqrt(r2))
    !write(*,*) '2',lu,lv
    hh1=(gamma/gamma1)*p1/r1+0.5_prec*(u1*u1+v1*v1)
    hh2=(gamma/gamma1)*p2/r2+0.5_prec*(u2*u2+v2*v2)  
    !!------------------------------
    !lhh=(hh1*sqrt(r1)+hh2*sqrt(r2))/(sqrt(r1)+sqrt(r2))     
    !lc=sqrt(gamma1*(lhh-0.5_prec*(lu*lu+lv*lv))) 
    !if(lhh-0.5_prec*(lu*lu+lv*lv) < 0)then
    !    lc = 0.0_prec
    !end if
    !a0=abs(lu*n1+lv*n2)+lc
    !!------------------------------
    !write(*,*)hh1-0.5_prec*(u1*u1+v1*v1)
    lc1=sqrt(gamma1*(hh1-0.5_prec*(u1*u1+v1*v1)))
    lc2=sqrt(gamma1*(hh2-0.5_prec*(u2*u2+v2*v2)))
    a0=max(abs(u1*n1+v1*n2)+lc1,abs(u2*n1+v2*n2)+lc2 )
    !
    call Func_ori_to_flu(ori1,f1,g1) 
    call Func_ori_to_flu(ori2,f2,g2) 
    !write(*,*)ori1(1),ori2(1),n1,n2
    do j=1,4
        fupw(j)=0.5_prec*((f1(j)+f2(j))*n1+(g1(j)+g2(j))*n2 - a0*(conser2(j)-conser1(j)))
        !write(*,*) 0.5_prec*((f1(j)+f2(j))*n1+(g1(j)+g2(j))*n2)
    end do 
    
    fupw(:) = fupw(:)*nn
    !if(isnan(fupw(1)))then
    !    write(*,*)'0',r1,r2
    !    write(*,*)'1',lhh,hh1
    !    write(*,*)'2',lhh,lu*lu+lv*lv,gamma1*(lhh-0.5_prec*(lu*lu+lv*lv))
    !    write(*,*)'3',fupw,lu,lv,lc
    !end if
    !write(*,*)a0
    !write(*,"(2F10.5)")ori1(1),ori2(1)
    !write(*,"(5F10.5)")f1(1),f2(1),g1(1),g2(1),fupw(1)
    !write(*,*)'0',fupw(:)
end subroutine laxf_flux
!    
!subroutine laxf_flux(var1,var2,norm,flr)
!
!
!
!n1=norm(1)
!n2=norm(2)
!n3=norm(3) 
!nnorm = sqrt(N1*N1+N2*N2+N3*N3)
!n1=norm(1)/nnorm
!n2=norm(2)/nnorm
!n3=norm(3)/nnorm
!
!lu=(u1*sqrt(r1)+u2*sqrt(r2))/(sqrt(r1)+sqrt(r2))
!lv=(v1*sqrt(r1)+v2*sqrt(r2))/(sqrt(r1)+sqrt(r2))
!lw=(w1*sqrt(r1)+w2*sqrt(r2))/(sqrt(r1)+sqrt(r2))
!
!hh1=(gamma/gamma1)*p1/r1+0.5_prec*(u1*u1+v1*v1+w1*w1)
!hh2=(gamma/gamma1)*p2/r2+0.5_prec*(u2*u2+v2*v2+w2*w2)  
!lhh=(hh1*sqrt(r1)+hh2*sqrt(r2))/(sqrt(r1)+sqrt(r2)) 
!  
!!lc=sqrt(gamma1*(lhh-0.5_prec*(lu*lu+lv*lv+lw*lw)))
!!a0=abs(lu*n1+lv*n2+lw*n3)+lc
!!write(*,*) 111,sqrt(gamma1*(hh1-0.5_prec*(u1*u1+v1*v1+w1*w1))),gamma1*(hh1-0.5_prec*(u1*u1+v1*v1+w1*w1)),gamma1*(hh1-0.5_prec*(u1*u1+v1*v1+w1*w1))
!lc1=sqrt(gamma1*(hh1-0.5_prec*(u1*u1+v1*v1+w1*w1)))
!lc2=sqrt(gamma1*(hh2-0.5_prec*(u2*u2+v2*v2+w2*w2)))
!a0=max(abs(u1*n1+v1*n2+w1*n3)+lc1,abs(u2*n1+v2*n2+w2*n3)+lc2 )
!
!call flux_fgh(var1,f1,g1,h1) 
!call flux_fgh(var2,f2,g2,h2) 
!     
!do j=1,5
!flr(j)=0.5_prec*((f1(j)+f2(j))*n1+(g1(j)+g2(j))*n2+(h1(j)+h2(j))*n3 - a0*(conser2(j)-conser1(j)))
!!flr(j)=0.5_prec*((f1(j)+f2(j))*n1+(g1(j)+g2(j))*n2+(h1(j)+h2(j))*n3)
!end do 
!
!flr(:) = flr(:)*nnorm
!
!end subroutine laxf_flux
!******************************************************************************************************************
!subroutine Roe_FDS(conser1,conser2,norm,ff)
!    use real_precision
!    use global_variables,only: gamma,gamma1,efix
!    implicit none
!    integer :: i,m,nl
!    real(prec) :: ff(1:4),df(1:4),ql(1:4),qr(1:4),qroe(1:4),dq(1:4),norm(2)
!    real(prec) :: l1,l3,l4,rm,um,vm,wm,cm,hm,l12,l32,l42
!    real(prec) :: nx,ny,nz,nt,cta,cgm,cgm1,rrl,rrl_1,ccgm,nn
!    real(prec) :: c2,c,dcta,cta1,gama2,deti,det2,rodctac
!    real(prec) :: rodcta,dp_c2,a1,a2,a3,a4,a5,a6,a7
!    real(prec) :: r1,u1,v1,p1,r2,u2,v2,p2,n1,n2,flr(4),f1(4),g1(4),f2(4),g2(4),conser1(4),conser2(4)
!    nx=norm(1)
!    ny=norm(2)
!!    nn = sqrt(nx*nx+ny*ny)
!!    nx=norm(1)/nn
!!    ny=norm(2)/nn
!
!
!    call get_origi_var(conser1,ql) 
!    call get_origi_var(conser2,qr) 
!    r1=ql(1)
!    u1=ql(2)   
!    v1=ql(3)
!    p1=ql(4)
!   
!    r2=qr(1)
!    u2=qr(2)   
!    v2=qr(3) 
!    p2=qr(4)
!    call flux_fgh(r1,u1,v1,p1,f1,g1) 
!    call flux_fgh(r2,u2,v2,p2,f2,g2)  
!    do i=1,4
!    flr(i)=(f1(i)+f2(i))*nx+(g1(i)+g2(i))*ny
!    end do
!
!!   the difference of variables between right and left sides
!    do m=1,4
!        dq(m) = qr(m) - ql(m)
!    enddo
!
!!   to calculate H at left and right sides
!    gama2 = gamma/gamma1
!    ql(4) = gama2*ql(4)/ql(1) + 0.5_prec*(ql(2)*ql(2)+ql(3)*ql(3))
!    qr(4) = gama2*qr(4)/qr(1) + 0.5_prec*(qr(2)*qr(2)+qr(3)*qr(3))
!
!!   to calculate density of Roe average
!    qroe(1)=sqrt(ql(1)*qr(1))
!
!!   to calculate velocity and H of Roe average
!    rrl = sqrt(qr(1)/ql(1))
!    rrl_1 = rrl+1
!    do m=2,4
!        qroe(m)=(ql(m)+qr(m)*rrl)/rrl_1
!    enddo
!
!    rm = qroe(1)
!    um = qroe(2)
!    vm = qroe(3)
!    hm = qroe(4)
!
!!   to calculate the speed of sound
!    c2 = (hm-0.5_prec*(um*um+vm*vm))*gamma1
!    if(c2.LT.0.0_prec)then
!        write(*,*)'sqrt(c2) is a math error in subroutine Roe_FDS'
!        stop
!    endif
!    c = sqrt(c2)
!
!    cta = um*nx+vm*ny+nt
!
!    cgm = max(sqrt(nx*nx + ny*ny ),1.0e-30_prec)
!    cgm1= 1.0_prec/cgm
!
!    nx = nx*cgm1
!    ny = ny*cgm1
!
!    ccgm = c*cgm
!
!    l1 = abs(cta-ccgm)
!    l3 = abs(cta      )
!    l4 = abs(cta+ccgm)
!
!!   Harten's entropy modification
!    l12   = l1*l1
!    l32   = l3*l3
!    l42   = l4*l4
!
!    deti  = efix*cgm  
!    det2  = deti*deti  
!
!    l1 = sqrt(l12+det2)
!    l3 = sqrt(l32+det2)
!    l4 = sqrt(l42+det2)
!
!    dcta= dq(2)*nx+dq(3)*ny
!
!    cta1= cta*cgm1
!
!    rodcta = qroe(1)*dcta
!    dp_c2  = dq(4)/c2
!
!    rodctac = rodcta/c
!
!    a1 = l3*(dq(1)-dp_c2)
!    a2 = l4*(dp_c2+rodctac)*0.5_prec
!    a3 = l1*(dp_c2-rodctac)*0.5_prec
!    a4 = a1+a2+a3
!    a5 = c*(a2-a3)
!    a6 = l3*(rm*dq(2)-nx*rodcta)
!    a7 = l3*(rm*dq(3)-ny*rodcta)
!
!
!    df(1) = a4
!    df(2) = um*a4+nx*a5+a6
!    df(3) = vm*a4+ny*a5+a7
!    df(4) = hm*a4+cta1*a5+um*a6+vm*a7-c2*a1/gamma1
!
!
!
!
!    do m=1,4
!        ff(m) = 0.5_prec*(flr(m)-df(m))
!    enddo
!
!!    ff(:) = ff(:)*nn
!
!
!end subroutine Roe_FDS
!
!
!
!
!subroutine flux_fgh(r,u,v,p,f,g)
!use real_precision
!use global_variables,only: gamma,gamma1
!implicit none
!real(prec) ::r,u,v,p,f(4),g(4)
!real(prec) ::va
!va=u*u+v*v
!f(1)=r*u
!f(2)=r*u*u+p
!f(3)=r*u*v
!f(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*u
!
!g(1)=r*v
!g(2)=r*u*v
!g(3)=r*v*v+p
!g(4)=(p*(gamma/gamma1)+0.5_prec*r*va)*v
!
!end subroutine flux_fgh
!
!subroutine get_conser_var(ruvwp,conser)
!use real_precision
!use global_variables,only: gamma,gamma1
!implicit none
!real(prec) ::r,u,v,p
!real(prec) ::ruvwp(4),conser(4)
!r=ruvwp(1)
!u=ruvwp(2)
!v=ruvwp(3)
!p=ruvwp(4)
!conser(1)=r
!conser(2)=r*u
!conser(3)=r*v
!conser(4)=p/gamma1+0.5_prec*r*(u*u+v*v) 
!
!end subroutine get_conser_var
!
!subroutine get_origi_var(conser,ruvwp)
!use real_precision
!use global_variables,only: gamma1
!real(prec) ::ruvwp(4),conser(4)
!ruvwp(1)=conser(1)
!ruvwp(2)=conser(2)/conser(1)
!ruvwp(3)=conser(3)/conser(1)
!ruvwp(4)=gamma1*(conser(4)-0.5_prec*(ruvwp(2)*ruvwp(2)+ruvwp(3)*ruvwp(3))*ruvwp(1))
!end subroutine get_origi_var
!
! 
!
!!!Osherflux
!subroutine osherflux(conser1,conser2,norm,flr)
!use real_precision
!use global_variables,only: gamma,gamma1,gama2
!implicit none
!real(prec),intrinsic::abs
!real(prec),intrinsic::max     
!real(prec),intrinsic::sqrt  
!real(prec),external::sign_fun
!real(prec) ::r1,u1,v1,p1,r2,u2,v2,p2,flr(4),var1(4),var2(4),norm(2),conser1(4),conser2(4)
!real(prec) ::n1,n2,nn
!real(prec) ::c1,lu1,lv1,c2,lu2,lv2
!real(prec) ::mc1,mr1,mu1,mv1,mp1,mlu1,mlv1   
!real(prec) ::mc2,mr2,mu2,mv2,mp2,mlu2,mlv2
!real(prec) ::sc1,sr1,su1,sv1,sp1,slu1,slv1
!real(prec) ::sc2,sr2,su2,sv2,sp2,slu2,slv2
!real(prec) ::sign1,sign2,sign3,sign4,sign5,sign6
!real(prec) ::fai(4),fai13,f1(4),f2(4),fm1(4),fm2(4),fs1(4),fs2(4)
!real(prec) ::sign_integral 
!real(prec) ::lamda2m1,lamda11,lamda42,lamda1m1,lamda4m2
!real(prec) ::maxdiff,ma1,ma2
!real(prec) ::f10(4),g10(4),f20(4),g20(4)
!integer :: j
!integer :: m,bdkind0,road
!integer :: integral_order,regu_or_irregu,calculation_method,direct
!
! bdkind0=0
!
!call get_origi_var(conser1,var1) 
!call get_origi_var(conser2,var2) 
!!基于原始变量的转换
!r1=var1(1)
!u1=var1(2)   
!v1=var1(3)
!p1=var1(4)
!   
!r2=var2(1)
!u2=var2(2)   
!v2=var2(3)
!p2=var2(4)
!
!n1=norm(1)
!n2=norm(2) 
!nn = sqrt(N1*N1+N2*N2)
!n1=norm(1)/nn
!n2=norm(2)/nn 
!
!
!integral_order=1 !1表示P升序；2表示O降序
!calculation_method=2 !求中间状态的方法：2表示先求密度和压强，3表示先求逆变速度和声速
!regu_or_irregu=2    !1表示固壁通过添加对称虚拟单元进行处理；2表示根据固壁u=0作特殊处理；3表示直接求得固壁边的通量
!
!ma1=((u1**2.0_prec+v1**2.0_prec)/(gamma*p1/r1))**0.5_prec
!ma2=((u2**2.0_prec+v2**2.0_prec)/(gamma*p2/r2))**0.5_prec
!
!
!    lu1 = u1*n1 + v1*n2 
!    lv1 = v1*n1 - u1*n2
!    
!    lu2 = u2*n1 + v2*n2 
!    lv2 = v2*n1 - u2*n2
!
!
!c1=sqrt(gamma*p1/r1)
!c2=sqrt(gamma*p2/r2)
!
!
!if(integral_order==2)then
!sign_integral=-1.0_prec
!else
!sign_integral=1.0_prec
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!P序!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!
!    if(calculation_method==1)then
!        !中间状态的第一种计算方法
!        fai(1) = 0.5_prec*gamma1*(lu2-lu1)-sign_integral*(c2+c1)
!        fai(2) = ((p1/p2)**(0.5_prec/gamma))*sqrt(r2/r1)
!        mr2  = ((-sign_integral*fai(1)/(c2*(1.0_prec+fai(2))))**(2.0/gamma1))*r2
!        mr1  = ((-sign_integral*fai(1)/(c1*(1.0_prec+1.0_prec/fai(2))))**(2.0_prec/gamma1))*r1 
!        mp2  = p1*((r1/mr1)**(-gamma))
!        mp1  = mp2
!        mc2  = sqrt(gamma*mp2/mr2)
!        mc1  = sqrt(gamma*mp1/mr1)
!        mlu2 = lu1+sign_integral*(2.0_prec/gamma1)*(c1-mc1)
!      
!        mlv2 = lv2
!        mlu1 = mlu2
!        mlv1 = lv1
!   else if(calculation_method==2)then
!        !中间状态的第二种计算方法  
!        fai(1) = lu2-lu1
!        fai(2) = c2+c1
!        fai(3) = p1/p2
!        fai(4) = r1/r2
!    
!        mr2  = (   abs( -sign_integral*(0.5_prec*gamma1*fai(1)-sign_integral*fai(2))/(c2*(1.0_prec+(fai(3)**(0.5/gamma))*((fai(4))**(-0.5_prec)))) )**(2.0_prec/gamma1)   )*r2
!        mr1  = (   abs( -sign_integral*(0.5_prec*gamma1*fai(1)-sign_integral*fai(2))/(c1*(1.0_prec+(fai(3)**(-0.5/gamma))*((fai(4))**(0.5_prec)))) )**(2.0_prec/gamma1)   )*r1
!        mp2  = p1*((r1/mr1)**(-gamma))
!        mp1  = mp2
!        mc2  = sqrt(gamma*mp2/mr2)
!        mc1  = sqrt(gamma*mp1/mr1)
!        mlu2 = lu1+sign_integral*(2.0_prec/gamma1)*(c1-mc1)
!             
!        mlv2 = lv2
!        mlu1 = mlu2
!        mlv1 = lv1
!   else 
!        !中间状态的第三种计算方法
!        fai(1) = p2/(r2**gamma)
!        fai(2) = lu2-sign_integral*(2.0_prec/gamma1)*c2
!        fai(3) = p1/(r1**gamma)
!        fai(4) = lu1+sign_integral*(2.0_prec/gamma1)*c1
!        fai13=1.0_prec+(fai(3)/fai(1))**(0.5_prec/gamma)
!    
!        mlu2=fai(2)+(fai(4)-fai(2))/fai13
!        mc2=0.5_prec*gamma1*sign_integral*(mlu2-fai(2))
!        mr2=((mc2*mc2)/(gamma*fai(1)))**(1.0_prec/gamma1)
!        mp2=((mc2*mc2)*mr2)/gamma
!    
!        mlu1=mlu2     
!        mc1=0.5_prec*gamma1*sign_integral*(fai(4)-mlu1)
!        mr1=((mc1*mc1)/(gamma*fai(3)))**(1.0_prec/gamma1)
!    
!        mp1=mp2
!        mlv2 = lv2 
!        mlv1 = lv1
!  end if
!
!    call inverse_trans(mlu1,mlv1,n1,n2,mu1,mv1)
!    call inverse_trans(mlu2,mlv2,n1,n2,mu2,mv2)
!
!
!     lamda11=lu1-sign_integral*c1  !lamda1在左端点取值
!     lamda1m1=mlu1-sign_integral*mc1
!     lamda42=lu2+sign_integral*c2  !lamda5在右端点取值
!     lamda4m2=mlu2+sign_integral*mc2
!     lamda2m1=mlu1
!
!    
! 
!      if(sign_fun(lamda11)*sign_fun(lamda1m1)<0.0_prec)then
!      slu1=(gamma1/gama2)*lu1+sign_integral*(2.0_prec/gama2)*c1
!      sc1=sign_integral*slu1
!      sr1=(((sc1/c1)*(sc1/c1))**(1.0_prec/gamma1))*r1
!      sp1=p1*((sr1/r1)**gamma)
!      slv1=lv1
!        call inverse_trans(slu1,slv1,n1,n2,su1,sv1)
!        call Comp_f(sr1,su1,sv1,sp1,slu1,n1,n2,fs1)   
!     
!        else
!        fs1=0.0_prec            
!      end if
!      
!        if(sign_fun(lamda4m2)*sign_fun(lamda42)<0.0_prec)then
!        slu2=(gamma1/gama2)*lu2-sign_integral*(2.0_prec/gama2)*c2
!        sc2=-sign_integral*slu2
!        sr2= (((sc2/c2)*(sc2/c2))**(1.0_prec/gamma1))*r2
!        sp2=p2*((sr2/r2)**gamma)
!        slv2=lv2
!        call inverse_trans(slu2,slv2,n1,n2,su2,sv2)
!        call Comp_f(sr2,su2,sv2,sp2,slu2,n1,n2,fs2)   
!
!        else 
!        fs2=0.0_prec        
!    end if
!      
!call Comp_f(r1,u1,v1,p1,lu1,n1,n2,f1) 
!call Comp_f(r2,u2,v2,p2,lu2,n1,n2,f2) 
!call Comp_f(mr1,mu1,mv1,mp1,mlu1,n1,n2,fm1) 
!call Comp_f(mr2,mu2,mv2,mp2,mlu2,n1,n2,fm2) 
!
!  if(bdkind0==2.and.regu_or_irregu==2)then
!        if(road==1)then
!        !u+c路径
!              if(lu1+c1<=0.0)then
!                    slu1=-(gamma1/gama2)*(2.0*c1/gamma1-u1)
!                    sc1=-slu1
!                    slv1=lv1
!                    sr1=(((sc1/c1)*(sc1/c1))**(1.0/gamma1))*r1
!                    sp1=p1*((sr1/r1)**gamma)    !(slu1*slu1)*sr1/gamma
!                    call inverse_trans(slu1,slv1,n1,n2,su1,sv1)
!                    call Comp_f(sr1,su1,sv1,sp1,slu1,n1,n2,fs1) 
!              end if
!        else
!        !u-c路径    
!                if(lu1-c1>0.0)then
!                    slu1=(gamma1/gama2)*(2.0*c1/gamma1+u1)
!                    sc1=slu1
!                    slv1=lv1
!                    sr1=(((sc1/c1)*(sc1/c1))**(1.0/gamma1))*r1
!                    sp1=p1*((sr1/r1)**gamma)    !(slu1*slu1)*sr1/gamma
!                    call inverse_trans(slu1,slv1,n1,n2,su1,sv1)
!                    call Comp_f(sr1,su1,sv1,sp1,slu1,n1,n2,fs1)   
!                end if  
!         end if           
!  end if
!  
!    sign1 =            1.0_prec     + sign_fun(lamda11)
!    sign2 = sign_fun(lamda1m1) - sign_fun(lamda11)
!    sign3 = sign_fun(lamda2m1) - sign_fun(lamda1m1)
!    sign4 = sign_fun(lamda4m2) - sign_fun(lamda2m1)
!    sign5 = sign_fun(lamda42)  - sign_fun(lamda4m2)
!    sign6 =              1.0_prec   - sign_fun(lamda42)
!    flr=0.0_prec
!
!    if(sign1/=0.0_prec)then 
!      do m=1,4
!         flr(m)= flr(m) + 0.5_prec*sign1*f1(m)
!      end do
!    end if
!    
!    if(sign2/=0.0_prec)then 
!      do m=1,4
!         flr(m)= flr(m) + 0.5_prec*sign2*fs1(m)
!      end do
!!      write(*,*)222
!    end if
!     
!    if(sign3/=0.0_prec)then 
!      do m=1,4
!         flr(m)= flr(m) + 0.5_prec*sign3*fm1(m)
!      end do
!    end if  
! 
!      
!    if(sign4/=0.0_prec)then 
!      do m=1,4
!         flr(m)= flr(m) + 0.5_prec*sign4*fm2(m)
!      end do
!    end if
!       
!    if(sign5/=0.0_prec)then 
!      do m=1,4
!         flr(m)= flr(m) + 0.5_prec*sign5*fs2(m)
!      end do
!    end if
! 
!        
!    if(sign6/=0.0_prec)then 
!      do m=1,4
!         flr(m)= flr(m) + 0.5_prec*sign6*f2(m)
!      end do
!    end if
!
!      if(bdkind0==2.and.regu_or_irregu==2)then
!        if(road==1)then
!        !u+c路径
!!              if(lamda51>0)then
!!                flr=f1        !fm1  !f1-fm1  ！0.0,ok
!!              else
!!                flr=fm1         !fm1+f1-fs1      !fs1-fm1  !fs1-fm1,ok 
!!              end if
!
!              if(lu1+c1>0)then
!                flr=f1    !f1+fm1-fs1 ok  ！fm1-fs1 
!              else
!                flr=fs1        !fm1 ok        !fm1-f1  
!              end if
!
!         else 
!         !u-c路径
!              if(lu1-c1>0)then
!                flr=f1+fm1-fs1    !f1+fm1-fs1 ok  ！fm1-fs1 
!              else
!                flr=fm1          !fm1 ok        !fm1-f1  
!              end if
!
!         end if
!  end if
!
!  if(bdkind0==2.and.regu_or_irregu==2.and.direct==1)then
!     flr(1)=0.0
!     flr(2)=mp1*n1
!     flr(3)=mp1*n2
!     flr(4)=0.0
!  end if
!
!
!!!entropy fix
!!!call get_conser_var(var1,conser1) 
!!!call get_conser_var(var2,conser2) 
!!!    flr=flr-maxdiff*(conser2-conser1)
!!!
!
!flr=flr*nn
!end subroutine osherflux