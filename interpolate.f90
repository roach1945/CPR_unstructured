!!!!!!!!Lagrange interpolate
function LaI_nPs(x,xk,func,n)
    !!n阶lagrange 
    use real_precision
    use global_var ,only:nsp
    implicit none
    integer :: n
    real(prec) :: LaI_nPs,x,temp,temp2
    real(prec),dimension(n) :: xk,func
    integer i,j
    LaI_nPs = 0.0_prec
    do i =1,n
        temp2 = 1.0_prec
        do j = 1,n
            if(i .NE. j)then
                temp = (x-xk(j))/(xk(i)-xk(j))
                temp2 = temp2*temp
            end if
        end do
        LaI_nPs = temp2*func(i)+LaI_nPs
    end do
    return 
end function LaI_nPs

subroutine LaI_nPs_arr(x,xk,func,n,res)
    !!n阶lagrange 
    use real_precision
    use global_var,only:nsp
    implicit none
    integer :: n
    real(prec) :: x,temp,temp2,xk(n),func(n,4),res(4)
    !real(prec),dimension(n) :: xk
    !real(prec),dimension(n,4) :: func
    !real(prec),dimension(4) :: LaI_nPs_array
    integer i,j
    !write(*,*) n,x!,xk,func,n
    res = 0.0_prec
    do i = 1,n
        temp2 = 1.0_prec
        do j = 1,n
            if(i .NE. j)then
                temp = (x-xk(j))/(xk(i)-xk(j))
                temp2 = temp2*temp              
            end if
        end do
        res(:)= temp2*func(i,:)+res(:)
    end do
    return 
end subroutine LaI_nPs_arr

function LaI_nPs_2D(x,xk,func,n)
    !2d Larange插值
    !用不到，先放着，和三角形上的好像也不一样。
    use real_precision
    use global_var ,only:nsp
    implicit none
    integer n
    real(prec) :: LaI_nPs_2D,x,temp,temp2
    real(prec),dimension(n) :: xk,func
    integer i,j
    
    
    
    LaI_nPs_2D = 0.0_prec
    
    return
end function LaI_nPs_2D

function LaI_nPs_deri(x,xk,func,n)
    !!n points Lagrange interpolation derivation n点拉格朗日插值多项式导数,用来求f_sub_kesi
    use real_precision 
    use global_var ,only:nsp    
    implicit none
    real(prec) :: LaI_nPs_deri,x,temp,numerator,numerator1,denominator
    real(prec),dimension(n) :: xk,func
    integer i,j,k,n
    LaI_nPs_deri = 0.0_prec
 
    do i =1,n
        temp = 1.0_prec    
        numerator = 0.0_prec
        denominator = 1.0_prec
        do j = 1,n
            numerator1 = 1.0_prec
            if(i .NE. j)then
                do k = 1,n
                    if(k .NE. i .AND. k .NE. j)then
                        numerator1 = numerator1*(x - xk(k))
                    end if                    
                end do
                numerator = numerator + numerator1
                denominator = denominator*(xk(i)-xk(j))
            end if
        end do

        temp = numerator/denominator
        LaI_nPs_deri = temp*func(i)+LaI_nPs_deri
    end do
    return    
    
end function LaI_nPs_deri

subroutine LaI_nPs_deri_arr(x,xk,func,n,res)
    !!n points Lagrange interpolation derivation n点拉格朗日插值多项式导数,用来求f_sub_kesi
    use real_precision 
    use global_var ,only:nsp    
    implicit none
    real(prec) :: res(4),x,temp,numerator,numerator1,denominator,xk(n),func(n,4)
    integer i,j,k,n
    res = 0.0_prec
 
    do i =1,n
        temp = 1.0_prec    
        numerator = 0.0_prec
        denominator = 1.0_prec
        do j = 1,n
            numerator1 = 1.0_prec
            if(i .NE. j)then
                do k = 1,n
                    if(k .NE. i .AND. k .NE. j)then
                        numerator1 = numerator1*(x - xk(k))
                    end if                    
                end do
                numerator = numerator + numerator1
                denominator = denominator*(xk(i)-xk(j))
            end if
        end do

        temp = numerator/denominator
        res = temp*func(i,:)+res
    end do
    !write(*,*)res
    return    
    
end subroutine LaI_nPs_deri_arr
function linear_inter(x,xk,func,n)
!线性插值
    use real_precision
    implicit none
    integer :: n,i
    real(prec) :: linear_inter,x
    real(prec),dimension(n) :: xk,func
    do i = 2,n
        if(x>=xk(i-1) .AND. x<=xk(i))then
            linear_inter = (func(i)-func(i-1))/(xk(i)-xk(i-1))*(x-xk(i))+func(i)
            !write(*,*)i,xk(i-1),xk(i)
            !write(*,*)x,linear_inter 
            exit
        elseif(x<xk(1))then
            linear_inter = (func(2)-func(1))/(xk(2)-xk(1))*(x-xk(2))+func(2)
            exit
        elseif(x>xk(n))then
            linear_inter = (func(n)-func(n-1))/(xk(n)-xk(n-1))*(x-xk(n))+func(n)
            exit
        end if
        !write(*,*)'--' 
    end do
           
    return
end function linear_inter  

!---分段线性函数插值-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function Linear_nPs(x,xk,func,n)
    !分段一次线性函数插值
    !n指数组里真实存在的数据
    use real_precision 
    implicit none
    integer :: n,i
    real(prec) :: x
    real(prec) :: Linear_nPs
    real(prec),dimension(n) :: xk,func 

    do i = 2,n
        if(x .GE. xk(i-1) .AND. x .LT. xk(i))then
            Linear_nPs = func(i-1)+(func(i)-func(i-1))/(xk(i)-xk(i-1))*(x-xk(i-1))
            !write(*,*) i,xk(i-1),xk(i),func(i-1),func(i),Linear_nPs
            exit
        end if
    end do       
    return
end function Linear_nPs
!---分段线性函数插值 End-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function LaI_nPs_Bary(x,xk,func,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Barycentric Lagrange interpolate 
!!!!!!!!!! 计算次数下降
    use real_precision
    implicit none
    integer n
    real(prec) :: LaI_nPs_Bary,x,temp,temp1,temp2
    real(prec),dimension(10) :: xk,func,omega
    integer i,j
    LaI_nPs_Bary = 0.0_prec
    temp1 = 0.0_prec
    temp2 = 0.0_prec
    omega = 1.0_prec
    do i =1,n
        do j = 1,n
            if(i .NE. j)then
                temp = 1.0_prec/(xk(i)-xk(j))
                omega(i) = omega(i)*temp
            end if
        end do
    end do

    do j = 1,n
        temp1 = temp1 + omega(j)*func(j)/(x-xk(j))
        temp2 = temp2 + omega(j)/(x-xk(j))
    end do
        
    LaI_nPs_Bary = temp1/temp2

    return 
end function LaI_nPs_Bary

!---非线性权WCNS插值-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function NonLinear_WCNS(x,xk,un)
    !非线性权WCNS插值
    !n指数组里真实存在的数据
    use real_precision 
    implicit none
    integer :: n,i,k
    real(prec) :: x
    real(prec) :: NonLinear_WCNS 
    real(prec),dimension(:) ::d(2),beta(2),IS(2),xk(-1:1),un(-1:1),omega(0:1)
    real(prec)::eps,beta_sum,tao,u1,u2,ux
    real(prec)::LaI_nPs
    real(prec) :: A(2,2),B(2)
    eps = 0.000001_prec
    d = 0.0_prec!线型系数
    u1 = LaI_nPs(x,xk(-1:0),un(-1:0),2)
    u2 = LaI_nPs(x,xk(0:1),un(0:1),2)
    ux = LaI_nPs(x,xk(-1:1),un(-1:1),3)
    A(1,:) = 1.0_prec
    A(2,1) = u1
    A(2,2) = u2
    B(1) = 1.0_prec
    B(2) = ux
    !write(*,*)A(1,:),A(2,:)
    call Gaussh(A,B,d,2)
    !write(*,*)'d',d

    !IS(1) = 0.25_prec*(3.0_prec*un(0)-4.0_prec*un(-1)+un(-2))**2 + (un(0)-2.0_prec*un(-1)+un(-2))**2
    !IS(2) = 0.25_prec*(un(2)-4.0_prec*un(1)+3.0_prec*un(0))**2 + (un(2)-2.0_prec*un(1)+un(0))**2
    
    IS(1)=(un(-1)-un(0))**2
    IS(2)=(un(0)-un(1))**2
    tao = abs(IS(1)-IS(2))

    beta_sum = 0.0_prec
    do k = 1,2
        beta(k) = d(k)*(4.0_prec+(tao/(eps+IS(k)))**2)
        beta_sum = beta_sum + beta(k)
    end do

    do k = 1,2
        omega(k-1) = beta(k)/beta_sum
    end do
   
    NonLinear_WCNS = omega(0)*u1 + omega(1)*u2
    !NonLinear_WCNS = d(1)*u1 + d(2)*u2
    return
end function NonLinear_WCNS
!---非线性权WCNS插值 End-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


function LSI_nPmOr(x,xk,func,m,n,ns)
!!!!!!!least square interpolate, n points, m order最小二乘多项式插值
    use real_precision

    implicit none
    integer i,j,k,m,n,power,ns
    real(prec) :: x,LSI_nPmOr
    real(prec),dimension(:) :: xk(ns),func(ns),B(m+1),C(m+1)
    real(prec),dimension(:,:) :: A(m+1,m+1)!m表示多项式阶数，m+1表示矩阵维度
    A = 0.0_prec
    C = 0.0_prec
    B = 0.0_prec
    LSI_nPmOr = 0.0_prec

!!!!!!!!!!  A    
    do i = 1,m+1
        do j = 1,m+1
            do k = 1,n
                power = 2*m+2-i-j
                A(i,j) = A(i,j)+xk(k)**power
            end do
            !write(*,*) i,j,A(i,j)
        end do
    end do
!!!!!!!!!!  B  
    do i = 1,m+1
        do k = 1,n
            power = m+1-i
            B(i) = B(i)+func(k)*(xk(k)**power)
            !write(*,*) 'B',i,B(i)
        end do
        !write(*,*) 'power',power
        !write(*,*) 'xk',xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7)
        !write(*,*) 'func',func(1),func(2),func(3),func(4),func(5),func(6),func(7)
        !write(*,*) 'func-',xk(7)*(func(7)-func(1)),xk(6)*func(6)-func(2),xk(5)*(5)-func(3),xk(4)*func(4)
        !write(*,*) 'func+',xk(7)*(func(7)-func(1))+xk(6)*func(6)-func(2)+xk(5)*(5)-func(3)+xk(4)*func(4)
        !write(*,*) 'B',i,B(i)
    end do
    call Gaussh(A,B,C,m+1)
    do i = 1,m+1
        power = m+1-i
        LSI_nPmOr = LSI_nPmOr + C(i)*(x**power)
    end do
    !write(*,*)C(1),C(2),C(3)
return

end function

subroutine Gaussh(A,B,C,n)
    use real_precision
    implicit none
    integer i,j,k,m,n
    real(prec),dimension(:,:) :: Ab(n,n+1),A(n,n)
    real(prec),dimension(:) :: B(n),C(n)
    real(prec) :: temp,bb,cc
    Ab(1:n,1:n) = A
    Ab(1:n,n+1) = B
 
    do k = 1,n-1
		! 检查Ab[k][k]是否为0，为0则与k+1行交换
		if(abs(Ab(k,k)) < 0.1e-15)then
            do m = 1,n+1 !有n行，每一行有n+1个元素 ，即Ab=n*(n+1)
				temp = Ab(k,m)
				Ab(k,m)=Ab(k+1,m)
				Ab(k+1,m)=temp	
			end do
		end if
	 !/*化为上三角矩阵 */
        do i = k,n-1
			bb=Ab(i+1,k)/Ab(k,k)
            do j = k,n+1
				Ab(i+1,j) = Ab(i+1,j)-bb*Ab(k,j)	
            end do
        end do	
	 	if(Ab(n,n)== 0.0)then
	 		!write(*,*) "没有唯一解！"
        end if
    end do

	!反向求解未知数
    if(Ab(n,n)== 0.0)then
        C(n) = 0.5_prec
    else
	    C(n) = Ab(n,n+1)/Ab(n,n)
    end if
    i = n-1
    do while(i>0)
		cc=0.0_prec

        do j = i+1,n
			cc=cc+Ab(i,j)*C(j)
        end do
		C(i)=(Ab(i,n+1)-cc)/Ab(i,i)
        i = i-1
	end do 

end subroutine Gaussh



function lim(Var,VarC,varMax,varMin)
    !!Ref: T.J. Barth, D.C. Jespersen, The design and application of upwind schemes on unstructured meshes, AIAA Paper 89-0366, 1989, pp. 1–12.
    use real_precision
    implicit none    
    real(prec) :: lim,Var,VarC,varMax,varMin,temp
    
    if(VarC > Var)then
        temp = (varMax-Var)/(VarC-Var)
        lim = min(1.0_prec,temp)
        !write(*,*)lim
    elseif(VarC < Var)then
        temp = (varMin-Var)/(VarC-Var)
        lim = min(1.0_prec,temp)
    elseif(Var == VarC)then
        lim = 1.0_prec
    endif
    return
end function lim

function vanLeerLimiter(d1,d2,u1,u2,u3)
    use real_precision
    
    implicit none
    real(prec) :: d1,d2,u1,u2,u3
    real(prec) :: du1,du2,theta,vanLeerLimiter
    
    du1 = (u2 - u1)/d1
    du2 = (u3 - u2)/d2

    if(du2 == 0)then
        theta = 1.0_prec
    else
        theta = du1/du2
    end if
    
    vanLeerLimiter = (abs(theta)+(theta))/(1.0_prec+abs(theta))
    !write(*,*)vanLeerLimiter
    return
end function vanLeerLimiter
!
!subroutine upwind_nonuniform_lr(num0,wcns_q,wcns_qlr)
!    use real_precision
!    use global_variables,only: trd,ebs,twone,dhw,stru_sp_kesi,stru_fp_kesi_w
!    implicit none
!    real(prec),intrinsic::sqrt  
!    integer :: i,j,m,s,num0
!    real(prec) :: wcns_q(1:5),wcns_qlr(1:2)
!     real(prec) :: u1(1:3) 
!      real(prec) :: maxu,minu,ua,ub,limua,limub,limu,du,ua0,ub0,dh1(0:6),omega1,omega2
!     real(prec) ::Ldx1,Ldx2,Rdx1,Rdx2,LRdx,du1,du2,du3,theta0,sig1,sig2
! 
!           !(7)vanleer-new
!    if(num0==1)then
!         Ldx1=stru_fp_kesi_w(num0) - (-2.0_prec+stru_sp_kesi(5))   
!      else
!         Ldx1=stru_fp_kesi_w(num0) - stru_sp_kesi(num0-1)
!       end if
!       Ldx2=stru_sp_kesi(num0) - stru_fp_kesi_w(num0)  
!      omega1=(1.0_prec/Ldx1)/(1.0_prec/Ldx1 + 1.0_prec/Ldx2)
!      omega2=(1.0_prec/Ldx2)/(1.0_prec/Ldx1 + 1.0_prec/Ldx2)
!      ua0=omega1*wcns_q(2)+omega2*wcns_q(3)
!      
!
!      Rdx1=stru_fp_kesi_w(num0+1) - stru_sp_kesi(num0)
!      if(num0==5)then
!        Rdx2=(2.0_prec+stru_sp_kesi(1)) - stru_fp_kesi_w(num0+1)    
!      else
!        Rdx2=stru_sp_kesi(num0+1) - stru_fp_kesi_w(num0+1)  
!      end if
!      omega1=(1.0_prec/Rdx1)/(1.0_prec/Rdx1 + 1.0_prec/Rdx2)
!      omega2=(1.0_prec/Rdx2)/(1.0_prec/Rdx1 + 1.0_prec/Rdx2)
!      ub0=omega1*wcns_q(3)+omega2*wcns_q(4)
!
!      du=(  (ub0-wcns_q(3))/(Rdx1*Rdx1)+ (wcns_q(3)-ua0)/(Ldx2*Ldx2)  )/( 1.0_prec/Rdx1 + 1.0_prec/Ldx2)
!      
!        du1=(wcns_q(3)-ua0)/Ldx2
!        du2=(ub0-wcns_q(3))/Rdx1  
!        
!        if(du2==0)then
!            theta0=1
!        else
!            theta0=du1/du2
!        end if
!        limu=(abs(theta0)+(theta0))/(1.0_prec+abs(theta0))
!  
!      wcns_qlr(1) =wcns_q(3) + limu*(stru_fp_kesi_w(num0)- stru_sp_kesi(num0))*du2
!      wcns_qlr(2) =wcns_q(3) + limu*(stru_fp_kesi_w(num0+1)- stru_sp_kesi(num0))*du2  
! 
!    end subroutine upwind_nonuniform_lr   


