module Lobatto_P_W
    
    !-----------------------------------------------------------------------------
    !
    !   方法：Gauss-Lobatto积分高斯点及权重的求解模块
    !   描述：可以求解任意阶Lobatto点和权值，对实现任意高阶CPR 方法铺垫。需要配合修正函数自动生成函数，还未添加。 
    !   参考：https://keisan.casio.com/exec/system/1280801905
    !   作者：gqShi 
    !   历史：修正函数还未实现任意阶求解 //2021.12.14
    !
    !-----------------------------------------------------------------------------    
    
    use real_precision
    use parameter_setting,only:nsp
     
    implicit none 
    
    integer :: n
    real(prec_dd),allocatable :: Lobatto_Point(:),Lobatto_Weight(:)
    real(ddouble_prec) :: eps
    
    contains
    
    subroutine get_pre
        implicit none
        if(prec_dd == ddouble_prec)then
            eps = 1E-34
        elseif(prec_dd == double_prec )then
            eps =  1E-15
        elseif(prec_dd == single_prec )then
            eps = 1E-8   
        end if
        
        n = nsp+1                                             !设置求解高斯点的个数
        allocate(Lobatto_Point(n),Lobatto_Weight(n))
    end subroutine get_pre
    
    function f_Legendre_n(x,p)                                   !定义 p-1 阶精度的函数f_Legendre_n(x)
        implicit none
        integer::i,p
        real(prec_dd) :: a(p),x,f_Legendre_n                 !a(p)代表p阶勒让德多项式
        a(1) = x                                    !1阶勒让德多项式
        a(2) = 1.5_prec*(x**2)-0.5_prec                 !2阶勒让德多项式
        if(p>=3)then
            do i = 3,p                                  !利用递推关系产生p阶勒让德多项式
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
            end do            
        end if
    
        f_Legendre_n = a(p)                                    !生成的 p 阶勒让德多项式 
    end function

    function f_Legendre_derivative(x,p)!生成n阶勒让德多项式的导数表达式
        implicit none
        integer :: i,p
        real(prec_dd) :: a(p), x, f_Legendre_derivative
        a(1) = x
        a(2) = 1.5_prec*x**2-0.5_prec

        if(p>=3)then
            do i = 3,p                                  !利用递推关系产生p阶勒让德多项式
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
            end do            
        end if
        f_Legendre_derivative = p*a(p-1)/(1-x**2)-p*x*a(p)/(1-x**2)         
    end function

    function f_dichotomy(a,b,p)                   !二分法求解函数的解
        implicit none
        real(prec_dd)::a,b,c,f_dichotomy                          
        integer :: p
        
        do while(.true.)
            c = (a+b)/2.0_prec
            if(f_Legendre_derivative(c,p-1)*f_Legendre_derivative(a,p-1) < 0)then
                b = c
            elseif(f_Legendre_derivative(c,p-1)*f_Legendre_derivative(a,p-1) > 0)then
                a = c
            else
                exit
            endif
            if((b-a) < eps)exit
            
        enddo
        f_dichotomy = c         
    end function
    
    subroutine over_post
        implicit none
        deallocate(Lobatto_Point,Lobatto_Weight)
    end subroutine over_post
        
    subroutine get_Lobatto_Point_Weight(Lobatto_Point,Lobatto_Weight)
        implicit none
        real(prec_dd) :: m,delta,Lobatto_Point(:),Lobatto_Weight(:)    !已经分配            
        integer :: i,j

        call get_pre

        j = 1                                                                                 
        Lobatto_Point(j) = -1.0_prec
        Lobatto_Weight(j) =  2.0_prec/(n*(n-1)*(f_Legendre_n(Lobatto_Point(j),n-1))**2)  
        m = -1.0001_prec                                                                                        
        do i = 1,20000                                                                                          
            if(f_Legendre_derivative(m,n-1)*f_Legendre_derivative(m+0.0001_prec ,n-1)<0)then               
                j = j+1                                                                                 
                Lobatto_Point(j) = f_dichotomy(m,m+0.0001_prec ,n)                                         
                                        
                Lobatto_Weight(j) =  2.0_prec/(n*(n-1)*(f_Legendre_n(Lobatto_Point(j),n-1))**2)            
            endif
             m = m+0.0001_prec                                 
        enddo
        j = j+1 
        Lobatto_Point(j) = 1.0_prec
        Lobatto_Weight(j) =  2.0_prec/(n*(n-1)*(f_Legendre_n(Lobatto_Point(j),n-1))**2)
        !write(*,*) ' '
        !
        !do j = 1,n         
        !    write(*,*)  'Lobatto点',Lobatto_Point(j)
        !end do
        !write(*,*) ' '
        !do j = 1,n         
        !    write(*,*)  'Lobatto权',Lobatto_Weight(j)
        !end do
        !

        call over_post
    end subroutine
end module Lobatto_P_W 