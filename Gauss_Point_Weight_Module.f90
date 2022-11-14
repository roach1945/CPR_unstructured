module Gauss_P_W

    !-----------------------------------------------------------------------------
    !
    !   方法：高斯—勒让德积分高斯点及权重的求解模块
    !   描述：可以求解任意阶Gauss点和权值，对实现任意高阶CPR 方法铺垫。需要配合修正函数自动生成函数，还未添加。 
    !   参考：http://bbs.fcode.cn/thread-219-1-1.html 
    !   作者：gqShi 
    !   历史：修正函数还未实现任意阶求解 //2021.12.14
    !
    !-----------------------------------------------------------------------------    
    
    use real_precision
    use parameter_setting,only:nsp

    
    implicit none 
    integer :: n                        
    real(prec_dd),allocatable :: Gauss_Point(:),Gauss_Weight(:)
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
        
        n = nsp                                             !设置求解高斯点的个数
        allocate(Gauss_Point(n),Gauss_Weight(n))
    end subroutine get_pre
    
    subroutine over_post
        implicit none
        deallocate(Gauss_Point,Gauss_Weight)
    end subroutine over_post
    
    function f_Legendre_n(x)                                !定义 n-1 阶精度的函数f_Legendre_n(x)
        implicit none
        integer::i
        real(prec_dd) :: a(n),x,f_Legendre_n                !a(n)代表n阶勒让德多项式
        a(1) = x                                            !1阶勒让德多项式
        a(2) = 1.5_prec*(x**2)-0.5_prec                     !2阶勒让德多项式

        if(n>=3)then
            do i = 3,n
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
            enddo            
        end if
        f_Legendre_n = a(n)                                 !生成的 n 阶勒让德多项式 
    end function

    function f_Legendre_nMinusOne(x) 
        implicit none
        integer::i
        real(prec_dd)::a(n), x, f_Legendre_nMinusOne
        a(1) = x
        a(2) = 1.5_prec*x**2-0.5_prec
        do i = 3,n-1
            a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
        enddo
        f_Legendre_nMinusOne = a(n-1)                                   !生成的（n-1)阶勒让德多项式 
    end function

    function f_Legendre_derivative(x)
        implicit none
        integer :: i
        real(prec_dd) :: a(n), x, f_Legendre_derivative
        a(1) = x
        a(2) = 1.5_prec*x**2-0.5_prec
        if(n>=3)then
            do i = 3,n
                a(i) = (2*i-1)*x*a(i-1)/i-(i-1)*a(i-2)/i
                !write(*,*)a(i)
            enddo            
        end if
        f_Legendre_derivative = n*a(n-1)/(1-x**2)-n*x*a(n)/(1-x**2)     !生成n阶勒让德多项式的导数表达式.第五递推公式
    end function

    function f_dichotomy(a,b)                                           !二分法求解函数的解
        implicit none
        real(prec_dd)::a,b,c,f_dichotomy                                !a,b是传递进来的划分好的有一个解存在的区间
        do while(.true.)
            c = (a+b)/2.0_prec
            if(f_Legendre_n(c)*f_Legendre_n(a) < 0)then
                b = c
            else
                a = c
            endif
            if((b-a) < eps)exit
            
        enddo
        f_dichotomy = c                                                     !f_dichotomy即是利用二分法求得的解
    end function

    subroutine get_Gauss_Point_Weight(Gauss_Point,Gauss_Weight)
        implicit none
        real(prec_dd) :: m,Gauss_Point(n),Gauss_Weight(n)                   !定义数组,大小n由module开始声明。
        integer :: i,j
        call get_pre

        j = 0                                                              !赋值控制循环变量的初值 
        m = -1.00001                                                       !设置计算域[-1，1] 的下限，即代替-1 
        do i = 1,20000                                                     !这个循环次数应该是由步长0.00001决 定,计算方法：200000=2/0.000001 
            if(f_Legendre_n(m)*f_Legendre_n(m+0.0001)<0)then               !从下限处开始往上逐步累加，!由步长0.00001说明最多求解10^5个解
                j = j+1                                                    !记录这是第几个解
                Gauss_Point(j) = f_dichotomy(m,m+0.0001)                   !调用二分法求解程序在分好的一小段上求解，将解存储在Gauss_Point（j)                                        
                Gauss_Weight(j) = 2.0_prec/(n*f_Legendre_nMinusOne(Gauss_Point(j))*f_Legendre_derivative(Gauss_Point(j)))     !高斯点的权重
            endif
             m = m+0.0001                                
        enddo
        !do j = 1,n         
        !    write(*,*)  '高斯点',Gauss_Point(j)
        !end do
        !write(*,*) ' '
        !do j = 1,n         
        !    write(*,*)  '高斯权',Gauss_Weight(j)
        !end do
        call over_post

    end subroutine

end module Gauss_P_W