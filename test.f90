subroutine test
    !测试模块
    !验证子程序正确性的附加模块，可去除
    use real_precision
    use parameter_setting
    use type_module
    use test_module
    use global_var
    implicit none
    
    character(len = 20)::filename  ! 文件名称
    character(len = 30)::str1,str2,str3     
    real(prec) :: A(2,2),B(2,2)
    integer i,j
    type(cell_record) cell
    type(side_record) side
    real(prec),dimension(:,:) :: cell_vertex(4,2)   !实际单元顶点，标准单元已知
    real(prec),dimension(:)   :: coor_C(2),coor_P(2)!标准单元，实际单元
    real(prec),dimension(:)   :: M_Jacobi(4),M_direct(4)          !标准单元,Jacobi矩阵
    real(prec) :: detJ
    real(prec) :: gl,gr,ruvpL(4),ruvpR(4),rhop_roe
    real(prec),external :: gl_sub_kesi, gr_sub_kesi
    real(prec) :: rho_L,u_L,v_L,p_L  !初始时刻，分隔面左侧参数
    real(prec) :: rho_R,u_R,v_R,p_R  !初始时刻，分隔面右侧参数
    
    !测试矩阵求逆
    A(1,1)=1.0_prec
    A(1,2)=1.0_prec
    A(2,1)=1.0_prec
    A(2,2)=1.0_prec
    !call Matrix_Inverse(A,2,B)
    !write(*,*)'b',B
    
    rho_L = 1.0_prec
    u_L = 0.0_prec
    v_L = 0.0_prec
    p_L = 1.0_prec   !初始时刻，分隔面左侧参数,sod激波管     
    rho_R = 0.125_prec
    u_R = 0.0_prec
    v_R = 0.0_prec
    p_R = 0.1_prec 
    ruvpL(1) = 1.0_prec
    ruvpL(2) = 0.0_prec
    ruvpL(3) = 0.0_prec
    ruvpL(4) = 1.0_prec   !初始时刻，分隔面左侧参数,sod激波管     
    ruvpR(1) = 1.0_prec
    ruvpR(2) = 0.0_prec
    ruvpR(3) = 0.0_prec
    ruvpR(4) = 1.0_prec 
    !
    !ruvpL(1) = 1.0_prec
    !ruvpL(2) = 1.0_prec
    !ruvpL(3) = 1.0_prec
    !ruvpL(4) = 1.0_prec   !初始时刻，分隔面左侧参数,sod激波管     
    !ruvpR(1) = 1.0_prec
    !ruvpR(2) = 1.0_prec
    !ruvpR(3) = 1.0_prec
    !ruvpR(4) = 1.0_prec 
    call RoeAverage(ruvpL,ruvpR,rhop_roe)
    
    !cell = cell_record(1.0) 
    !cell.nodes = 1
    !side.nodes = 2
    !write(*,*)cell.nodes(1)
    !write(*,*)side.nodes
    !m = 10 
    !n = 4
    !call test_sub1
    !do i =1,m
    !    allocate(arrayset(i).array(n))
    !    arrayset(i).array = i
    !end do
    !write(*,*)arrayset(1).array
    cell_vertex(1,1) = 0.0_prec;    cell_vertex(1,2) = 0.0_prec
    cell_vertex(2,1) = 2.0_prec;    cell_vertex(2,2) = 0.0_prec
    cell_vertex(3,1) = 1.0_prec;    cell_vertex(3,2) = 2.0_prec
    cell_vertex(4,1) = 0.0_prec;    cell_vertex(4,2) = 1.0_prec
    coor_C(1) = 0.9_prec
    coor_C(2) = 0.8_prec
    
    !cell_vertex(4,1) = 0.0_prec;    cell_vertex(4,2) = 0.0_prec
    !cell_vertex(3,1) = 2.0_prec;    cell_vertex(3,2) = 0.0_prec
    !cell_vertex(2,1) = 1.0_prec;    cell_vertex(2,2) = -2.0_prec
    !cell_vertex(1,1) = 0.0_prec;    cell_vertex(1,2) = -1.0_prec
    !coor_C(1) = 0.9_prec
    !coor_C(2) = -0.8_prec
    !call quad_C2Phy(cell_vertex,coor_C,coor_P)
    !call solve_Jacobi(cell_vertex,coor_C,M_Jacobi,M_direct,detJ)
    !write(*,*) detJ
    !write(*,*) M_Jacobi
    !write(*,*) M_direct
    !call sort_qua(cell_vertex)
    !write(*,*) ''
    !do j =1,4
    !    write(*,"(2F10.5)")cell_vertex(j,:)
    !end do
    
    !do j = 1,nsp
    !    gl = gl_sub_kesi(SPs(j))
    !    gr = gr_sub_kesi(SPs(j))
    !    write(*,*) gl,gr
    !end do
end subroutine test


subroutine StringSplit(InStr,delimiter,StrArray,nsize)
!----------------------------------------------
!---将字符串InStr进行分割,结果放入StrArray中
!---delimiter::分隔符号,例如';,,' 使用;和,分割字符串
!---nsize:分割数目
!---吴徐平2011-04-29(wxp07@qq.com)
!----------------------------------------------
implicit none
character(len = *) , Intent( IN ) :: InStr
character(len = *)  , Intent( IN ) :: delimiter
character(len = LEN(InStr)),dimension(LEN(InStr)),Intent( OUT ) :: StrArray
integer, Intent( OUT ) :: nsize ! Effective Size of StrArray
integer:: i,j ! loop variable
integer:: istart ! split index for Start Position
nsize=0
istart=1

do i=1,LEN(InStr)
	do j=1,LEN(delimiter)
		if (InStr(i:i) == delimiter(j:j)) then
			if (istart == i) then
			istart=i+1 ! ---可防止分隔符相连的情况
			end if
			if (istart<i) then
				nsize=nsize+1
				StrArray(nsize)=InStr(istart:i-1)
				istart=i+1
			end if
		end if
	end do
end do
! ---匹配最后一个子字符串
if (nsize>0) then
	if (istart<LEN(InStr)) then
		nsize=nsize+1
		StrArray(nsize)=InStr(istart:LEN(InStr))
	end if
end if
! ---如果无可分割的子字符串,则包含整个字符串为数组的第一元素
if ( (nsize<1) .AND. (LEN(TRIM(InStr)) > 0 )) then
		nsize=1
		StrArray(1)=InStr
end if
end subroutine StringSplit
!
!=============================================================
subroutine StrReplace(InStr,OldChar,NewChar,OutStr)
!------------------------------------------------------------
!---将字符串InStr中的字符串OldChar替换成NewChar
!---结果放入字符串OutStr中
!---吴徐平2013-07-20(wxp07@qq.com)
!------------------------------------------------------------
implicit none
character(len = *) , Intent( IN ) :: InStr
character(len = *) , Intent( IN ) :: OldChar
character(len = LEN(OldChar)) , Intent( IN ) ::NewChar
character(len = LEN(InStr)) , Intent( INOUT ) :: OutStr
integer :: i  ! loop variable
OutStr=InStr
i=INDEX(OutStr,OldChar)
do while(i>0)
	OutStr(i:i+LEN(OldChar)-1)=NewChar
	i=INDEX(OutStr,OldChar)
end do
end subroutine StrReplace
!------------------------------------------------------------
