Module SparseMatrix
  Use typedef
  Implicit None

  Type SparseMatrixNode
     Integer :: colno
     Real(dp) :: Data
     Type(SparseMatrixNode), Pointer :: next
  End Type SparseMatrixNode
  Type SparseMatrixLine 
     Integer :: numnodes
     Type(SparseMatrixNode), Pointer :: first
  End Type SparseMatrixLine
  Type SparseMatrixType
     Integer :: nlines
     Integer :: ncols
     Type(SparseMatrixLine), Pointer :: lines(:)
  End Type SparseMatrixType

  Type(SparseMatrixNode), Pointer, Save :: pnew

Contains

  
  Function SM_allocate(matrix, nl, nc)
    Use typedef
    Implicit None
    Integer :: SM_allocate
    Type(SparseMatrixType), Intent(INOUT) :: matrix
    Type(SparseMatrixNode), Pointer :: pnode
    Integer, Intent(in) :: nl, nc
    Integer :: err, i

    SM_allocate=1
    matrix%nlines=nl
    matrix%ncols=nc
    Allocate(matrix%lines(nl), stat=err)
    If(err/=0) Return
    Do i=1,nl
       ! allocate the diagonal element at position 1
       Allocate(pnode, stat=err)
       If(err/=0) Return
       matrix%lines(i)%first=>pnode
       matrix%lines(i)%numnodes=1
       pnode%colno=i
       pnode%data=0.0
       Nullify(pnode%next)
    Enddo
    SM_allocate=0
    Nullify(pnew)
  End Function SM_allocate

  
  Subroutine SM_delete(matrix)
    Use typedef
    Implicit None
    Type(SparseMatrixType), Intent(INOUT) :: matrix
    Integer :: i
    Type(SparseMatrixNode), Pointer :: p1, p2

    Do i=1,matrix%nlines
       p1=>matrix%lines(i)%first
       Do While(Associated(p1))
          p2=>p1%next
          Deallocate(p1)
          p1=>p2
       Enddo
    Enddo
    Deallocate(matrix%lines)
    matrix%nlines=0
    matrix%ncols=0
    Nullify(matrix%lines)
  End Subroutine SM_delete

  
  Function SM_numberOfElements(matrix)
    Implicit None
    Integer :: SM_numberOfElements
    Type(SparseMatrixType), Intent(IN) :: matrix
    SM_numberOfElements=Sum(matrix%lines(:)%numnodes)
  End Function SM_numberOfElements

  
  Function SM_set(matrix, nl, nc, value, replace)
    Implicit None
    Integer :: SM_set
    Type(SparseMatrixType), Intent(INOUT) :: matrix
    Integer, Intent(IN) :: nl, nc
    Real(dp), Intent(IN) :: value
    Logical, Intent(IN) :: replace
    Integer ::  err
    Type(SparseMatrixNode), Pointer :: p1, p2

    ! set the pointers
    p1=>matrix%lines(nl)%first
    If(nl==nc) Then
       ! diagonal element is always the first element in the line 
       If(replace) Then
          p1%data=value
       Else
          p1%data=p1%data+value
       Endif
    Else 
       ! insert in list
       SM_set=1
       If(.Not.Associated(pnew)) Then
          Allocate(pnew, stat=err)
          If(err/=0) Return
       Endif
       ! set data
       pnew%colno=nc
       pnew%data=value
       Do
          p2=>p1%next
          If(.Not.Associated(p2)) Then
             ! insert at bottom
             p1%next=>pnew
             matrix%lines(nl)%numnodes=matrix%lines(nl)%numnodes+1
             Nullify(pnew%next)
             Nullify(pnew)
             Exit
          Else If(p2%colno==nc) Then
             ! node already exists
             If(replace) Then
                p2%data=value
             Else
                p2%data=p2%data+value
             Endif
             Exit
          Else If(p2%colno>nc) Then
             ! insert between two nodes
             pnew%next=>p2
             p1%next=>pnew
             matrix%lines(nl)%numnodes=matrix%lines(nl)%numnodes+1
             Nullify(pnew)
             Exit
          Endif
          p1=>p2
       Enddo
    Endif
    SM_set=0
  End Function SM_set

  
  Function SM_get(matrix, nl, nc, valid)
    Implicit None
    Real(dp) :: SM_get
    Type(SparseMatrixType), Intent(IN) :: matrix
    Integer, Intent(IN) :: nl, nc
    Logical, Intent(OUT) :: valid
    Type(SparseMatrixNode), Pointer :: p

    If(nl==nc) Then 
       valid=.True.
       SM_get=matrix%lines(nl)%first%data
    Else
       p=>matrix%lines(nl)%first%next
       Do While(Associated(p))
          If(p%colno==nc) Then
             valid=.True.
             SM_get=p%data
             Return
          Else If(p%colno>nc) Then
             Exit
          Endif
          p=>p%next
       Enddo
       valid=.False.
       SM_get=0.0
    Endif
  End Function SM_get

  
  Function SM_get_line_pointer(matrix, nl)
    Implicit None
    Type(SparseMatrixNode), Pointer :: SM_get_line_pointer
    Type(SparseMatrixType), Intent(IN) :: matrix
    Integer, Intent(IN) :: nl

    SM_get_line_pointer=>matrix%lines(nl)%first
  End Function SM_get_line_pointer


  Subroutine SM_write(mat)
    Use typedef
    Implicit None
    Type(SparseMatrixType), Intent(IN) :: mat
    Type(SparseMatrixNode), Pointer :: p
    Integer :: i

    Write(*,*) 'Printing sparse matrix with ',mat%nlines,' lines and ',mat%ncols,' columns'
    Do i=1,mat%nlines
       Write(*,*) 'line ',i,' has ',mat%lines(i)%numnodes,' elements'
       p=>mat%lines(i)%first
       Do While(Associated(p))
          Write(*,*) '   column=',p%colno,' data=',p%data
          p=>p%next
       Enddo
    Enddo
  End Subroutine SM_write


  Subroutine SM_multiply(mat, x, b)
    Use Typedef
    Implicit None
    Type(SparseMatrixType), Intent(IN) :: mat
    Real(dp), Intent(IN) :: x(:)
    Real(dp), Intent(OUT) :: b(:)
    Integer :: i
    Type(SparseMatrixNode), Pointer :: p

    Do i=1,mat%nlines
       p=>mat%lines(i)%first
       b(i)=p%data*x(i)
       p=>p%next
       Do While(Associated(p))
          b(i)=b(i)+p%data*x(p%colno)
          p=>p%next
       Enddo
    Enddo
  End Subroutine SM_multiply


  Subroutine SM_multiply_transpose(mat, x, b)
    Use Typedef
    Implicit None
    Type(SparseMatrixType), Intent(IN) :: mat
    Real(dp), Intent(IN) :: x(:)
    Real(dp), Intent(OUT) :: b(:)
    Integer :: i,j
    Type(SparseMatrixNode), Pointer :: p

    Do i=1,mat%nlines
       b(i)=0.0
    Enddo
    Do i=1,mat%nlines
       p=>mat%lines(i)%first
       Do While(Associated(p))
          j=p%colno
          b(j)=b(j)+p%data*x(i)
          p=>p%next
       Enddo
    Enddo
  End Subroutine SM_multiply_transpose

  
  Subroutine SM_divide_by_diagonal(mat, b, x)
    Implicit None
    Type(SparseMatrixType), Intent(IN) :: mat
    Real(dp), Intent(in) ::b(:)
    Real(dp), Intent(out) :: x(:)
    Integer :: i

    do i=1,mat%nlines
       x(i)=b(i)/mat%lines(i)%first%data
    Enddo
  End Subroutine SM_divide_by_diagonal

End Module SparseMatrix
