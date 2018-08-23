
!------------------------------------------------------------------------------------------------------------------------------------------------
! ichannel = 1
!
!  
!  rotation 1: Phi(k, k+q; k'-k) -> Phi(k1, k1'; q1) 
!
do i = 1, Nt      ! k
    i1 = i         ! k1 = k
    do k1 = 1, Nb  ! q1
        call index_operation_FaddB(Index_Fermionic(i), index_bosonic(inode*Nb+k1), ComIdx1)  ! k'=k+q1
        if (ComIdx1%iw >=1 .and. ComIdx1%iw <= Nf) then
            j = list_index_Fermionic(ComIdx1)
            do k = 1, Nb
                call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx2)
                j1 = list_index_Fermionic(ComIdx2)   ! k1' = k + q

                if (ComIdx2%iw <= Nf .and. ComIdx2%iw >= 1) then
                    F_d(i, j, k) = F_d(i, j, k) - 0.5d0*mat(i1, j1, k1)
                    F_m(i, j, k) = F_m(i, j, k) - 0.5d0*mat(i1, j1, k1)
                else
                    F_d(i, j, k) = F_d(i, j, k) - 0.5d0*Kernel(CHANNEL_D, index_fermionic(i1), ComIdx2, index_bosonic(inode*Nb+k1))
                    F_m(i, j, k) = F_m(i, j, k) - 0.5d0*Kernel(CHANNEL_D, index_fermionic(i1), ComIdx2, index_bosonic(inode*Nb+k1))
                end if
            end do
        end if
    end do

    ! Time reversal symmetric part: Phi(k, k+q; k'-k) -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1) 
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k1 = 1, Nb
        if (Index_Bosonic(inode*Nb+k1)%iw > 1) then
            call index_operation_FaddB(ComIdx1, index_bosonic(inode*Nb+k1), ComIdx2)  ! -k+q1
            call index_operation_MinusF(ComIdx2, ComIdx2, ComIdx3) ! k'=k-q1
            if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                do k = 1, Nb
                    call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx4) ! k+q
                    call index_operation_MinusF(ComIdx4, ComIdx4, ComIdx5)
                    j1 = list_index_Fermionic(ComIdx5)  
                    if (ComIdx5%iw <= Nf .and. ComIdx5%iw >= 1) then
                        F_d(i, j, k) = F_d(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                        F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                    else
                        F_d(i, j, k) = F_d(i, j, k) - 0.5d0*conjg(Kernel(CHANNEL_D, index_fermionic(i1), ComIdx5, index_bosonic(inode*Nb+k1)))
                        F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(Kernel(CHANNEL_D, index_fermionic(i1), ComIdx5, index_bosonic(inode*Nb+k1)))
                    end if
                end do
            end if
        end if
    end do
end do

!
! rotation 2: Phi(k, q-k'; k'-k) -> Phi(k1, k1'; q1)
!
do i = 1, Nt   ! k
    i1 = i      ! k1 = k
    do k1 = 1, Nb   ! q1
        call index_operation_FaddB(index_fermionic(i), index_bosonic(inode*Nb+k1), ComIdx1)  ! k'=k+q1
        if (ComIdx1%iw >=1 .and. ComIdx1%iw <= Nf) then
            j = list_index_Fermionic(ComIdx1)
            do k = 1, Nb  ! q
                call index_operation_MinusF(ComIdx1, ComIdx1,ComIdx2) ! -k'
                call index_operation_FaddB(ComIdx2, index_bosonic(id*Nb+k), ComIdx3) ! q-k'
                if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then
                    j1 = list_index_Fermionic(ComIdx3)
                    F_s(i, j, k) = F_s(i, j, k) + 0.5d0*mat(i1, j1, k1)
                    F_t(i, j, k) = F_t(i, j, k) - 0.5d0*mat(i1, j1, k1)
                else
                    F_s(i, j, k) = F_s(i, j, k) + 0.5d0*Kernel(CHANNEL_D, index_fermionic(i1), ComIdx3, index_bosonic(inode*Nb+k1))
                    F_t(i, j, k) = F_t(i, j, k) - 0.5d0*Kernel(CHANNEL_D, index_fermionic(i1), ComIdx3, index_bosonic(inode*Nb+k1))
                end if
            end do
        end if
    end do

    ! Time reversal symmetric part: Phi(k, q-k'; k'-k) -> Phi(-k1, -k1'; -q1)=conjg(Phi(k1, k1'; q1))
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k1 = 1, Nb ! q1
        if (index_Bosonic(inode*Nb+k1)%iw > 1) then
            call index_operation_FaddB(ComIdx1, index_bosonic(inode*Nb+k1), ComIdx2) ! -k+q1
            call index_operation_MinusF(ComIdx2, ComIdx2, ComIdx3) ! k'=k-q1
            if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                do k = 1, Nb
                    call index_operation_MinusF(ComIdx3, ComIdx3, ComIdx4) ! -k'
                    call index_operation_FaddB(ComIdx4, index_bosonic(id*Nb+k), ComIdx5) ! q-k'
                    call index_operation_MinusF(ComIdx5, ComIdx5, ComIdx4) ! k1' = k'-q 
                    if (ComIdx4%iw >=1 .and. ComIdx4%iw <= Nf) then
                        j1 = list_index_Fermionic(ComIdx4)
                        F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                        F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                    else
                        F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(Kernel(CHANNEL_D, index_fermionic(i1), ComIdx4, index_bosonic(inode*Nb+k1)))
                        F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(Kernel(CHANNEL_D, index_fermionic(i1), ComIdx4, index_bosonic(inode*Nb+k1)))
                    end if
                end do
            end if
        end if
    end do
end do

!
! rotation 3: Phi(k, k'; q-k-k') -> Phi(k1, k1'; q1)
!
do i = 1, Nt   ! k
    i1 = i      ! k1 = k
    do k1 = 1, Nb   ! q1
        call index_operation_FaddB(index_fermionic(i), index_bosonic(inode*Nb+k1), ComIdx1)  ! k+q1
        call index_operation_MinusF(ComIdx1, ComIdx1, ComIdx2) ! -k-q1
        do k = 1, Nb  ! q            
            call index_operation_FaddB(ComIdx2, index_bosonic(id*Nb+k), ComIdx3) ! k'=q-k-q1            
            if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                    !j1 = j!list_index_Fermionic(ComIdx3)
                F_s(i, j, k) = F_s(i, j, k) + 0.5d0*mat(i, j, k1)
                F_t(i, j, k) = F_t(i, j, k) + 0.5d0*mat(i, j, k1)
            end if
        end do
    end do

    ! Time reversal symmetric part: Phi(k, k'; q-k-k') -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1)            
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k1 = 1, Nb ! q1
        if (index_Bosonic(inode*Nb+k1)%iw > 1) then
            call index_operation_FaddB(ComIdx1, index_bosonic(inode*Nb+k1), ComIdx2)  ! -k+q1          
            do k = 1, Nb  ! q            
                call index_operation_FaddB(ComIdx2, index_bosonic(id*Nb+k), ComIdx3) ! k'=q-k+q1
                if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                    j = list_index_Fermionic(ComIdx3)                   
                    call index_operation_MinusF(ComIdx3, ComIdx3, ComIdx4) ! -k'
                    if (ComIdx4%iw >=1 .and. ComIdx4%iw <= Nf) then
                        j1 = list_index_Fermionic(ComIdx4)
                        F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                        F_t(i, j, k) = F_t(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                    end if
                end if
            end do !Nb 
        end if !q1%iw>0
    end do !Nb
end do !Nt   

!------------------------------------------------------------------------------------------------------------------------------------------------
! ichannel = 2
!
!  
!  rotation 1: Phi(k, k+q; k'-k) -> Phi(k1, k1'; q1) 
!
do i = 1, Nt      ! k
    i1 = i         ! k1 = k
    do k1 = 1, Nb  ! q1
        call index_operation_FaddB(Index_Fermionic(i), index_bosonic(inode*Nb+k1), ComIdx1)  ! k'=k+q1
        if (ComIdx1%iw >=1 .and. ComIdx1%iw <= Nf) then
            j = list_index_Fermionic(ComIdx1)
            do k = 1, Nb
                call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx2)
                j1 = list_index_Fermionic(ComIdx2)   ! k1' = k + q

                if (ComIdx2%iw <= Nf .and. ComIdx2%iw >= 1) then
                    F_d(i, j, k) = F_d(i, j, k) - 1.5d0*mat(i1, j1, k1)
                    F_m(i, j, k) = F_m(i, j, k) + 0.5d0*mat(i1, j1, k1)
                else
                    F_d(i, j, k) = F_d(i, j, k) - 1.5d0*Kernel(CHANNEL_M, index_fermionic(i1), ComIdx2, index_bosonic(inode*Nb+k1))
                    F_m(i, j, k) = F_m(i, j, k) + 0.5d0*Kernel(CHANNEL_M, index_fermionic(i1), ComIdx2, index_bosonic(inode*Nb+k1))
                end if
            end do
        end if
    end do

    ! Time reversal symmetric part: Phi(k, k+q; k'-k) -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1) 
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k1 = 1, Nb
        if (Index_Bosonic(inode*Nb+k1)%iw > 1) then
            call index_operation_FaddB(ComIdx1, index_bosonic(inode*Nb+k1), ComIdx2)  ! -k+q1
            call index_operation_MinusF(ComIdx2, ComIdx2, ComIdx3) ! k'=k-q1
            if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                do k = 1, Nb
                    call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx4) ! k+q
                    call index_operation_MinusF(ComIdx4, ComIdx4, ComIdx5)
                    j1 = list_index_Fermionic(ComIdx5)  
                    if (ComIdx5%iw <= Nf .and. ComIdx5%iw >= 1) then
                        F_d(i, j, k) = F_d(i, j, k) - 1.5d0*conjg(mat(i1, j1, k1))
                        F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                    else
                        F_d(i, j, k) = F_d(i, j, k) - 1.5d0*conjg(Kernel(CHANNEL_M, index_fermionic(i1), ComIdx5, index_bosonic(inode*Nb+k1)))
                        F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(Kernel(CHANNEL_M, index_fermionic(i1), ComIdx5, index_bosonic(inode*Nb+k1)))
                    end if
                end do
            end if
        end if
    end do
end do

!
! rotation 2: Phi(k, q-k'; k'-k) -> Phi(k1, k1'; q1)
!
do i = 1, Nt   ! k
    i1 = i      ! k1 = k
    do k1 = 1, Nb   ! q1
        call index_operation_FaddB(index_fermionic(i), index_bosonic(inode*Nb+k1), ComIdx1)  ! k'=k+q1
        if (ComIdx1%iw >=1 .and. ComIdx1%iw <= Nf) then
            j = list_index_Fermionic(ComIdx1)
            do k = 1, Nb  ! q
                call index_operation_MinusF(ComIdx1, ComIdx1,ComIdx2) ! -k'
                call index_operation_FaddB(ComIdx2, index_bosonic(id*Nb+k), ComIdx3) ! q-k'
                if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then
                    j1 = list_index_Fermionic(ComIdx3)
                    F_s(i, j, k) = F_s(i, j, k) - 1.5d0*mat(i1, j1, k1)
                    F_t(i, j, k) = F_t(i, j, k) - 0.5d0*mat(i1, j1, k1)
                else
                    F_s(i, j, k) = F_s(i, j, k) - 1.5d0*Kernel(CHANNEL_M, index_fermionic(i1), ComIdx3, index_bosonic(inode*Nb+k1))
                    F_t(i, j, k) = F_t(i, j, k) - 0.5d0*Kernel(CHANNEL_M, index_fermionic(i1), ComIdx3, index_bosonic(inode*Nb+k1))
                end if
            end do
        end if
    end do

    ! Time reversal symmetric part: Phi(k, q-k'; k'-k) -> Phi(-k1, -k1'; -q1)=conjg(Phi(k1, k1'; q1))
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k1 = 1, Nb ! q1
        if (index_Bosonic(inode*Nb+k1)%iw > 1) then
            call index_operation_FaddB(ComIdx1, index_bosonic(inode*Nb+k1), ComIdx2) ! -k+q1
            call index_operation_MinusF(ComIdx2, ComIdx2, ComIdx3) ! k'=k-q1
            if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                do k = 1, Nb
                    call index_operation_MinusF(ComIdx3, ComIdx3, ComIdx4) ! -k'
                    call index_operation_FaddB(ComIdx4, index_bosonic(id*Nb+k), ComIdx5) ! q-k'
                    call index_operation_MinusF(ComIdx5, ComIdx5, ComIdx4) ! k1' = k'-q 
                    if (ComIdx4%iw >=1 .and. ComIdx4%iw <= Nf) then
                        j1 = list_index_Fermionic(ComIdx4)
                        F_s(i, j, k) = F_s(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                        F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                    else
                        F_s(i, j, k) = F_s(i, j, k) - 1.5d0*conjg(Kernel(CHANNEL_M, index_fermionic(i1), ComIdx4, index_bosonic(inode*Nb+k1)))
                        F_t(i, j, k) = F_t(i, j, k) - 0.5d0*conjg(Kernel(CHANNEL_D, index_fermionic(i1), ComIdx4, index_bosonic(inode*Nb+k1)))
                    end if
                end do
            end if
        end if
    end do
end do

!
! rotation 3: Phi(k, k'; q-k-k') -> Phi(k1, k1'; q1)
!
do i = 1, Nt   ! k
    i1 = i      ! k1 = k
    do k1 = 1, Nb   ! q1
        call index_operation_FaddB(index_fermionic(i), index_bosonic(inode*Nb+k1), ComIdx1)  ! k+q1
        call index_operation_MinusF(ComIdx1, ComIdx1, ComIdx2) ! -k-q1
        do k = 1, Nb  ! q            
            call index_operation_FaddB(ComIdx2, index_bosonic(id*Nb+k), ComIdx3) ! k'=q-k-q1            
            if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                !j1 = j!list_index_Fermionic(ComIdx3)
                F_s(i, j, k) = F_s(i, j, k) - 1.5d0*mat(i, j, k1)
                F_t(i, j, k) = F_t(i, j, k) + 0.5d0*mat(i, j, k1)
            end if
        end do
    end do

    ! Time reversal symmetric part: Phi(k, k'; q-k-k') -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1)            
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k1 = 1, Nb ! q1
        if (index_Bosonic(inode*Nb+k1)%iw > 1) then
            call index_operation_FaddB(ComIdx1, index_bosonic(inode*Nb+k1), ComIdx2)  ! -k+q1          
            do k = 1, Nb  ! q            
                call index_operation_FaddB(ComIdx2, index_bosonic(id*Nb+k), ComIdx3) ! k'=q-k+q1
                if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                    j = list_index_Fermionic(ComIdx3)                   
                    call index_operation_MinusF(ComIdx3, ComIdx3, ComIdx4) ! -k'
                    if (ComIdx4%iw >=1 .and. ComIdx4%iw <= Nf) then
                        j1 = list_index_Fermionic(ComIdx4)
                        F_s(i, j, k) = F_s(i, j, k) - 1.5d0*conjg(mat(i1, j1, k1))
                        F_t(i, j, k) = F_t(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                    end if
                end if
            end do !Nb 
        end if !q1%iw>0
    end do !Nb
end do !Nt

!------------------------------------------------------------------------------------------------------------------------------
! ichannel = 3
!
do i = 1, Nt   ! k
    do k = 1, Nb   ! q
        call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx1)  ! k+q
        call index_operation_MinusF(ComIdx1, ComIdx1, ComIdx2) ! -k-q
        do k1 = 1, Nb  ! q1            
            call index_operation_FaddB(ComIdx2, index_bosonic(inode*Nb+k1), ComIdx3) ! k'=q1-k-q            
            if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                F_d(i, j, k) = F_d(i, j, k) + 0.5d0*mat(i, j, k1) 
                F_m(i, j, k) = F_m(i, j, k) - 0.5d0*mat(i, j, k1)
            end if
        end do
    end do

! Time reversal symmetric part: Phi(k, k'; q-k-k') -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1)            
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k = 1, Nb ! q
        call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx2)  ! k+q
        do k1 = 1, Nb  ! q1            
            if (index_Bosonic(inode*Nb+k1)%iw > 1) then  
                call index_operation_FaddB(ComIdx2, index_bosonic(inode*Nb+k1), ComIdx3) ! -k'=q1+k+q
                if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                    j1 = list_index_Fermionic(ComIdx3)              
                    call index_operation_MinusF(ComIdx3, ComIdx3, ComIdx4) ! k'=-q1-k-q
                    if (ComIdx4%iw >=1 .and. ComIdx4%iw <= Nf) then
                        j = list_index_Fermionic(ComIdx4)
                        F_d(i, j, k) = F_d(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                        F_m(i, j, k) = F_m(i, j, k) - 0.5d0*conjg(mat(i1, j1, k1))
                    end if
                end if
            end if
        end do
    end do
end do

!------------------------------------------------------------------------------------------------------------------------------
! ichannel = 4
!
do i = 1, Nt   ! k
    do k = 1, Nb   ! q
        call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx1)  ! k+q
        call index_operation_MinusF(ComIdx1, ComIdx1, ComIdx2) ! -k-q
        do k1 = 1, Nb  ! q1            
            call index_operation_FaddB(ComIdx2, index_bosonic(inode*Nb+k1), ComIdx3) ! k'=q1-k-q            
            if (ComIdx3%iw >=1 .and. ComIdx3%iw <= Nf) then
                j = list_index_Fermionic(ComIdx3)
                F_d(i, j, k) = F_d(i, j, k) + 1.5d0*mat(i, j, k1)
                F_m(i, j, k) = F_m(i, j, k) + 0.5d0*mat(i, j, k1)
            end if
        end do
    end do

! Time reversal symmetric part: Phi(k, k'; q-k-k') -> Phi(-k1, -k1'; -q1) = conjg(k1, k1'; q1)            
    call index_operation_MinusF(Index_Fermionic(i), index_Fermionic(i), ComIdx1) ! -k
    i1 = list_index_Fermionic(ComIdx1)  ! k1 = -k
    do k = 1, Nb ! q
        call index_operation_FaddB(index_fermionic(i), index_bosonic(id*Nb+k), ComIdx2)  ! k+q
        do k1 = 1, Nb  ! q1            
            if (index_Bosonic(inode*Nb+k1)%iw > 1) then  
                call index_operation_FaddB(ComIdx2, index_bosonic(inode*Nb+k1), ComIdx3) ! -k'=q1+k+q
                if (ComIdx3%iw >= 1 .and. ComIdx3%iw <= Nf) then
                    j1 = list_index_Fermionic(ComIdx3)              
                    call index_operation_MinusF(ComIdx3, ComIdx3, ComIdx4) ! k'=-q1-k-q
                    if (ComIdx4%iw >=1 .and. ComIdx4%iw <= Nf) then
                        j = list_index_Fermionic(ComIdx4)
                        F_d(i, j, k) = F_d(i, j, k) + 1.5d0*conjg(mat(i1, j1, k1))
                        F_m(i, j, k) = F_m(i, j, k) + 0.5d0*conjg(mat(i1, j1, k1))
                    end if
                end if
            end if
        end do
    end do
end do
