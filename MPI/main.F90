!fortran 非阻塞通信例子
!0号进程数组test_array(1,2,3,4,5),1号进程数组test_array(6,7,8,9,10),使用mpi非阻塞通信交换不同进程数组的后两个值
!编译 mpif90 main.F90 -o main
!运行 mpirun -n 2 ./main
        module grist_mpi
          include 'mpif.h' 
          contains
              
            function mpi_rank(comm) result(res)
            implicit none
                integer, optional :: comm
                    integer :: ierr
                        integer :: res 
                            
                            if(present(comm)) then
                                   call MPI_COMM_RANK(Comm, res, ierr)
            else
                   call MPI_COMM_RANK(MPI_COMM_WORLD, res, ierr)
            end if
                
              end function

                function mpi_size(comm) result(res)
            implicit none
                integer, optional :: comm
                    integer :: ierr
                        integer :: res 
                            
                            if(present(comm)) then
                                   call MPI_COMM_SIZE(Comm, res, ierr)
            else
                   call MPI_COMM_SIZE(MPI_COMM_WORLD, res, ierr)
            end if
                
              end function

                
              end module grist_mpi

        
        module test
                use grist_mpi
        contains
                subroutine test_mpi()
                        implicit none
                        integer :: MPI_COMMM_WORLD,ierr
                        real(8),target  :: test_array(5)
                        real(8),pointer :: q(:)
                        real(8) :: data_send(2),data_recv(2)
                       integer,allocatable :: reqs_send(:), reqs_recv(:)
                       integer,allocatable :: status(:,:)
                        call mpi_init(ierr)

                        
                        if(mpi_rank()==0)then
                                test_array=(/1.0,2.0,3.0,4.0,5.0/)
                        else
                                test_array=(/6.0,7.0,8.0,9.0,10.0/)
                        endif
                        
                        allocate(reqs_send(1))
                        allocate(reqs_recv(1))
                        q=>test_array(4:5)
                        data_send(1)=test_array(4)
                        data_send(2)=test_array(5)

                        if(mpi_rank()==0) then
               call MPI_Isend(data_send(1),size(data_send),MPI_REAL8,1,&
                             102,MPI_COMM_WORLD,reqs_send(1),ierr)
                        else
               call MPI_Isend(data_send(1),size(data_send),MPI_REAL8,0,&
                             102,MPI_COMM_WORLD,reqs_send(1),ierr)
                        endif

                        if(mpi_rank()==0) then
               call MPI_Irecv(data_recv(1),size(data_recv),MPI_REAL8,1,&
                             102,MPI_COMM_WORLD,reqs_recv(1),ierr)
                        else
               call MPI_Irecv(data_recv(1),size(data_recv),MPI_REAL8,0,&
                             102,MPI_COMM_WORLD,reqs_recv(1),ierr)
                        endif

                        allocate(status(MPI_STATUS_SIZE,1))
                        call MPI_WaitAll(1,reqs_send,status,ierr)
                        call MPI_WaitAll(1,reqs_recv,status,ierr)

                        q(1)=data_recv(1)
                        q(2)=data_recv(2)
                        print*,mpi_rank(),test_array
                        call mpi_finalize(ierr)
                
                end subroutine
        end module test


        program main
                use test
                implicit none
                call test_mpi()
        end program main
