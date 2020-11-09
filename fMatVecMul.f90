!Caveats:
!	- type of tag_MPI and rank_MPI assumed to be same
!	- number of process must be >=2
!	- Master receives from any source that has completed a computation
!	- Multiple fast processes still have to wait in a queue to get their next job.
!Checked for:
!	- nprocs <= rows,cols <= nprocs
use mpi
Implicit none
!definitions for MPI
integer ierr_MPI, rank_MPI, nprocs_MPI,&
	tag_MPI,status_MPI(MPI_STATUS_SIZE),&
	req_MPI, count_MPI,& 
	Itags_MPI, Dtags_MPI

!definitions for fMatmul
integer I,J, cols, rows,&
	nRowSent
double precision, allocatable :: a_mat(:,:),b_vec(:),&
				 c_vec(:),&
				 a_row(:)
double precision tmp_buffer


cols = 2
rows = 2
allocate(a_mat(rows,cols),b_vec(rows))
CALL RANDOM_NUMBER(a_mat)
call random_number(b_vec)
!a_mat(1,1) = 1.1D0
!a_mat(1,2) = 1.2D0
!a_mat(1,3) = 1.3D0
!a_mat(2,1) = 1.4D0
!a_mat(2,2) = 1.5D0
!a_mat(2,3) = 1.6D0
!a_mat(3,1) = 1.7D0
!a_mat(3,2) = 1.8D0
!a_mat(3,3) = 1.9D0

!b_vec(1) = 1.1D0
!b_vec(2) = 2.1D0
!b_vec(3) = 3.1D0


!Initialize MPI
call MPI_INIT(ierr_MPI)

!Setup communicator size
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs_MPI,ierr_MPI)

!Setup rank
call MPI_COMM_RANK(MPI_COMM_WORLD, rank_MPI, ierr_MPI)

!Setup MPI based parameters
Dtags_MPI = 3 !The number of different tags sent by each slave process

!------------------------------------Main code block starts-----------------------------
   !allocate variables
        
   !check for process type (Master/Slave)
   if (rank_MPI .eq. 0) THEN
     !Master process---starts     
       write(*,1002)rank_MPI
       1002 format("Master's Process id =",i3)
       !allocatate variables for master
       allocate(c_vec(rows))
       !write the matrix A
         write(*,*)"The matrix A = "
         Do I = 1,rows
           Do J = 1, cols
             write(*,1001,ADVANCE='NO')a_mat(I,J)
             1001 format(f5.2," , ")
           enddo
           write(*,*)
         enddo !-loop for printing A
       !broadcast vector b
       call MPI_BCAST(b_vec, rows, MPI_DOUBLE_PRECISION,&
                   0, MPI_COMM_WORLD, ierr_MPI)
       write(*,*)"<Master> The b vector broadcast completed"                   
       !Initialize tasks to min(nprocs-1,rows) slaves
       nRowSent = 0 !Counter for sent rows
       Do Itags_MPI = 1,min(nprocs_MPI-1,rows)
         !Send the rows
         tag_MPI = Itags_MPI !-tag for which row is to be computed
         call MPI_SEND(a_mat(Itags_MPI,:),cols,MPI_DOUBLE_PRECISION,Itags_MPI,& 
                    tag_MPI, MPI_COMM_WORLD, ierr_MPI)
         nRowSent = nRowSent + 1
         write(*,*)"<Master> Sending initial rows"  
       enddo !-loop for task initialization
       write(*,*)"<Master> All initial rows sent"  
       !Terminate extra slaves (i.e. slave_pid > rows)
       Do Itags_MPI = rows+1,nprocs_MPI-1
         tag_MPI = 0 !-tag for terminating a slave
         call MPI_SEND(MPI_BOTTOM,0,MPI_DOUBLE_PRECISION,Itags_MPI,& 
                    tag_MPI, MPI_COMM_WORLD, ierr_MPI)
       enddo !-loop for terminating extra slaves
       write(*,*)"<Master> The extra slaves (if any) terminated"  
       !-------------------------------------------------Receive from slaves starts
       write(*,*)"<Master> Starting to receive from slaves"  
       Do Itags_MPI = 1,rows
         write(*,1004)Itags_MPI
         1004 format("<Master> Starting receipt of ",i3,"-th row")
         !Receive computed/updated sum_mat from slaves         
         call MPI_RECV(tmp_buffer,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,& 
                      MPI_ANY_TAG,MPI_COMM_WORLD,status_MPI,ierr_MPI) 
         write(*,1005)Itags_MPI
         1005 format("<Master> Received ",i3,"-th row successfully")     
         !Update c_vec
         c_vec(status_MPI(MPI_TAG)) = tmp_buffer         
         
         !Send the slave another unsent row if any, else terminate the slave
         if(nRowSent .lt. rows)THEN
           call MPI_SEND(a_mat(nRowSent+1,:),cols,MPI_DOUBLE_PRECISION,&
                    status_MPI(MPI_SOURCE),& 
                    nRowSent+1, MPI_COMM_WORLD, ierr_MPI)
           nRowSent = nRowSent + 1
         else
           tag_MPI = 0 !-tag for terminating a slave
           call MPI_SEND(MPI_BOTTOM,0,MPI_DOUBLE_PRECISION,&
                    status_MPI(MPI_SOURCE),& 
                    tag_MPI, MPI_COMM_WORLD, ierr_MPI)
         endif  !-if for another unsent row checkings         
       endDo !-loop for collecting all 1<=tags<=nprocs-1
       !-------------------------------------------------Receive from slaves stops
       !Master process---stops
       
     
   else ! 1<=ranks<=nprocs-1
     !-------------------------------------------------Send from slaves starts
     !Slave process
       write(*,1003)rank_MPI
       1003 format("Process id =",i3)
       !allocate slave specific variables
       allocate(a_row(cols))
       nRowSent = 0	!# of rows computed and sent
       
       !get the b_vec
       call MPI_BCAST(b_vec, rows, MPI_DOUBLE_PRECISION,&
                   0, MPI_COMM_WORLD, ierr_MPI)
       
       write(*,1006)rank_MPI
         1006 format("<Slave",i3,"> Received broadcast of b vector")
       write(*,1007)rank_MPI
         1007 format("<Slave",i3,"> Starting to receive rows of A")
  
       !compute elements of c_vec                  
         !get the row(s) of a_mat
           call MPI_RECV(a_row,cols,MPI_DOUBLE_PRECISION,0,& 
                      MPI_ANY_TAG,MPI_COMM_WORLD,status_MPI,ierr_MPI)           
           !store the tag
           tag_MPI = status_MPI(MPI_TAG)
           write(*,1008)rank_MPI,tag_MPI
           1008 format("<Slave",i3,"> Received ",i3,"-th row of A. Now checking for tag")
         !check for termination
         Do while(tag_MPI .ne. 0)     
           write(*,1009)rank_MPI,tag_MPI
           1009 format("<Slave",i3,"> Non-zero tag. Computing value for ",i3,"-th row of A")     
           !compute the dot product
           Do I = 1,cols
             a_row(I) = a_row(I)*b_vec(I)
           enddo !-loop for dot product
           write(*,1010)rank_MPI,tag_MPI
           1010 format("<Slave",i3,"> Sending ",i3,"-th row of c_vec.") 
           !send computed row element
           call MPI_SEND(sum(a_row),1,MPI_DOUBLE_PRECISION,0,& 
                    tag_MPI, MPI_COMM_WORLD, ierr_MPI) 
           nRowSent = nRowSent + 1
           write(*,1011)rank_MPI,tag_MPI
           1011 format("<Slave",i3,"> Sent ",i3,"-th row of c_vec. Waiting on next~~~")
           !Receive next task
           call MPI_RECV(a_row,cols,MPI_DOUBLE_PRECISION,0,& 
                      MPI_ANY_TAG,MPI_COMM_WORLD,status_MPI,ierr_MPI)           
           !store the tag
           tag_MPI = status_MPI(MPI_TAG)
         enddo             
       !-------------------------------------------------Send from slaves stops 
       !deallocate slave specific allocations
       deallocate(a_row)
       write(*,1012)rank_MPI,nRowSent 
       1012 format("<Slave",i3,"> computed a total of ",i3," rows.Terminating....")      
   endif !if-check for process type
!-------------------------------------Main code block stops-----------------------------
!Finialize MPI
call MPI_FINALIZE(ierr_MPI)	!No MPI calls after this line

!write(*,*)"Rank = ",rank_MPI
if(rank_MPI .eq. 0)THEN
  write(*,*)"The product vector = "
  Do I = 1,rows
    write(*,*)c_vec(I)
  enddo !-loop for printing A
endif

end
