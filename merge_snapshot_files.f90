!=====================================================
SUBROUTINE MERGE_SNAPSHOT_FILES(NPROC,NSNAP,BASE_NAME)
!=====================================================
!
IMPLICIT NONE
!
INTEGER :: IPROC,NPROC,NSNAP,N
INTEGER :: I,J,K,L
INTEGER :: I1,I2,J1,J2,K1,K2,X1,X2,Y1,Y2,Z1,Z2
REAL, DIMENSION(:,:,:), ALLOCATABLE :: X,Y,Z,UX,UY,UZ
CHARACTER :: FILENAME*100
CHARACTER(LEN=*) :: BASE_NAME
LOGICAL EXISTS
!
!==========================
!  MERGE SNAPSHOTS FILES  !
!==========================
!
WRITE(FILENAME,"(a,I3.3,a,I3.3)")TRIM(BASE_NAME),NSNAP,"_",0
INQUIRE(FILE=TRIM(FILENAME),EXIST=EXISTS)
IF(EXISTS)THEN
   !
   ! READ ALL FILE HEADER AND DETERMINE GRID DIMENSIONS
   !
   DO IPROC = 0,NPROC-1
      WRITE(FILENAME,"(A,I3.3,A,I3.3)")TRIM(BASE_NAME),NSNAP,"_",IPROC
      OPEN(26,FILE=TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM')
      READ(26)N,I1,I2,J1,J2,K1,K2
      IF(IPROC==0)THEN
         X1=I1 ; X2=I2 ; Y1=J1 ; Y2=J2 ; Z1=K1 ; Z2=K2 
      ELSE
         X1 = MIN(X1,I1)
         X2 = MAX(X2,I2)
         Y1 = MIN(Y1,J1)
         Y2 = MAX(Y2,J2)
         Z1 = MIN(Z1,K1)
         Z2 = MAX(Z2,K2)
      ENDIF
      CLOSE(26)
   ENDDO
   !
   ALLOCATE(X(x1:x2,y1:y2,z1:z2))
   ALLOCATE(Y(x1:x2,y1:y2,z1:z2))
   ALLOCATE(Z(x1:x2,y1:y2,z1:z2))
   !
   ALLOCATE(UX(x1:x2,y1:y2,z1:z2))
   ALLOCATE(UY(x1:x2,y1:y2,z1:z2))
   ALLOCATE(UZ(x1:x2,y1:y2,z1:z2))
   !
   ! MERGE FILES, READ/WRITE DATA
   !
   DO IPROC = 0,NPROC-1
      WRITE(FILENAME,"(A,I3.3,A,I3.3)")TRIM(BASE_NAME),NSNAP,"_",IPROC
      OPEN(26,FILE=TRIM(FILENAME),STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM')
      READ(26)N,I1,I2,J1,J2,K1,K2
      !      
      DO L = 1,N
         !
         READ(26)I        , J        , K
         READ(26)X(I,J,K) , Y(I,J,K) , Z(I,J,K)
         READ(26)UX(I,J,K), UY(I,J,K), UZ(I,J,K)
         !         
      ENDDO
      !      
      CLOSE(26, STATUS='DELETE')
      !
   ENDDO
   !
   ! METERS TO KILOMETERS
   ! 
   x = X/6371000.d0
   y = Y/6371000.d0
   z = Z/6371000.d0
   !
   WRITE(FILENAME,'(A,I3.3,A)')TRIM(BASE_NAME),NSNAP,'.vts' 
   CALL WRITE_VTS_SNAPSHOT(X1,X2,Y1,Y2,Z1,Z2,X,Y,Z,UX,UY,UZ,TRIM(FILENAME),27)
   !
   DEALLOCATE(X,Y,Z,UX,UY,UZ)
   !
ENDIF
!
RETURN
!
!==================================
END SUBROUTINE MERGE_SNAPSHOT_FILES
!==================================

!=================================================================================
SUBROUTINE WRITE_VTS_SNAPSHOT(x1,x2,y1,y2,z1,z2,x,y,z,ux,uy,uz,file_name,unit_vtk)
!=================================================================================
!
integer :: E_IO, Unit_VTK
integer :: i, x1, x2, y1, y2, z1, z2, offset
integer :: dim_arrays
character(1), parameter:: end_rec = char(10)
character*20 :: cx1, cx2, cy1, cy2, cz1, cz2, coffset
character(len=*) :: file_name
real, dimension(x1:x2,y1:y2,z1:z2) :: x, y, z, ux, uy, uz
!
write(cx1,'(I10)')x1
write(cx2,'(I10)')x2
write(cy1,'(I10)')y1
write(cy2,'(I10)')y2
write(cz1,'(I10)')z1
write(cz2,'(I10)')z2
!
open(unit       = Unit_VTK,        &
     file       = trim(file_name), &
     form       = 'UNFORMATTED',   &
     access     = 'STREAM',        &
     action     = 'WRITE',         &
     status     = 'REPLACE',       &
     convert    = 'LITTLE_ENDIAN', &
     iostat     = E_IO)
!
offset = 0
!
! WRITE HEADER
!
write(unit=Unit_VTK,iostat=E_IO)'<?xml version = "1.0"?>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<VTKFile type = "StructuredGrid" version="0.1" byte_order="LittleEndian">'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)'<StructuredGrid WholeExtent="'  &
                                        //trim(adjustl(cx1))//' '&
                                        //trim(adjustl(cx2))//' '&
                                        //trim(adjustl(cy1))//' '&
                                        //trim(adjustl(cy2))//' '&
                                        //trim(adjustl(cz1))//' '&
                                        //trim(adjustl(cz2))//'">'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)'<Piece Extent="'           &
                                   //trim(adjustl(cx1))//' '&
                                   //trim(adjustl(cx2))//' '&
                                   //trim(adjustl(cy1))//' '&
                                   //trim(adjustl(cy2))//' '&
                                   //trim(adjustl(cz1))//' '&
                                   //trim(adjustl(cz2))//'">'//end_rec
!
write(unit=Unit_VTK,iostat=E_IO)'<PointData>'//end_rec
write(coffset,'(I10)')offset
write(unit=Unit_VTK,iostat=E_IO)'<DataArray type="Float32" Name="displacement vector" NumberOfComponents="3"'&
&,' format="appended" offset="'//trim(adjustl(coffset))//'" />'//end_rec
offset = offset+sizeof(ux)+sizeof(uy)+sizeof(uz)+sizeof(i)
write(unit=Unit_VTK,iostat=E_IO)'</PointData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<CellData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</CellData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<Points>'//end_rec
write(coffset,'(I10)')offset
write(unit=Unit_VTK,iostat=E_IO)'<DataArray type="Float32" NumberOfComponents="3"'&
&,' format="appended" offset="'//trim(adjustl(coffset))//'" />'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</Points>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</Piece>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</StructuredGrid>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'<AppendedData encoding="raw">'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'_'
!
! FREE DISK SPACE FOR DATA 
!
dim_arrays = sizeof(ux)+sizeof(uy)+sizeof(uz)
write(unit=Unit_VTK,iostat=E_IO)dim_arrays
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=Unit_VTK,iostat=E_IO)ux(i,j,k),uy(i,j,k),uz(i,j,k)
      enddo
   enddo
enddo
dim_arrays = sizeof(x)+sizeof(y)+sizeof(z)
write(unit=Unit_VTK,iostat=E_IO)dim_arrays
do k = z1,z2
   do j = y1,y2
      do i = x1,x2
         write(unit=Unit_VTK,iostat=E_IO)x(i,j,k),y(i,j,k),z(i,j,k)
      enddo
   enddo
enddo
write(unit=Unit_VTK,iostat=E_IO)end_rec
!
! WRITE FOOTER
!
write(unit=Unit_VTK,iostat=E_IO)'</AppendedData>'//end_rec
write(unit=Unit_VTK,iostat=E_IO)'</VTKFile>'//end_rec
!
close(unit=Unit_VTK)
!
RETURN
!
!================================
END SUBROUTINE WRITE_VTS_SNAPSHOT
!================================

