!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht, Stefan Grimme
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!=========================================================================================!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!=========================================================================================!
!===============================================================================!
!
!  This is the CCEGEN routine, used for clustering a conformational ensemble
!  and determine representative structures for the molecule.
!
!  A principal component analysis (PCA) is performed, and the generated data
!  is clustered.
!
!
!  On Input:   pr - printout boolean
!              fname - name of the ensemble file
!
!==============================================================================!
subroutine CCEGEN(env,pr,fname)
  use crest_parameters,idp => dp
  use crest_data
  use zdata
  use strucrd
  use utilities
  implicit none
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  type(timer) :: ctimer
  logical,intent(in)   :: pr
  character(len=*),intent(in) :: fname
  type(zmolecule) :: zmol
  type(zequal) :: groups
  type(zequal) :: subgroups
  integer,allocatable :: inc(:)
  logical :: heavyonly
  integer :: i,j,k,l,ich,c
  integer :: nat,nall
  real(wp) :: dum,dum2
  type(coord),allocatable :: mols(:)

  character(len=:),allocatable :: measuretype

  !>--- SVD params
  integer :: ntaken
  integer :: nallnew
  real(wp),allocatable :: xyznew(:,:,:)
  !integer,allocatable  :: atnew(:)
  real(wp),allocatable :: measure(:,:) !this is what is passed to the SVD
  integer :: mn,mm                     !--> measure(mn,mm)
  real(wp),allocatable :: pc(:)        !the principal components
  real(wp),allocatable :: pcvec(:,:)   !the principal component eigenvectors
  real(wp),allocatable :: pcdum(:,:)
  integer :: nbnd,ndied
  integer,allocatable :: diedat(:,:)   !atoms spanning relevant dihedral angles
  real(wp),allocatable :: diedr(:)
  real(wp) :: pcsum
  real(wp) :: pcthr
  real(wp) :: pcmin
  integer :: pccap
  integer :: npc
  real(wp),allocatable :: geo(:,:)
  integer,allocatable :: na(:),nb(:),nc(:)

  !>--- CLUSTERING params
  character(len=:),allocatable :: clusteralgo
  integer :: nclust                    !number of clusters
  integer :: nclustiter                !iterator for nclust
  integer :: nclustmin,nclustmax
  integer,allocatable :: member(:)     !track cluster correspondence
  real(ap),allocatable :: p(:),q(:)
  real(sp),allocatable :: dist(:)
  real(ap),allocatable :: centroid(:,:)
  integer(idp) :: ndist,klong
  real(ap) :: eucdist !this is a function
  real(wp) :: DBI,pSF,SSRSST,SSRSSTthr
  real(wp) :: csthr
  integer :: ncb,ancb
  real(wp),allocatable :: eclust(:)
  integer,allocatable :: clustbest(:),ind(:)
  real(wp),allocatable :: statistics(:,:)
  integer,allocatable :: clust_sizes(:)
  logical,allocatable :: extrema(:,:)
  logical :: autolimit
  real(wp) :: fraclimit

  !>--- printout and params
  real(wp) :: emin,erel

  call ctimer%init(20)
  if (pr) then
    call largehead('Principal Component Analysis (PCA) and Clustering')
    write (stdout,'(1x,a,a)') 'Input file: ',trim(fname)
  end if

!=========================================================!
! set threads
!=========================================================!
  call cregen_setthreads(stdout,env,pr)

!=========================================================!
! Prepare a coordinate ensemble for the clustering
!=========================================================!
  !>--- 0. Set defaults, read ensemble
  call rdensemble(fname,nall,mols)
  if (nall < 1) then
    error stop "Ensemble is empty! must stop"
  end if
  nat = mols(1)%nat
  if (nall == 1) then
    if (pr) then
      write (stdout,*) 'Only one structure in ensemble!'
      write (stdout,*) 'Write structure to ',clusterfile,' and skip PCA parts'
    end if
    open (newunit=ich,file=clusterfile)
    call mols(1)%append(ich)
    close (ich)
    deallocate (mols)
    return
  end if

  heavyonly = .true.
  !measuretype = 'dihedral'
  measuretype = env%pcmeasure
  clusteralgo = 'kmeans'
  !pcthr = 0.85d0  !PCs must add up to this amount of "representability"
  pcthr = env%pcthr
  pcmin = env%pcmin
  !csthr = 0.80d0  !threshold for SSR/SST to select a suitable cluster count
  csthr = env%csthr
  !pccap = 100   !a cap for the number of principal components used in the clustering
  pccap = env%pccap
  autolimit = .true. !if the ensemble is very large, take only a fraction to speed up things
  !(only for predefined clustersizes with "-cluster N" )
  fraclimit = 0.25d0 !if autolimit=true, 1/4 of the ensemble is taken

  !>--- override pcthr and csthr based on clustering level
  if (env%maxcluster == 0) then
    call clustleveval(env,nclustmax,csthr,SSRSSTthr,pcthr)
  end if

  !>--- 1. topology for reference structure
  if (env%wbotopo) then
    env%wbofile = 'wbo'
  else
    env%wbofile = 'none given'
  end if
  call simpletopo(nat,mols(1)%at,mols(1)%xyz,zmol,pr,.false.,env%wbofile)
  allocate (inc(zmol%nat),source=0)

!===========================================================!
  if (measuretype .ne. 'dihedral') then
!===========================================================!
    !>--- 2. read nuclear equivalencies
    if (pr) then
      write (stdout,*)
      call smallhead('READING NUCLEAR EQUIVALENCIES')
    end if
    call readequals('anmr_nucinfo',zmol,groups)
    if (pr) then
      call groups%prsum(6) !--- print summary to screen
      write (stdout,'(1x,a)') 'Unlisted nuclei (groups) are unique.'
    end if

    !>--- 3. distribute groups into subgroups basedon topology
    if (pr) then
      write (stdout,*)
      call smallhead('ANALYZING EQUIVALENCIES')
    end if
    call distsubgr(zmol,groups,subgroups,inc,pr)

    !>--- 4. Equivalent atoms must be excluded in clustering to reduce noise
    if (pr) then
      write (stdout,*)
      call smallhead('DETERMINE ATOMS TO INCLUDE IN PCA')
    end if
    call excludeFromRMSD(zmol,inc)
    if (sum(inc) == 0) then
      if (pr) then
        write (stdout,*) 'WARNING: No atoms included in PCA'
        write (stdout,*) 'Including more atoms ...'
      end if
      inc = 1
      !-- for this exclude the equivalent atoms from anmr_nucinfo directly
      do i = 1,groups%ng
        if (groups%grp(i)%nm > 1) then
          write (stdout,*) groups%grp(i)%mem
          do j = 1,groups%grp(i)%nm
            k = groups%grp(i)%mem(j)
            inc(k) = 0
          end do
        end if
      end do
    end if
    !>-- exclude user set atoms
    if (env%pcaexclude) then
      call excludeSelected(zmol,inc,env%atlist)
    end if
    !>-- exclude H atoms
    if (heavyonly) then
      call excludeLight(zmol,inc)
    end if
    if (pr) then
      do i = 1,zmol%nat
        if (inc(i) == 1) then
          write (stdout,'(1x,a,a,i0,a,5x,a)') zmol%zat(i)%el,'(',i,')','taken'
        end if
      end do
    end if
    ntaken = sum(inc)
    !>--- if we have too few, include the heav atoms at least
    if (ntaken <= 3) then
      do i = 1,zmol%nat
        if (zmol%at(i) /= 1) then
          inc(i) = 1
          if (pr) write (stdout,'(1x,a,a,i0,a,5x,a)') zmol%zat(i)%el,'(',i,')','taken'
        end if
      end do
    end if
    ntaken = sum(inc)

    call zmol%deallocate

    !>-- for very large ensemble files limit the clustering
    if (autolimit) then
      if ((env%nclust /= 0).and.(env%nclust*100 < nall)) then
        dum = float(nall)*fraclimit
        dum2 = float(env%nclust)
        nallnew = nint(max(dum,dum2))
      else
        nallnew = nall
      end if
    else
      nallnew = nall
    end if

    !>--- 5. Transfer the relevant atoms to a new array
    allocate (xyznew(3,ntaken,nallnew))
    do i = 1,nallnew
      k = 0
      do j = 1,nat
        if (inc(j) == 1) then
          k = k+1
          xyznew(:,k,i) = mols(i)%xyz(:,j)
        end if
      end do
    end do

!===================================================!
  else !measuretype=='dihedral'
!===================================================!

    !-- for very large ensemble files limit the clustering
    if (autolimit) then
      if ((env%nclust /= 0).and.(env%nclust*100 < nall)) then
        dum = float(nall)*fraclimit
        dum2 = float(env%nclust)
        nallnew = nint(max(dum,dum2))
      else
        nallnew = nall
      end if
    else
      nallnew = nall
    end if

    inc = 1
    ntaken = sum(inc)

    call zmol%countbonds()
    nbnd = zmol%nb
    allocate (diedat(4,zmol%nb),source=0)
    call getdiederatoms(zmol,zmol%nat,inc,nbnd,diedat,ndied)
    ntaken = ndied
!==================================================!
  end if
!==================================================!

!===================================================================================================!
! do the SVD to get the principal components
!===================================================================================================!
  if (ntaken > 3) then !> all of this only makes sense if we have something to compare
    call ctimer%start(1,'PCA')
    if (pr) then
      write (stdout,*)
      call smallhead('PRINCIPAL COMPONENT ANALYSIS')
    end if
    mm = nallnew
    select case (measuretype)
      !==========================================================================!
    case ('cma','CMA','cmadist')
      if (pr) then
        write (stdout,'(1x,a)') 'Using CMA DISTANCES as descriptors:'
      end if
      !>-- all structures should have been shifted to the CMA by CREGEN
      !>   therefore assume the CMA is at (0,0,0).
      !>   somewhat robust measure, but provides less information.
      mn = min(ntaken,mm)
      allocate (measure(mn,mm),pc(mn),pcvec(mm,mn))
      do i = 1,mm
        do j = 1,mn
          measure(j,i) = xyznew(1,j,i)**2+ &
       &                 xyznew(2,j,i)**2+ &
       &                 xyznew(3,j,i)**2
          measure(j,i) = sqrt(measure(j,i))
        end do
      end do
      !==========================================================================!
    case ('cartesian','coords')
      if (pr) then
        write (stdout,'(1x,a)') 'Using CARTESIAN COORDINATES as descriptors:'
      end if
      !>-- all Cartesian components of the selected atoms
      !>   REQUIRES PERFECT ALIGNMENT(!), hence not very robust
      mn = min(ntaken*3,mm)
      allocate (measure(mn,mm),pc(mn),pcvec(mm,mn))
      do i = 1,mm
        l = 0
        do j = 1,ntaken
          do k = 1,3
            l = l+1
            if (l > mn) exit
            measure(l,i) = xyznew(k,j,i)
          end do
          if (l > mn) exit
        end do
      end do
      !==========================================================================!
    case default !case( 'zmat','zmatrix' )
      if (pr) then
        write (stdout,'(1x,a)') 'Using ZMATRIX as descriptors (sin/cos of dihedrals):'
      end if
      !>-- dihedral angles, sin/cos transformed for periodicity
      l = ntaken-3 !>--- first three dihedral angles are zero
      mn = min(mm,2*l) !>--- two descriptors per dihedral, no more than structures
      if (mn < 2) then !> we need at least 1 dihedral angle (2 descriptors)
        if (pr) then
          write (stdout,*) "Not enough descriptors for PCA!"
          return
        end if
      end if
      allocate (measure(mn,mm),pc(mn),pcvec(mm,mn))
      allocate (geo(3,ntaken),source=0.0d0)
      allocate (na(ntaken),nb(ntaken),nc(ntaken))
      do i = 1,mm
        na = 0; nb = 0; nc = 0
        geo = 0.0d0
        call xyzint(xyznew(1:3,1:ntaken,i),ntaken,na,nb,nc,radtodeg,geo)
        do j = 1,mn/2
          k = j+3
          dum = geo(3,k)*degtorad
          measure(2*j-1,i) = sin(dum)
          measure(2*j,i) = cos(dum)
        end do
      end do
      deallocate (nc,nb,na,geo)
      !=========================================================================!
    case ('dihedral')
      mn = min(mm,2*ntaken) !>--- two descriptors per dihedral (sin/cos)
      allocate (measure(mn,mm),diedr(ndied))
      if (pr) then
        write (stdout,'(1x,a)') 'Using DIHEDRAL ANGLES as descriptors (sin/cos transformed):'
        do i = 1,ntaken
          write (stdout,'(1x,a,4i6)') 'Atoms: ',diedat(1:4,i)
        end do
        write (stdout,*)
      end if
      do i = 1,mm
        call calc_dieders(mols(i),ndied,diedat,diedr)
        do j = 1,min(ntaken,mn/2)
          measure(2*j-1,i) = sin(diedr(j))
          measure(2*j,i) = cos(diedr(j))
        end do
      end do
      if (allocated(diedat)) deallocate (diedat)
      if (allocated(diedr)) deallocate (diedr)
      allocate (pc(mn),pcvec(mm,mn))
      !=====================================================!
    end select
    !=====================================================!
    if (pr) then
      write (stdout,*)
      write (stdout,'(1x,a,i0,a,i0,a)') 'Performing SVD for ', &
&          mm,' structures and ',mn,' props'
    end if
    !>--- do the SVD  ! MM must not be smaller than MN !
    call SVD_to_PC(measure,mm,mn,pc,pcvec,.false.)
!=========================================================================================!
    call ctimer%stop(1)
  else
    write (stdout,*) 'There are not enough descriptors for a PCA!'
    write (stdout,*) 'Taking all structures as representative and writing ',clusterfile
    open (newunit=ich,file=clusterfile)
    do i = 1,nall
      call mols(i)%append(ich)
    end do
    close (ich)
    return
  end if
!========================================================================================!
  if (allocated(measure)) deallocate (measure)
  if (allocated(xyznew)) deallocate (xyznew)
  if (allocated(inc)) deallocate (inc)

  !>--- normalize PC eigenvalues
  pcsum = sum(pc)
  pc = pc/pcsum
  !>--- get the contributing principal components
  pcsum = 0.0d0
  npc = 0
  do i = 1,mn
    if (pc(i) < pcmin) exit
    pcsum = pcsum+pc(i)
    npc = npc+1
    if (pcsum .ge. pcthr) exit
  end do
  !npc = max(npc,2) !>-- at least 2 principal components should be used
  npc = min(npc,pccap)
  pcsum = 0.0d0
  do i = 1,npc
    pcsum = pcsum+pc(i)
  end do

  if (pr) then
    i = min(100,MM)
    k = min(npc,6)
    write (stdout,*)
    call smallhead('EIGENVECTORS AND NORMALIZED EIGENVALUES OF SVD ANALYIS')
    call PRMAT(6,pcvec,i,k,'Eigenvectors of principal components')
    write (stdout,'(1x,a,i0,a)') 'NOTE: eigenvectors are only shown for the first ',i,' structures'
    write (stdout,'(1x,a,i0,a)') '      and the first ',k,' contributing principal components.'
    write (stdout,*)

    write (stdout,*) mn,'principal component eigenvalues (normalized)'
    write (stdout,*) pc
    !call PRMAT(6,pc,mn,1,'Principal components (eigenvalues)')
    write (stdout,*)
    write (stdout,'(1x,a,i0,a,f6.2,a)') 'The first ',npc,' components account for a total of ',100.d0*pcsum,'% of the'
    write (stdout,'(1x,a)') 'ensembles unique structural features and are used for the clustering'
  end if

  !>--- use some less memory and rearrange the eigenvectors (COLUMNS AND ROWS ARE SWAPPED)
  !>    untaken PCs are not considered further
  allocate (pcdum(npc,mm))
  do i = 1,mm
    pcdum(1:npc,i) = pcvec(i,1:npc)
  end do
  call move_alloc(pcdum,pcvec)  !>THIS CHANGES THE SHAPE OF pcvec (COLUMNS AND ROWS ARE SWAPPED)

!=========================================================!
! do the Clustering
!=========================================================!

  if (pr) then
    write (stdout,*)
    call smallhead('CLUSTERING ANALYSIS OF PRINCIPAL COMPONENTS')
  end if

  allocate (member(mm),source=0)

  !>--- get Euclidean distances (packed matrix) between all structures
  !> ndist = (mm*(mm+1))/2 ! overflows for large ensembles
  ndist = mm
  ndist = ndist*(mm+1)
  ndist = ndist/2
  allocate (dist(ndist),source=0.0_sp)
  allocate (p(npc),q(npc),source=0.0_ap)

  do i = 1,mm
    p(1:npc) = pcvec(1:npc,i)
!$OMP PARALLEL PRIVATE ( j, klong, q, dum ) &
!$OMP SHARED ( i, dist, npc, p, pcvec )
!$OMP DO
    do j = 1,i
      q(1:npc) = pcvec(1:npc,j)
      dum = eucdist(npc,p,q)
      klong = lina(i,j)
      dist(klong) = real(dum,sp)
    end do
!$OMP END DO
!$OMP END PARALLEL
  end do

  !>-- NOTE
  !>-- different clustering algorithms exist, but what is common
  !>-- among them is, that no optimal number of clusters is known
  !>-- at the beginning. The lower bound for the number of
  !>-- clusters is the number of investigated PCs, the upper
  !>-- bound is the number of structures
  if (pr) then
    select case (clusteralgo)
    case ('means','kmeans')
      write (stdout,'(1x,a)') 'Using a MEANS cluster algorithm.'
    end select
    write (stdout,'(1x,a)') 'For a good review of cluster algorithms see'
    write (stdout,'(1x,a)') 'JCTC, 2007, 3, 2312 (doi.org/10.1021/ct700119m)'
    write (stdout,*)
    write (stdout,'(1x,a)') 'DBI = Davies-Bouldin index'
    write (stdout,'(1x,a)') 'pSF = pseudo F-statistic'
    write (stdout,'(1x,a)') 'SSR/SST = ratio of explained and unexplained variation'
    write (stdout,*)
    write (stdout,'(1x,a8,4x,a14,4x,a14,4x,a14)') 'Nclust','DBI','pSF','SSR/SST'
    write (stdout,'(1x,a8,4x,a14,4x,a14,4x,a14)') '------','-------------','-------------','-------------'
  end if

!-----------------------------------------------------------------------!
! Cluster evaluation settings
!-----------------------------------------------------------------------!
  if (env%maxcluster == 0) then
    !nclustmax=100  !some random default value
    call clustleveval(env,nclustmax,csthr,SSRSSTthr,pcthr) ! defaults
    nclustmax = min(mm,nclustmax)
    !SSRSSTthr=0.90  !exit if this value is reached  for SSR/SST
  else
    nclustmax = max(2,env%maxcluster)  !no less than 2 clusters
    nclustmax = min(mm,env%maxcluster) !there cannot be more clustes than structures.
  end if
  if (env%nclust == 0) then
    nclustmin = 1
  else
    !>-- predefined number of clusters
    nclust = min(mm,env%nclust)
    nclustmin = nclust
    nclustmax = nclust
  end if

  allocate (statistics(3,nclustmax),source=0.0d0)
  allocate (clust_sizes(nclustmax),source=0)
  CLUSTERSIZES: do nclustiter = nclustmin,nclustmax

    !>-- regular case: test continuous cluster sizes
    nclust = nclustiter
    !>-- special case: test cluster sizes incrementally (good for large ensembles)
    if (env%clustlev >= 10) then
      dum = float(mm)/float(nclustmax)
      dum2 = dum*float(nclustiter)
      nclust = nint(dum2)
    end if
    clust_sizes(nclustiter) = nclust

    allocate (centroid(npc,nclust),source=0.0_ap)
    centroid = 0.0_ap

    select case (clusteralgo)
    case ('means','kmeans')
      call ctimer%start(2,'k-Means clustering')
      call kmeans(nclust,npc,mm,centroid,pcvec,ndist,dist,member)
      call ctimer%stop(2)
    end select

    call ctimer%start(3,'statistics')
    call cluststat(nclust,npc,mm,centroid,pcvec,member,DBI,pSF,SSRSST)
    if (pr) then
      write (stdout,'(1x,i8,4x,f14.6,4x,f14.6,4x,f14.6)') nclust,DBI,pSF,SSRSST
    end if
    call ctimer%stop(3)
    deallocate (centroid)

    statistics(1,nclustiter) = DBI
    statistics(2,nclustiter) = pSF
    statistics(3,nclustiter) = SSRSST

    if (nclust == env%nclust) exit
    if (SSRSST > SSRSSTthr) exit
  end do CLUSTERSIZES
  if (allocated(centroid)) deallocate (centroid)

  write (stdout,*)
  if (env%nclust == 0) then
    if (pr) then
      write (stdout,'(1x,a,i0,a)') 'Ensemble checked up to a partitioning into ',nclust,' clusters.'
      write (stdout,'(1x,a)') 'Local MINIMA of the DBI indicate adequate cluster counts.'
      write (stdout,'(1x,a)') 'Local MAXIMA of the pSF indicate adequate cluster counts.'
      write (stdout,'(1x,a)') 'Higher SSR/SST vaules indicate more distinct clusters.'
      write (stdout,'(1x,a)') 'Analyzing statistical values ...'
    end if
    k = min(nclustiter,nclustmax) !> last completed iteration index
    allocate (extrema(2,k))
    call ctimer%start(3,'statistics')
    call statanal(k,nclustmax,statistics,extrema,pr,clust_sizes)
    if (pr) call statwarning(fname)
    !>-- determine a suggested cluster size (smallest suggested cluster with good SSR/SST)
    do i = 2,k
      if ((extrema(1,i).or.extrema(2,i)).and.(statistics(3,i) > csthr)) then
        nclust = clust_sizes(i)
        exit
      end if
    end do
    call ctimer%stop(3)
    deallocate (extrema)
    if (pr) then
      write (stdout,*)
      write (stdout,'(1x,a,f4.2,a,i0)') 'Suggested (SSR/SST >',csthr,') cluster count: ',nclust
    end if
    !>-- calculate the determined partition into clusters again for final file
    allocate (centroid(npc,nclust),source=0.0_ap)
    select case (clusteralgo)
    case ('means','kmeans')
      call ctimer%start(2,'k-Means clustering')
      call kmeans(nclust,npc,mm,centroid,pcvec,ndist,dist,member)
      call ctimer%stop(2)
    end select
    deallocate (centroid)
  else
    if (pr) then
      write (stdout,'(1x,a,i0,a)') 'Ensemble partitioning into ',nclust,' clsuters.'
    end if
  end if
  deallocate (statistics,clust_sizes)

  deallocate (q,p,dist)
  !>-- finally, assign a representative structure to each group (based on lowest energy)
  !>-- and write the new ensemble file
  call PCA_grpwrite(nclust,npc,mm,pcvec,member)

  !ncb=maxval(member,1) !--total number of cluster
  ncb = nclust
  ancb = ncb           !>--actual number of clusters

  if (ancb .le. 1) return

  if (pr) then
    write (stdout,*)
    write (stdout,'(1x,a)') 'Representative structures'
    write (stdout,'(1x,a6,1x,a6,3x,a6,1x,a16,1x,a16)') 'Nr.','conf.','clust.','Etot/Eh','Erel/ kcal/mol'
    write (stdout,'(1x,a6,1x,a6,3x,a6,1x,a16,1x,a16)') '---','-----','------','-------','--------------'
  end if
  allocate (eclust(ncb),source=0.0d0)
  allocate (clustbest(ncb),ind(ncb),source=0)
  !>--- initialize eclust and clustbest
  iiincb: do i = 1,ncb
    do j = 1,mm
      if (member(j) == i) then
        eclust(i) = mols(j)%energy
        clustbest(i) = j
        cycle iiincb
      end if
    end do
  end do iiincb
  !>--- then, look for the lowest energy
  do i = 1,ncb
    ind(i) = i
    c = 0
    do j = 1,mm
      if (member(j) == i) then
        c = c+1
        if (mols(j)%energy < eclust(i)) then
          eclust(i) = mols(j)%energy
          clustbest(i) = j
        end if
      end if
    end do
    !>-- if there are clusters without structures
    if (c == 0) then
      clustbest(i) = -1
      ancb = ancb-1
    end if
  end do
  call qsort(eclust,1,ncb,ind)
  emin = minval(eclust,1)
  open (newunit=ich,file=clusterfile)
  do i = 1,ncb
    k = clustbest(ind(i))
    if (k > 0) then
      dum = mols(k)%energy
      call mols(k)%append(ich)
      if (pr) then
        erel = (dum-emin)*autokcal
        write (stdout,'(1x,i6,1x,i6,3x,i6,1x,f16.8,1x,f16.4)') i,k,member(k),dum,erel
      end if
    end if
  end do
  close (ich)
  if (pr) then
    write (stdout,'(/,1x,a)') '(The "clust." column refers to the cluster "ID")'
    write (stdout,*)
    write (stdout,'(1x,a,a,a,i0,a)') 'File ',clusterfile,' written with ',ancb,' representative structures.'
    if (ancb < ncb) then
      write (stdout,'(1x,a,i0,a)') '(',ncb-ancb,' clusters discarded due to cluster merge)'
    end if
  end if
  if (allocated(mols)) deallocate (mols)

  if (pr) then
    write (stdout,*)
    !write(*,'(1x,a)') 'Timings:'
    !call eval_sub_timer(ctimer)
    call ctimer%write(stdout,'PCA/k-Means clustering')
  end if
  call ctimer%clear()
  return
end subroutine CCEGEN

!=======================================================================================!
! set clustering level defaults
!=======================================================================================!
subroutine clustleveval(env,maxclust,csthr,SSRSSTthr,pcthr)
  !*********************************************************
  !* Set clustering level defaults for maxclust, csthr,    *
  !* SSRSSTthr, and pcthr based on the clustering level.   *
  !*********************************************************
  use crest_parameters,idp => dp
  use crest_data
  implicit none
  type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
  integer :: clev
  integer :: maxclust
  real(wp) :: csthr
  real(wp) :: SSRSSTthr
  real(wp) :: pcthr

  SSRSSTthr = 0.90  !exit if this value is reached  for SSR/SST

  clev = env%clustlev
  if (env%clustlev >= 10) then  !for incremental modes
    clev = env%clustlev-10
  end if

  select case (clev)
  case (-1) !-- loose
    maxclust = 25
    csthr = 0.80d0
    pcthr = 0.80d0
  case (1)  !-- tight
    maxclust = 400
    if (env%clustlev >= 10) maxclust = 50
    csthr = 0.85d0
    pcthr = 0.90d0
  case (2)  !-- vtight
    maxclust = 400
    if (env%clustlev >= 10) maxclust = 100
    csthr = 0.9d0
    pcthr = 0.95d0
    SSRSSTthr = 0.92d0
  case default !-- normal
    maxclust = 100
    if (env%clustlev >= 10) maxclust = 25
    csthr = 0.80d0
    pcthr = 0.85d0
  end select

  return
end subroutine clustleveval

!=======================================================================================!
! write a file with the 2 most contributing principal components of each structure
! and the cluster to which the structure belongs
!=======================================================================================!
subroutine PCA_grpwrite(nclust,npc,mm,pcvec,member)
  use crest_parameters,idp => dp
  implicit none
  integer,intent(in) :: nclust  ! number of required centroids
  integer,intent(in) :: npc,mm
  real(wp),intent(in) :: pcvec(npc,mm)
  integer,intent(in) :: member(mm)  ! membership for each structure
  integer :: ich,i
  open (newunit=ich,file='cluster.order')
  write (ich,'(4x,i0,4x,i0,4x,i0)') mm,nclust,npc
  if (npc > 1) then
    do i = 1,mm
      write (ich,'(i8,1x,f16.8,1x,f16.8,1x,i8)') i,pcvec(1,i),pcvec(2,i),member(i)
    end do
  else
    do i = 1,mm
      write (ich,'(i8,1x,f16.8,1x,i8)') i,pcvec(1,i),member(i)
    end do
  end if
  close (ich)
  return
end subroutine PCA_grpwrite

!=======================================================================================!
! Exclude light - exclude H atoms in the inc array
!=======================================================================================!
subroutine excludeLight(zmol,inc)
  use crest_parameters,idp => dp
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: inc(zmol%nat)
  integer :: i
  do i = 1,zmol%nat
    if (zmol%at(i) == 1) then
      inc(i) = 0
    end if
  end do
  return
end subroutine excludeLight

!=======================================================================================!
! Exclude Specified Atoms
!=======================================================================================!
subroutine excludeSelected(zmol,inc,atlist)
  use crest_parameters,idp => dp
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: inc(zmol%nat)
  character(len=*) :: atlist   !a string containing atom numbers, needs parsing
  integer :: i,ncon
  integer,allocatable :: inc2(:)
  allocate (inc2(zmol%nat),source=0)
  call parse_atlist_new(atlist,ncon,zmol%nat,zmol%at,inc2)
  do i = 1,zmol%nat
    if (inc2(i) == 1) inc(i) = 0
  end do
  deallocate (inc2)
  return
end subroutine excludeSelected

!=======================================================================================!
! Perform a single value decomposition (SVD) and get the principal components
!
! X = U*sig*V^(T)
!
! The eigenvalues saved in sig are the principal components
! The SVD only works if M >= N
!
!=======================================================================================!
subroutine svd_to_pc(measure,m,n,sig,U,pr)
  use crest_parameters,idp => dp
  implicit none
  integer :: n,m
  real(wp) :: measure(n,m)
  real(wp) :: sig(n)
  real(wp) :: U(m,n)
  integer :: i,j,info,lwork
  real(wp),allocatable :: mean(:),tmp(:)
  integer,allocatable :: ind(:)
  real(wp),allocatable :: X(:,:),V(:,:),work(:)
  integer,allocatable :: iwork(:)
  logical :: pr
  if (pr) then
    write (stdout,*) m,' mesaurements'
    write (stdout,*) n,' props'
  end if
  allocate (mean(n),ind(m),tmp(m))
  lwork = max(2*M+N,6*N+2*N*N)
  allocate (X(m,n),V(n,n),iwork(m+3*n),work(lwork))
  mean = 0.0d0
  do i = 1,m
    do j = 1,n
      mean(j) = mean(j)+measure(j,i)
    end do
  end do
  mean = mean/float(m)
  if (pr) write (stdout,*) mean
  do i = 1,m
    do j = 1,n
      X(i,j) = (mean(j)-measure(j,i))
    end do
  end do
  if (pr) then
    call PRMAT(6,X,m,n,'X')
  end if
!>--- LAPACKs' DGEJSV
  call DGEJSV('C','U','V','N','N','N', &
 &              m,n,X,m,sig,U,m,V,n,                &
 &              WORK,LWORK,IWORK,INFO)
  if (pr) then
    write (stdout,*) info
    write (stdout,*) sig
    call PRMAT(6,U**2,M,N,'U')
    call PRMAT(6,V,N,N,'V')
  end if
  deallocate (work,iwork,V,X,tmp,ind,mean)
  return
end subroutine svd_to_pc

!======================================================================!
! calculate the Euclidian distance between two points p and q
!======================================================================!
function eucdist(ndim,p,q) result(dist)
  use crest_parameters,idp => dp
  implicit none
  real(ap) :: dist
  integer :: ndim
  real(ap) :: p(ndim)
  real(ap) :: q(ndim)
  integer :: i
  dist = 0.0d0
  do i = 1,ndim
    dist = dist+(q(i)-p(i))**2
  end do
  dist = sqrt(dist)
  return
end function eucdist

!======================================================================!
! K-means clustering algorithm
!
! determine a position of cluster centroids iteratively for a given
! number of centroids.
!
!======================================================================!
subroutine kmeans(nclust,npc,mm,centroid,pcvec,ndist,dist,member)
  use crest_parameters,idp => dp
  implicit none
  integer,intent(in) :: nclust  ! number of required centroids
  integer,intent(in) :: npc,mm
  real(wp),intent(in) :: pcvec(npc,mm)
  integer(idp),intent(in) :: ndist
  real(sp),intent(in)   :: dist(ndist)
  integer,intent(inout) :: member(mm)  ! membership for each structure
  real(ap),intent(inout):: centroid(npc,nclust)
  integer,allocatable :: refmember(:)
  integer :: iter
  integer,parameter :: maxiter = 300

  if (nclust .le. 1) return !no singular clusters!

  allocate (refmember(mm),source=0)

  !>-- determine seeds for the centroids (i.e., initial positions)
  call kmeans_seeds(nclust,npc,mm,centroid,pcvec,ndist,dist)

  do iter = 1,maxiter
    !>-- determine cluster membership for all structures
    !>   (by shortest Euc. distance to the respective centroid)
    member = 0 !reset
    call kmeans_assign(nclust,npc,mm,centroid,pcvec,member)

    !>-- check if memberships changed w.r.t. previous memberships
    if (all(member == refmember)) then
      exit
    else
      refmember = member
    end if
    !>-- update centroids if necessary
    call kmeans_recenter(nclust,npc,mm,centroid,pcvec,member)
  end do

  deallocate (refmember)
  return
end subroutine kmeans
!===================================================================!
! determine cluster seeds for the K-means algo
!===================================================================!
subroutine kmeans_seeds(nclust,npc,mm,centroid,pcvec,ndist,dist)
  use crest_parameters,idp => dp
  use utilities
  implicit none
  integer :: nclust,npc,mm
  real(ap) :: centroid(npc,nclust)
  real(wp) :: pcvec(npc,mm)
  integer(idp) :: ndist
  real(sp) :: dist(ndist)
  real(sp) :: ddum
  integer(idp) :: k,kiter
  integer :: i,j,l,c
  real(wp) :: distsum,maxdistsum
  real(ap) :: eucdist !this is a function
  real(ap),allocatable :: p(:),q(:)
  integer,allocatable :: taken(:)

  !>--- first two centroids are located at the most apart points
  !>    in the PC space
  ddum = 0.0_sp
  do kiter = 1,ndist
    if (dist(kiter) > ddum) then
      ddum = dist(kiter)
      k = kiter
    end if
  end do
  call revlin(k,j,i) !> reverse of the lin function, get indices i and j

  centroid(1:npc,1) = pcvec(1:npc,i)
  centroid(1:npc,2) = pcvec(1:npc,j)

  if (nclust .le. 2) return !>-- if we only need two centroids, return

  !>-- If more centroids are needed, search for one point which has the maximal sum of
  !>-- the distances between the already determined centroids and itself.
  allocate (p(npc),q(npc),taken(nclust))
  taken = 0
  taken(1) = i
  taken(2) = j
  do i = 3,nclust
    maxdistsum = 0.0d0
    c = 0
!$OMP PARALLEL PRIVATE ( l, q, p, j, distsum ) &
!$OMP SHARED ( i, centroid, npc, mm, maxdistsum, pcvec, c, taken )
!$OMP DO
    do j = 1,mm
      distsum = 0.0d0
      p(1:npc) = pcvec(1:npc,j)
      do l = 1,i-1
        q(1:npc) = centroid(1:npc,l)
        distsum = distsum+eucdist(npc,p,q)
      end do
      !$OMP CRITICAL
      if (distsum .gt. maxdistsum) then
        if (.not.any(taken == j)) then
          maxdistsum = distsum
          c = j
          taken(i) = c
        end if
      end if
      !$OMP END CRITICAL
    end do
!$OMP END DO
!$OMP END PARALLEL
    if (c == 0) then
      exit
    else
      centroid(1:npc,i) = pcvec(1:npc,c)
    end if
  end do

  deallocate (taken,q,p)

  return
end subroutine kmeans_seeds
!===================================================================!
! assign structures as members to a centroid
!===================================================================!
subroutine kmeans_assign(nclust,npc,mm,centroid,pcvec,member)
  use crest_parameters,idp => dp
  implicit none
  integer :: nclust,npc,mm
  real(ap) :: centroid(npc,nclust)
  real(wp) :: pcvec(npc,mm)
  integer :: member(mm)
  real(ap) :: eucdist  !this is a function
  integer :: i,j,c
  real(ap),allocatable :: centdist(:)
  real(ap),allocatable :: p(:),q(:)

  allocate (centdist(nclust),source=0.0_ap)
  allocate (p(npc),q(npc))
!$OMP PARALLEL PRIVATE ( i, j, p, q, c, centdist ) &
!$OMP SHARED ( mm, nclust, member, centroid, npc, pcvec )
!$OMP DO
  do i = 1,mm
    p(1:npc) = pcvec(1:npc,i)
    do j = 1,nclust
      q(1:npc) = centroid(1:npc,j)
      centdist(j) = eucdist(npc,p,q)
    end do
    c = minloc(centdist,1)
    member(i) = c
  end do
!$OMP END DO
!$OMP END PARALLEL
  deallocate (q,p)
  deallocate (centdist)
  return
end subroutine kmeans_assign

!===================================================================!
! re-center centroids for given (sorted) structures
!===================================================================!
subroutine kmeans_recenter(nclust,npc,mm,centroid,pcvec,member)
  use crest_parameters,idp => dp
  implicit none
  integer :: nclust,npc,mm
  real(ap) :: centroid(npc,nclust)
  real(wp) :: pcvec(npc,mm)
  integer :: member(mm)
  integer :: i,j,c
  real(ap),allocatable :: p(:),q(:)

  allocate (p(npc),q(npc))
  do i = 1,nclust
    c = 0
    p = 0.0d0
    do j = 1,mm
      if (member(j) == i) then
        c = c+1
        p(1:npc) = p(1:npc)+pcvec(1:npc,j)
      end if
    end do
    if (c > 0) then
      p = p/float(c)
      centroid(1:npc,i) = p(1:npc)
    else
      p = 999.9d0
    end if
  end do
  deallocate (q,p)

  return
end subroutine kmeans_recenter

!======================================================================!
! calculate statistical values for the given cluster size
! Values to compute:
! DBI - the Davies-Bouldin index
! pSF - the "pseudo-F statistic"
! SSR/SST ratio
!======================================================================!
subroutine cluststat(nclust,npc,mm,centroid,pcvec,member,DBI,pSF,SSRSST)
  use crest_parameters,idp => dp
  implicit none
  integer,intent(in) :: nclust  ! number of required centroids
  integer,intent(in) :: npc,mm
  real(wp),intent(in) :: pcvec(npc,mm)
  integer,intent(in) :: member(mm)  ! membership for each structure
  real(ap),intent(in):: centroid(npc,nclust)
  real(wp),intent(out) :: DBI,pSF,SSRSST
  real(wp) :: SSE,SSR,SST
  real(ap),allocatable :: p(:),q(:)
  real(wp),allocatable :: compact(:)
  real(wp),allocatable :: DBmat(:,:)
  real(ap) :: eucdist !this is a function
  real(wp) :: d,Rij,maxDB,weight
  integer :: i,c,k,c2

  DBI = 0.0d0
  pSF = 0.0d0
  SSRSST = 0.0d0

  if (nclust < 2) return

  allocate (p(npc),q(npc))

  !>-- Sum of squares error
  SSE = 0.0d0
  do c = 1,nclust
    p(1:npc) = centroid(1:npc,c)
    do i = 1,mm
      if (member(i) == c) then
        q(1:npc) = pcvec(1:npc,i)
        d = eucdist(npc,p,q)
        SSE = SSE+d**2
      end if
    end do
  end do
  SSE = SSE

  !>-- Total sum of squares
  SST = 0.0d0
  p = 0.0d0
  do c = 1,nclust
    weight = real(count(member(:) == c,1),wp)/real(mm,wp)
    p(1:npc) = p(1:npc)+centroid(1:npc,c)*weight
  end do
  !p = p/float(nclust)
  do i = 1,mm
    q(1:npc) = pcvec(1:npc,i)
    d = eucdist(npc,p,q)
    SST = SST+d**2
  end do
  SST = SST

  !>-- Sum of squares regression
  SSR = SST-SSE

  !>-- SSR/SST ratio
  SSRSST = SSR/SST

  !>-- pseudo-F statistic
  if (nclust > 1) then
    pSF = (SSR/(float(nclust)-1.0d0))
    if (mm == nclust) then
      pSF = 0.0d0
    else
      pSF = pSF/(SSE/(float(mm)-float(nclust)))
    end if
  else
    pSF = 0.0d0
  end if

  !>-- Davies-Bouldin index (DBI)
  allocate (compact(nclust),source=0.0d0) !cluster compactness
  do c = 1,nclust
    p(1:npc) = centroid(1:npc,c)
    k = 0
    do i = 1,mm
      if (member(i) == c) then
        k = k+1
        q(1:npc) = pcvec(1:npc,i)
        d = eucdist(npc,p,q)
        compact(c) = compact(c)+d
      end if
    end do
    if (k > 0) then
      compact(c) = compact(c)/float(k)
    else
      compact(c) = 0
    end if
  end do
  allocate (DBmat(nclust,nclust),source=0.0d0)
  do c = 1,nclust
    p(1:npc) = centroid(1:npc,c)
    do c2 = 1,nclust
      if (c2 == c) cycle
      q(1:npc) = centroid(1:npc,c2)
      d = eucdist(npc,p,q)
      Rij = (compact(c)+compact(c2))/d
      DBmat(c,c2) = Rij
    end do
  end do
  do c = 1,nclust
    maxDB = maxval(DBmat(:,c),1)
    DBI = DBI+maxDB
  end do
  DBI = DBI/float(nclust)
  deallocate (DBmat)
  deallocate (compact)

  deallocate (q,p)
  return
end subroutine cluststat

!==============================================================!
! analyze the statistical values DBI and pSF to get the
! respective extrema
!==============================================================!
subroutine statanal(n,nmax,statistics,extrema,pr,clust_sizes)
  use crest_parameters
  implicit none
  integer :: n,nmax
  real(wp) :: statistics(3,nmax)
  logical,intent(inout) :: extrema(2,n)
  logical :: pr
  integer,intent(in),optional :: clust_sizes(n)
  real(wp) :: last,next,current
  integer :: i,csize

  extrema = .false.
!>--- identify local extrema of the DBI and pSF
  do i = 2,n-1
    !>-- DBI
    last = statistics(1,i-1)
    next = statistics(1,i+1)
    current = statistics(1,i)
    if ((current < last).and.(current < next)) then
      extrema(1,i) = .true.
    end if
    !>-- pSF
    last = statistics(2,i-1)
    next = statistics(2,i+1)
    current = statistics(2,i)
    if ((current > last).and.(current > next)) then
      extrema(2,i) = .true.
    end if
  end do
  !>--- boundary check: last cluster count (one-sided comparison)
  if (n >= 2) then
    if (statistics(1,n) < statistics(1,n-1)) extrema(1,n) = .true.
    if (statistics(2,n) > statistics(2,n-1)) extrema(2,n) = .true.
  end if

  if (pr) then
    write (stdout,*)
    write (stdout,'(1x,a,/)') 'Suggestions for cluster sizes:'
    do i = 1,n
      if (extrema(1,i).or.extrema(2,i)) then
        csize = i
        if (present(clust_sizes)) csize = clust_sizes(i)
        if (extrema(1,i).and.extrema(2,i)) then
          write (stdout,'(1x,i8,''*'',3x,a,f8.4)') csize,'SSR/SST',statistics(3,i)
        else
          write (stdout,'(1x,i8,4x,a,f8.4)') csize,'SSR/SST',statistics(3,i)
        end if
      end if
    end do
    write (stdout,'(/,1x,a)') 'Cluster counts marked with a star (*) are reasonable'
    write (stdout,'(1x,a)') 'suggestions according to BOTH the DBI and pSF.'
  end if

  return
end subroutine statanal

!==============================================================!
! print a warning regarding the nature of the cluster partitioning
!==============================================================!
subroutine statwarning(fname)
  use crest_parameters
  implicit none
  character(len=*) :: fname
  write (stdout,*)
  write (stdout,'(1x,a)') '!---------------------------- NOTE ----------------------------!'
  write (stdout,'(2x,a)') 'The partitioning of data (the ensemble) into clusters'
  write (stdout,'(2x,a)') 'of similar characteristics (structures) is ARBITRARY'
  write (stdout,'(2x,a)') 'and depends on many criteria (e.g. choice of PCs).'
  write (stdout,'(2x,a)') 'The selected cluster count is the smallest reasonable'
  write (stdout,'(2x,a)') 'number of clusters that can be formed according to'
  write (stdout,'(2x,a)') 'the DBI and pSF values for the given data.'
  write (stdout,*)
  write (stdout,'(2x,a)') 'If other cluster sizes are desired, rerun CREST with'
  write (stdout,'(2x,3a)') '"crest --sort ',trim(fname),' --cluster <number of clusters>"'
  write (stdout,*)
  write (stdout,'(2x,a)') 'Other default evaluation settings can be chosen with the'
  write (stdout,'(2x,a)') 'keywords "loose","normal", and "tight" as <level> via'
  write (stdout,'(2x,3a)') '"crest --sort ',trim(fname),' --cluster <level>"'
  write (stdout,'(1x,a)') '!--------------------------------------------------------------!'
end subroutine statwarning

!====================================================================!
subroutine getdiederatoms(zmol,nat,inc,nb,diedat,ndied)
  use crest_parameters,idp => dp
  use zdata
  use strucrd
  implicit none
  type(zmolecule) :: zmol
  integer :: nat
  integer :: inc(nat) !contains 1 (=include) or 0 (=ignore)
  integer :: nb
  integer :: diedat(4,nb)
  integer,intent(out) :: ndied
  integer :: a,b,c,d
  integer :: i,j,k

  ndied = 0
  do i = 1,nb
    a = zmol%bondpairs(1,i)
    b = zmol%bondpairs(2,i)
    if (inc(a) == 0) cycle          !ignored by user?
    if (inc(b) == 0) cycle          !ignored by user?
    if (zmol%zat(a)%nei == 1) cycle !terminal atom?
    if (zmol%zat(b)%nei == 1) cycle !terminal atom?
    if (zmol%methyl(a)) cycle !methyl C?
    if (zmol%methyl(b)) cycle !methyl C?
    !>-- passed all checks, so let's get atoms
    !>-- a neighbour for a
    do j = 1,zmol%zat(a)%nei
      c = zmol%zat(a)%ngh(j)
      if (c == b) then
        cycle
      else
        exit
      end if
    end do
    !>-- a neighbour for b
    do k = 1,zmol%zat(b)%nei
      d = zmol%zat(b)%ngh(k)
      if (d == a) then
        cycle
      else
        exit
      end if
    end do
    ndied = ndied+1
    !>the bond is between a and b
    !>c is a neighbour of a, d is a neighbour of b
    diedat(2,ndied) = a
    diedat(3,ndied) = b
    diedat(1,ndied) = c
    diedat(4,ndied) = d
  end do

  return
end subroutine getdiederatoms

subroutine calc_dieders(mol,ndied,diedat,diedr)
  !*****************************************************
  !* Calculate dihedral angles for selected atom       *
  !* quartets. Results are in radians (-pi, pi).       *
  !*****************************************************
  use crest_parameters,idp => dp
  use strucrd
  implicit none
  type(coord),intent(in) :: mol
  integer,intent(in) :: ndied
  integer,intent(in) :: diedat(4,ndied)
  real(wp),intent(out) :: diedr(ndied)
  integer :: i

  diedr = 0.0_wp
  do i = 1,ndied
    !> diedat: (1)=neighbour of a, (2)=a, (3)=b, (4)=neighbour of b
    diedr(i) = mol%dihedral(diedat(1,i),diedat(2,i), &
    &                       diedat(3,i),diedat(4,i))
  end do

  return
end subroutine calc_dieders

