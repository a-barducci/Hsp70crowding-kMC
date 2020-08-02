program kmc
  implicit none
  integer i,j,k
  integer it
  integer nstot,nbtot,nbfree,rpen
  integer, allocatable:: occ(:)
  real(8), allocatable:: koff(:),offrate(:),onrate(:)
  real(8):: kon,koff_start,pen,npen
  real(8):: t,dt,tmax,twrite
  real(8):: br,ur,rr,comp
  real(8):: avocc,avkd,avfree,avocc2
  real(8):: sumocc,sumfree,sumocc2
  real(8):: trate,tonrate,toffrate
  real(8):: ran1,ran2,a
  logical:: debug,dtreal

  open(unit=10,file='input')
  open(unit=77,file='extended.out')

  debug=.false.  ! extra dump 
  dtreal=.false.  ! if dtreal=.false. the time update is equal to the inverse of the sum of all the rates at time t 

  !! parameters
  read (10,*)
  read (10,*) nbtot ! total number of binders
  read (10,*) 
  read (10,*) nstot ! total number of receptor sites
  read (10,*)
  read (10,*) tmax  ! simulation time
  read (10,*)
  read (10,*) twrite ! stride for output
  read (10,*)
  read (10,*) koff_start ! basal value of koff
  read (10,*)
  read (10,*) kon ! kon
  read (10,*)
  read (10,*) npen ! energetic penalty for crowding
  read (10,*) 
  read (10,*) rpen ! range of the penalty (number of neighboring sites affected by binding at site j)
  
  write (*,*) 'Nbinders ',nbtot
  write (*,*) 'Nsites ',nstot 
  write (*,*) 'Tmax ', tmax
  write (*,*) 'Koff ', koff_start
  write (*,*) 'Kon ', kon
  write (*,*) 'penalty kBT ', npen 
  write (*,*) 'penalty range ',rpen


  
  call random_seed()
  allocate(occ(nstot))
  allocate(koff(nstot),offrate(nstot),onrate(nstot))
  t=0.0
  it=1

  ! calculate penalty factor on rates, assuming kbT=0.593 kcal/mol
  pen=exp(npen/0.593)
  
!  write(*,*) pen

  occ(:)=0
  koff(:)=koff_start
  nbfree=nbtot
  sumocc=0
  sumfree=0
  sumocc2=0

  do !! main KMC cycle
     
     if (t.ge.tmax) exit

     !! computing all the possible binding/unbinding events at time t

     !! cycling on receptor sites
     do j=1,nstot
        !! if site j is free, it can be bound with a rate corresponding to kon multiplied by the number of free binders (nbfree)
        if (occ(j)==0) then 
           onrate(j)=kon*nbfree 
           offrate(j)=0.0
        else
        !! if site j is occupied, an unbinding event can occur with koff(j)
           onrate(j)=0.0
           offrate(j)=koff(j)
        endif
     enddo
     
     tonrate=sum(onrate(:))    ! cumulative binding rate
     toffrate=sum(offrate(:))  ! cumulative unbinding rate
     
     write (77,*) t,nbfree,sum(occ(:)),nstot, tonrate, toffrate

     trate= tonrate+toffrate   ! overall cumulative rate

     if (dtreal) then

     !! calculate stochastic time update according to the standard kMC procedure 

        call random_number(ran2)
        dt=(1/trate)*log(1/(1.0-ran2))

     else

     !! use deterministic time update

        dt=(1/trate)

     endif
     
     !! update cumulative counters

     sumocc=sumocc+sum(occ(:))*dt
     sumocc2=sumocc2+sum(occ(:))**2*dt
     sumfree=sumfree+nbfree*dt
     
     !! pick random binding/unbinding event according to kmc algorithm

     call random_number(ran1)

     rr=(1.0-ran1)*trate
     
     if (rr.le.tonrate) then ! binding event
        
        ! which site?

        do j=1,nstot
           br=sum(onrate(1:j))
           if (rr.le.br) exit
        enddo

        if (occ(j)==1) write (*,*) 'error binding' ! sanity check: site j should be free
        
        ! update site j, which is now occupied
        
        occ(j)=1
        
        ! check for neighbouring occupied sites and apply crowding penalty by increasing the unbinding rate constants both for this binder and the neighbours 

        do k=j-rpen,j+rpen
           if ((k.ne.j).and.(k.ge.1).and.(k.le.nstot)) then
              if (occ(k)==1) then
                 koff(k)=koff(k)*pen
                 koff(j)=koff(j)*pen
              endif
           endif
        enddo
        
        if (debug) write(78,*) maxval(koff(:)),log(maxval(koff(:))/koff_start),pen**log(maxval(koff(:))/koff_start)
        
        ! update the number of free binders
        nbfree=nbfree-1

        if (debug)  write(*,*) 'binding',j,nbfree

     else ! unbinding event
        
        ! which site ?
        comp=(rr-tonrate)
        do j=1,nstot
           ur=sum(offrate(1:j))
           if (comp.le.ur) exit
        enddo

        if (occ(j)==0) write (*,*) 'error unbinding' ! sanity check: site j should be occupied
        
        ! update site j, which is now free

        occ(j)=0

        ! update number of free binders
        nbfree=nbfree+1
        
        ! remove crowding penalty
        koff(j)=koff_start
        
        do k=j-rpen,j+rpen
           if ((k.ne.j).and.(k.ge.1).and.(k.le.nstot)) then
              if (occ(k)==1) then
                 koff(k)=koff(k)/pen
              endif
           endif
        enddo
      endif
     
     !! update time
     t=t+dt
!
     !! calculate time averages
     
     avfree=sumfree/t
     avocc2=sumocc2/t
     avocc=sumocc/t
     
     !! check if it is time to dump
     if (t.ge.it*twrite) then
        it=it+1
!       write time, average number of free binders and bound sites, and estimates for Kd        
        write(*,*) t, avfree, avocc, avfree*(nstot-avocc)/avocc, (avfree*(nstot-avocc)+(avocc2-avocc**2))/avocc
     endif

  enddo
  write(*,*) 'final averages'  
  write(*,*) t, avfree, avocc, avfree*(nstot-avocc)/avocc, (avfree*(nstot-avocc)+(avocc2-avocc**2))/avocc 
  write (*,*) 'bound fraction'
  write (*,*) avocc/nbtot

end program kmc
      
