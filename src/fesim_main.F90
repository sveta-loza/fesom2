!=============================================================================!
!
!                 Finite Volume Sea-ice Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    

program main
  use fesim_module

  integer nsteps

  call fesim_init(nsteps)
  call fesim_runloop(nsteps)
  call fesim_finalize

end program main
