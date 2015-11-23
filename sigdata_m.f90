module sigdata_m

private
public ovar,zs,zsi_m,lsm_m,sigdataalloc,sigdatadealloc

real, dimension(:), allocatable, save :: ovar,zs
real, dimension(:), allocatable, save :: zsi_m,lsm_m

contains

subroutine sigdataalloc(il)

implicit none

integer, intent(in) :: il
integer jl,ifull

jl=6*il
ifull=il*jl

allocate(ovar(ifull),zs(ifull))
allocate(zsi_m(ifull),lsm_m(ifull))

return
end subroutine sigdataalloc

subroutine sigdatadealloc

implicit none

deallocate(ovar,zs)
deallocate(zsi_m,lsm_m)

return
end subroutine sigdatadealloc

end module sigdata_m
