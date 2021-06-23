!**********************************************************************
!                               WriteMain
!**********************************************************************
subroutine WriteMain()
  use gblcom
  use utils
  implicit none
  integer :: lun
  integer ios, ispace(num_files)
  character*80 file_name*80, version_string*20
  !version_string_rd must be declared with same lenth in rpd

  version_string="main Version 2.1"
  file_name=trim(savedir)//'/main'
  if(MakeBackupCopy(trim(file_name))) continue
  open(newunit=lun, file=file_name, access="stream",action="write",&
       form="unformatted", iostat=ios)
  if(ios.ne.0) then
     call ceprint("#Error in routine wpd opening file "//trim(file_name)//"#")
     return
  endif
  ispace(:)=0               ! integer dummy array
  write(lun,iostat=ios) version_string
  write(lun,iostat=ios) ndet, rectat,  nevtat, kpa, npd, num_files,&
       tcalfact, trigger_mode, stab_ignore_single(1:ndet),&
       stab_timu(1:ndet), stab_timo(1:ndet), stab_limu(1:ndet),&
       stab_limo(1:ndet), stab_live_removed(1:ndet),&
       kam(1:ndet), koff_limit(1:ndet), kofi(1:ndet), koff(1:ndet),&
       kmaf(1:ndet), gauge(1:kpa,1:ndet), xyl(1:kpa,1:ndet),&
       xyh(1:kpa,1:ndet), xyl_complete(1:kpa,1:ndet),&
       xyh_complete(1:kpa,1:ndet), xyl_selected(1:kpa,1:ndet),&
       xyh_selected(1:kpa,1:ndet),& 
       ispace(1:num_files), flist(1:num_files), ptx,&
       ftx,& ! originally num_flags istead of num_flags+1 were saved 
       stab_usecs_file(1:ndet, 1:num_files),&
       stab_dead_secs_file(1:ndet, 1:num_files), kpr,&
       triggers(1:ndetm), recsum(1:num_files),&
       ! added in version 2.1
       t_file_begin(:), t_file_end(:), time_to_add(:),tsum

  close(unit=lun)
  if(ios.ne.0) then
     call ceprint("Error in routine wpd writing file "//trim(file_name)//"#")
  endif
  return
end subroutine WriteMain